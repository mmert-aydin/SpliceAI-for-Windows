"""The reference score database in the window (spec 5.7, 5.8): its status,
the guided build (Route 2), import, and the cache of genes computed live.

RefdbManager keeps track of the two databases (hg19, hg38) and runs builds,
one at a time, each as a process of its own (spliceai_pipeline.refdb_build,
or the exe's --build-refdb) steered through files in the database folder
(refdb.py: progress.json, control.txt, priority.txt). The builds still to do
are kept in the settings, so one interrupted by closing the program, a crash
or a restart continues at the next start. RefdbDialog shows all of it,
RefdbWizard is the guided build, and MissingDatabasePanel is what the
Sequence-to-cDNA tab shows while there is no database yet.
"""
import os
import shutil
import subprocess
import sys
import threading

from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox, QDialog, QFileDialog, QFrame, QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton,
    QRadioButton, QSpinBox, QVBoxLayout, QWidget, QWizard, QWizardPage,
)

from spliceai_pipeline import refdb, transcripts
from spliceai_pipeline.app_paths import app_root_dir
from spliceai_pipeline.cli import format_elapsed
from spliceai_pipeline.mane import find_cached_summary

from . import config, theme
from .cpu_meter import ProcessMeter
from .help_button import make_help_button
from .spliceai_setup import PACKAGES_DIR, is_frozen, is_spliceai_installed
from .i18n import t

BUILDS = ("hg19", "hg38")

# Measured on this program's own builds (2026-10-27, Ryzen 5 5600): a whole
# build stores 0.854 (hg38) and 0.858 (hg19) bytes per base as float16, and
# while a chromosome is being built its finished genes wait as float32 .npz
# files of 7.4 bytes per base (APP, RUNX1), until the chromosome is packed.
BYTES_PER_BASE = 0.86
WORK_BYTES_PER_BASE = 7.4
# A whole build runs slower than the benchmark (refdb_build --benchmark, about
# 600,000 bases of chromosome 21): on the same PC 37,950 (hg38) and 37,260
# (hg19) bases/s against the benchmark's 42,080.
BENCHMARK_TO_BUILD = 0.89
# Without a benchmark: this program's own whole hg38 build ran at 37,950 bases/s
# (Ryzen 5 5600). Loading the models the first time adds about 10 s.
MEASURED_BUILD_SPEED = 37_950
MODEL_LOAD_SECONDS = 10


def live_seconds(bases):
    """About how long computing a gene of `bases` stored bases live takes on
    this PC: its benchmark speed (from the wizard) if measured, else the speed
    measured here."""
    path = os.path.join(refdb.default_root(), "benchmark.json")
    speed = MEASURED_BUILD_SPEED
    try:
        speed = refdb.read_json(path)["bases_per_second"] * BENCHMARK_TO_BUILD
    except (OSError, ValueError, KeyError):
        pass
    return bases / speed + MODEL_LOAD_SECONDS

BUILD_PROCESS_FLAGS = (subprocess.CREATE_NO_WINDOW | subprocess.BELOW_NORMAL_PRIORITY_CLASS
                       if sys.platform == "win32" else 0)

DATABASE_HELP_TEXT = (
    t("The Sequence to cDNA tab requires SpliceAI's acceptor and donor probabilities at every base of the "
      "reference genome's genes. As predicting them takes SpliceAI several minutes per gene, they are computed "
      "once and stored in a reference score database: one per genome build (hg19, hg38), about 1.1 GB each, in "
      "SpliceAI_reference_data\\spliceai-refdb in the user folder, next to the reference genomes.\n"
      "\n"
      "When the program is installed from the USB drive, Setup copies the database. Otherwise, it can be obtained "
      "in two ways:\n"
      "  • From a folder: a database folder obtained separately (for example spliceai-refdb\\hg19 on the USB "
      "drive) is installed after being checked.\n"
      "  • Compute it on this computer: SpliceAI computes the database from the reference genome. All genes take "
      "many hours per build (about 9 on a 6-core computer); meanwhile the program remains usable (the VCF tab is "
      "not affected), and the computation can be paused or the program closed, as it resumes at the next start.\n"
      "\n"
      "Before use, each database is checked: it must be for the selected build, must have been produced with the "
      "same SpliceAI model files as those installed, and every file must match its checksum. A database computed "
      "on this computer is also compared with a small validation set of known values and is not used if it does "
      "not match.\n"
      "\n"
      "Genes not yet in the database are computed live when needed (a minute or two each), cached and marked "
      "\"computed live\".\n"
      "\n"
      "As the database is derived from SpliceAI's models, SpliceAI's terms apply to it: academic and "
      "non-commercial use only.")
)

BUILD_OPTIONS_HELP_TEXT = (
    t("CPU threads: the number of processor threads the computation may use. \"All\" is fastest; fewer threads "
      "leave more of the computer free for other work (the computation already runs at low priority, so the "
      "program and other applications remain responsive in either case).\n"
      "\n"
      "Keep this computer awake: prevents Windows from going to sleep while the computation runs (the screen may "
      "still turn off). This is released when the computation finishes, is paused or stopped, or the program "
      "closes.\n"
      "\n"
      "A Windows Update restart interrupts the computation. No data are lost: it resumes from the last completed "
      "gene the next time the program starts.")
)

CACHE_HELP_TEXT = (
    t("Genes not (yet) in the database are computed live when an analysis needs them, using the same method as the "
      "database, and stored here so that each is computed only once. They are marked \"computed live\". Clearing "
      "the cache only means that they will be computed again.")
)

VALIDATION_HELP_TEXT = (
    t("The validation set contains SpliceAI's values for about 50 genes per build, computed when this program's "
      "own databases were built and checked. After a computation on this computer, the genes it shares with the "
      "validation set are compared with these values; a database that does not match is not used. The validation "
      "set is a separate file of a few MB (it may not be included with the program for licensing reasons); it can "
      "be selected here if available.")
)


def builder_argv(args):
    """The command that runs refdb_build with `args`: the exe's own
    --build-refdb when installed, the module when running from source."""
    if is_frozen():
        return [sys.executable, "--build-refdb", *args]
    return [sys.executable, "-W", "ignore", "-m", "spliceai_pipeline.refdb_build", *args]


def refseq_dir():
    """The folder with the curated RefSeq tables (transcripts.default_refseq_dir)."""
    return transcripts.default_refseq_dir()


def mane_dir(chosen):
    """The MANE folder the window uses, else installer/thirdparty/mane of a source checkout."""
    for folder in (chosen, os.path.join(app_root_dir(), "installer", "thirdparty", "mane")):
        if folder and find_cached_summary(folder):
            return folder
    return chosen


def validation_set_path(build, chosen=""):
    """The validation set to verify a local build with: the file chosen in the
    window, one shipped with the program, or the one this program's own
    database build wrote; None if there is none."""
    name = f"validation-set-{build}.h5"
    for path in (chosen, os.path.join(app_root_dir(), "validation-sets", name),
                 os.path.join(getattr(sys, "_MEIPASS", app_root_dir()), "validation-sets", name),
                 os.path.join(refdb.default_root(), "validation", build, name)):
        if path and os.path.exists(path):
            return path
    return None


def installed_models():
    """(model SHA-256 dict or None, problem or None): None, None when SpliceAI isn't installed."""
    from spliceai_pipeline import predict

    if not is_spliceai_installed():
        return None, None
    try:
        return predict.verify_models(), None
    except predict.ModelMismatchError as exc:
        return None, str(exc)


def scan(build, models, building=False):
    """What there is for `build`: a dict with state "missing", "refused" (with
    the reason), "building", "ready" (all genes) or "partial", and for a
    usable one its manifest, gene count, size and whether it was verified."""
    folder = os.path.join(refdb.default_root(), build)
    status = {"build": build, "folder": folder, "state": "missing", "message": ""}
    try:
        manifest = refdb.check(folder, build, models, files=not building)
    except refdb.DatabaseRefused as exc:
        return dict(status, state="refused", message=str(exc))
    except (OSError, ValueError) as exc:
        return dict(status, state="refused", message=t("The {build} database in {folder} can't be read: {error}",
                                                       build=build, folder=folder, error=exc))
    if manifest is None:
        return status
    verified = next((r.get("passed") for r in (manifest.get("validation"), manifest.get("local_verification"))
                     if r and r.get("passed")), None)
    complete = manifest.get("scope") == "all genes" and not manifest.get("unpacked_genes")
    return dict(status, state="building" if building else ("ready" if complete else "partial"), manifest=manifest,
                genes=manifest.get("gene_count", 0) + len(manifest.get("unpacked_genes") or []),
                size=refdb.folder_size(folder), verified=verified)


def describe(status):
    """One line for a status, for the dialog and the top bar's tooltip."""
    build, state = status["build"], status["state"]
    if state == "checking":
        return f"{build}: " + t("checking...")
    if state == "missing":
        return f"{build}: " + t("not on this computer yet")
    if state == "refused":
        return f"{build}: " + t("not used -- {why}", why=t(status["message"]))
    m = status["manifest"]
    parts = [t("{n:,} genes (all)", n=status["genes"]) if state == "ready" else t("{n:,} genes", n=status["genes"]),
             f"{status['size'] / 1e9:.2f} GB", t("updated {date}", date=m.get("updated", "?")[:10])]
    if status["verified"]:
        parts.append(t("verified"))
    elif state != "building":
        parts.append(t("not verified"))
    what = {"ready": "✓ " + t("ready"), "partial": t("partly built"), "building": t("being built")}[state]
    return f"{build}: {what} -- " + ", ".join(parts)


def progress_text(p):
    """A build's progress.json in words: where it is, how far overall, the time left."""
    if not p:
        return t("Starting (loading SpliceAI's models)...")
    parts = []
    if p.get("chrom"):
        parts.append(f"chr{p['chrom']}: " + t("{done:,} of {total:,} genes", done=p["chrom_genes_done"],
                                               total=p["chrom_genes_total"]))
    total, done = p.get("bases_total") or 0, p.get("bases_done") or 0
    if total:
        parts.append(t("overall {pct:.1f}%", pct=100 * done / total))
    # The speed of this run, once it has run long enough to say.
    if p.get("seconds_this_run", 0) > 30 and p.get("bases_this_run") and total > done:
        parts.append(t("about {time} left", time=format_elapsed((total - done) * p["seconds_this_run"]
                                                               / p["bases_this_run"])))
    text = {"paused": t("Paused. "), "stopped": t("Stopped. ")}.get(p.get("state"), "") + "; ".join(parts)
    return text + (f" -- now {p['gene']}" if p.get("gene") else "")


def _size(n):
    return f"{n / 1e9:,.1f} GB" if n >= 1e8 else f"{n / 1e6:,.1f} MB"


def estimate(bases, bench, largest_chrom_bases):
    """(seconds or None, bytes on disk at the end, bytes free needed meanwhile)."""
    seconds = bases / (bench * BENCHMARK_TO_BUILD) if bench else None
    final = bases * BYTES_PER_BASE
    return seconds, final, final + largest_chrom_bases * WORK_BYTES_PER_BASE


def _run_in_thread(fn, done):
    """Runs fn() off the window's thread; done(result, error) back on it."""
    relay = _Relay()
    relay.finished.connect(lambda result, error: (done(result, error), relay.deleteLater()), Qt.QueuedConnection)

    def target():
        try:
            result, error = fn(), None
        except Exception as exc:  # noqa: BLE001 -- shown to the user
            result, error = None, exc
        try:
            relay.finished.emit(result, error)
        except RuntimeError:        # the window closed meanwhile: nobody to tell
            pass

    threading.Thread(target=target, daemon=True).start()
    return relay


class _Relay(QObject):
    finished = Signal(object, object)


class RefdbManager(QObject):
    """The databases and the builds, for the whole window. fasta_for(build)
    and mane_dir_for() give the window's current reference FASTA and MANE
    folder."""

    changed = Signal()
    notice = Signal(str, str)   # (title, text), when a build ends or a verification is done

    def __init__(self, fasta_for, mane_dir_for, parent=None):
        super().__init__(parent)
        self.fasta_for, self.mane_dir_for = fasta_for, mane_dir_for
        self.status = {b: {"build": b, "state": "checking", "message": ""} for b in BUILDS}
        self.models = self.model_problem = None
        self.proc = self.job = self.meter = self._log = None
        self.progress, self.stats, self.errors = {}, {}, {}
        self._stopping, self._polls, self._size0 = None, 0, 0
        self._relays, self._readers = [], {}
        self.importing = False
        self.timer = QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self._poll)

    # --- status -----------------------------------------------------------------------
    def refresh(self):
        building = self.job["build"] if self.proc else None

        def work():
            models, problem = installed_models()
            return models, problem, {b: scan(b, models, building=b == building) for b in BUILDS}

        def done(result, error):
            if error is None:
                self.models, self.model_problem, self.status = result
            self.changed.emit()

        self._relays.append(_run_in_thread(work, done))

    def usable(self, build):
        return self.status[build]["state"] in ("ready", "partial", "building")

    def reader(self, build):
        """A refdb.Reader for the build's database, or None; the same one
        until its manifest changes (a full manifest is about 10 MB)."""
        if not self.usable(build):
            return None
        folder = self.status[build]["folder"]
        path = os.path.join(folder, "manifest.json")
        stamp = (folder, os.path.getmtime(path) if os.path.exists(path) else None)
        cached = self._readers.get(build)
        if cached is None or cached[0] != stamp:
            self._readers[build] = (stamp, refdb.Reader(folder))
        return self._readers[build][1]

    def cache_dir(self):
        return os.path.join(refdb.default_root(), "live-cache")

    def clear_cache(self):
        shutil.rmtree(self.cache_dir(), ignore_errors=True)

    def validation_set(self, build):
        return validation_set_path(build, config.load()["refdb_validation_sets"].get(build, ""))

    def validation_genes(self, build):
        """The genes of the build's validation set (none without one)."""
        path = self.validation_set(build)
        if not path:
            return set()
        import h5py

        with h5py.File(path, "r") as vs:
            return set(vs.keys())

    def set_validation_set(self, build, path):
        settings = config.load()
        settings["refdb_validation_sets"] = dict(settings["refdb_validation_sets"], **{build: path})
        config.save(settings)

    # --- the queue of builds ------------------------------------------------------------
    def jobs(self):
        return list(config.load()["refdb_jobs"])

    def _save_jobs(self, jobs):
        settings = config.load()
        settings["refdb_jobs"] = jobs
        config.save(settings)

    def add_jobs(self, jobs, front=False):
        """Queues builds ({"build", "genes": None for all or a list, "threads",
        "keep_awake"}) and starts the first if nothing is running. front: right
        after the running one."""
        queue = self.jobs()
        at = (1 if self.proc else 0) if front else len(queue)
        self._save_jobs(queue[:at] + list(jobs) + queue[at:])
        self.start_next()
        self.changed.emit()

    def prerequisites(self, build):
        """What keeps a build of `build` from starting, or None."""
        fasta = self.fasta_for(build)
        if not fasta or not os.path.isfile(fasta):
            return t("No {build} reference FASTA is set (Settings > Advanced settings).", build=build)
        if not os.path.exists(fasta + ".fai"):
            return t("The {build} reference FASTA has no index ({name}.fai) next to it.", build=build,
                     name=os.path.basename(fasta))
        if not is_spliceai_installed():
            return t("SpliceAI isn't installed (the button at the top of the window).")
        if self.model_problem:
            return self.model_problem
        if not os.path.exists(transcripts.table_path(build, refseq_dir())):
            return t("The curated RefSeq table for {build} isn't there ({folder}).", build=build, folder=refseq_dir())
        if not find_cached_summary(mane_dir(self.mane_dir_for())):
            return t("The MANE Select summary file isn't there (Settings > Advanced settings).")
        return None

    def start_next(self):
        if self.proc is not None:
            return
        queue = self.jobs()
        if not queue:
            return
        job = queue[0]
        build = job["build"]
        problem = self.prerequisites(build)
        if problem:
            self.errors[build] = t("The {build} build can't start: {problem}", build=build, problem=problem)
            self.changed.emit()
            return
        self.errors.pop(build, None)
        folder = os.path.join(refdb.default_root(), build)
        os.makedirs(folder, exist_ok=True)
        for name in ("control.txt", "progress.json"):
            if os.path.exists(os.path.join(folder, name)):
                os.remove(os.path.join(folder, name))
        args = ["--build", build, "--fasta", self.fasta_for(build), "--out", folder,
                "--refseq-dir", refseq_dir(), "--mane-dir", mane_dir(self.mane_dir_for())]
        if job.get("genes"):
            genes_file = os.path.join(folder, "requested-genes.txt")
            with open(genes_file, "w", encoding="utf-8") as fh:
                fh.write("\n".join(job["genes"]) + "\n")
            args += ["--genes-file", genes_file, "--ignore-unknown"]
        if job.get("threads"):
            args += ["--threads", str(job["threads"])]
        if not job.get("keep_awake", True):
            args.append("--no-keep-awake")
        log_path = os.path.join(folder, "builder-output.log")
        self._log = open(log_path, "a", encoding="utf-8")
        env = dict(os.environ, SPLICEAI_CLI_LOG=log_path, TF_CPP_MIN_LOG_LEVEL="2")
        # From source the module runs in a plain Python, which has to be told where SpliceAI is.
        env["PYTHONPATH"] = os.pathsep.join(x for x in (os.environ.get("PYTHONPATH"), str(PACKAGES_DIR)) if x)
        # Installed, the exe writes its own output to SPLICEAI_CLI_LOG (run_app.py).
        self.proc = subprocess.Popen(builder_argv(args), cwd=app_root_dir(), env=env, stdin=subprocess.DEVNULL,
                                     stdout=subprocess.DEVNULL if is_frozen() else self._log,
                                     stderr=subprocess.STDOUT, creationflags=BUILD_PROCESS_FLAGS)
        self.job, self.progress, self.stats, self._stopping, self._polls = job, {}, {}, None, 0
        self.meter = ProcessMeter(self.proc.pid)
        self._size0 = refdb.folder_size(folder)
        self.timer.start()
        self.refresh()

    def _control(self, state):
        if self.proc is not None:
            with open(os.path.join(refdb.default_root(), self.job["build"], "control.txt"), "w", encoding="utf-8") as fh:
                fh.write(state)

    def pause(self):
        self._control("pause")

    def resume(self):
        self._control("run")

    def cancel(self):
        """Stops the running build and takes it off the queue; what it finished stays usable."""
        if self.proc is not None:
            self._stopping = "cancel"
            self._control("stop")

    def shutdown(self, wait=20):
        """The window is closing: stops the build and keeps it queued, to continue at the next start."""
        if self.proc is None:
            return
        self._stopping = "exit"
        self._control("stop")
        try:
            self.proc.wait(wait)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            self.proc.wait(10)
        self._ended(self.proc.returncode)

    def request_gene(self, build, gene):
        """"Build this gene now": next in the running build of that genome
        build, or a build of just this gene ahead of the others."""
        if self.proc is not None and self.job["build"] == build:
            with open(os.path.join(refdb.default_root(), build, "priority.txt"), "a", encoding="utf-8") as fh:
                fh.write(gene + "\n")
        else:
            last = (self.jobs() or [{}])[0]
            self.add_jobs([{"build": build, "genes": [gene], "threads": last.get("threads", 0),
                            "keep_awake": last.get("keep_awake", True)}], front=True)

    def _poll(self):
        if self.proc is None:
            return
        folder = os.path.join(refdb.default_root(), self.job["build"])
        try:
            # start_next() removed any older one, so this is the running build's.
            self.progress = refdb.read_json(os.path.join(folder, "progress.json"))
        except (OSError, ValueError):
            pass
        # From source, the venv's python.exe only launches the real one: measure
        # the process the build itself names.
        pid = self.progress.get("pid")
        if pid and pid != self.meter.pid:
            self.meter.close()
            self.meter = ProcessMeter(pid)
        cpu, ram = self.meter.sample()
        if cpu is not None:
            self.stats.update(cpu=cpu, ram=ram)
        if self._polls % 5 == 0:
            self.stats.update(written=max(0, refdb.folder_size(folder) - self._size0),
                              free=shutil.disk_usage(refdb.default_root()).free)
        self._polls += 1
        code = self.proc.poll()
        if code is not None:
            self._ended(code)
        self.changed.emit()

    def _ended(self, code):
        job, stopping, progress = self.job, self._stopping, self.progress
        self.timer.stop()
        self.meter.close()
        self._log.close()
        self.proc = self.job = self.meter = self._log = None
        self.progress = {}
        build = job["build"]
        folder = os.path.join(refdb.default_root(), build)
        log = os.path.join(folder, "builder-output.log")
        queue = self.jobs()
        if stopping == "exit" or (code == 0 and progress.get("state") == "stopped" and stopping is None):
            pass        # stays queued; continues at the next start
        elif stopping == "cancel":
            self._save_jobs([j for j in queue if j != job])
            self.notice.emit(t("Build stopped"), t("The {build} build was stopped. The genes it finished are kept "
                                                   "and used; build again to continue it.", build=build))
        elif code == 0 and progress.get("state") == "finished":
            self._save_jobs([j for j in queue if j != job])
            self._verify(build, folder)
        else:
            self.errors[build] = t("The {build} build stopped unexpectedly (exit code {code}). Its finished genes "
                                   "are kept; it continues at the next start. Details: {log}", build=build, code=code,
                                   log=log)
            self.notice.emit(t("Build stopped unexpectedly"), self.errors[build])
        self.refresh()
        if stopping != "exit" and build not in self.errors:
            self.start_next()

    def _verify(self, build, folder):
        set_path = self.validation_set(build)
        if set_path is None:
            self.notice.emit(t("Build finished"), t("The {build} database is built. It could not be verified: no "
                                                    "validation set was found (Reference database > Validation set).",
                                                    build=build))
            return

        def done(record, error):
            if error is not None:
                text = t("The {build} database is built, but verifying it failed: {error}", build=build, error=error)
            elif record["passed"] is None:
                text = t("The {build} database is built. None of the validation set's genes are in it, so it "
                         "could not be verified.", build=build)
            else:
                try:
                    refdb.record_verification(folder, record)
                except OSError as exc:      # the folder went meanwhile (deleted, replaced by an import)
                    self.notice.emit(t("Build finished"), t("The {build} database's verification couldn't be "
                                                            "recorded: {error}", build=build, error=exc))
                    self.refresh()
                    return
                text = (t("The {build} database is built and verified: {n} genes match the validation set "
                          "(largest difference {diff:.6f}).", build=build, n=record["genes"],
                          diff=record["max_abs_diff"]) if record["passed"] else
                        t("The {build} database is built but failed its verification, so it won't be used: ",
                          build=build) + "; ".join(record["problems"][:3]))
            self.notice.emit(t("Build finished"), text)
            self.refresh()

        self._relays.append(_run_in_thread(lambda: refdb.verify_against_set(folder, set_path), done))

    def import_database(self, src, on_progress, done):
        """Imports a database folder (refdb.import_database) off the window's thread; done(dest, error). One at a
        time (self.importing): two would share the <build>.importing folder."""
        self.importing = True

        def finished(dest, error):
            self.importing = False
            done(dest, error)
            self.refresh()

        self._relays.append(_run_in_thread(lambda: refdb.import_database(src, self.models, on_progress=on_progress),
                                           finished))


def database_folder_of(path):
    """The database folder a dropped file or folder stands for (a folder with
    manifest.json, or that file itself), or None."""
    folder = os.path.dirname(path) if os.path.basename(path).lower() == "manifest.json" else path
    return folder if os.path.isfile(os.path.join(folder, "manifest.json")) else None


def help_button(parent, title, text):
    return make_help_button(parent, title, text)


class DatabaseActions(QWidget):
    """The two ways to get a database, side by side: take it from a folder
    (the USB drive's spliceai-refdb) or compute it on this computer (the
    wizard); the import's progress below them."""

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self.manager = manager
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.import_button = QPushButton(t("From a folder (USB drive)..."))
        self.import_button.setObjectName("runButton")
        self.import_button.setToolTip(t("Install a database folder obtained separately (e.g. spliceai-refdb\\hg19 on "
                                      "the USB drive) -- or drag the folder here."))
        self.import_button.clicked.connect(self._on_import)
        row.addWidget(self.import_button)
        self.build_button = QPushButton(t("Compute it on this computer..."))
        self.build_button.setToolTip(t("Compute it here with SpliceAI: many hours per genome build."))
        self.build_button.clicked.connect(lambda: (RefdbWizard(manager, self.window()).exec(), manager.refresh()))
        row.addWidget(self.build_button)
        row.addWidget(help_button(self, t("Reference score database"), DATABASE_HELP_TEXT))
        row.addStretch(1)
        layout.addLayout(row)
        progress = QHBoxLayout()
        self.bar = QProgressBar()
        self.bar.setRange(0, 1000)
        progress.addWidget(self.bar, stretch=1)
        layout.addLayout(progress)
        self.label = QLabel()
        self.label.setWordWrap(True)
        self.label.setStyleSheet(theme.STATUS_MUTED)
        layout.addWidget(self.label)
        self._busy(False)

    def _busy(self, on, text=""):
        for w in (self.import_button, self.build_button):
            w.setEnabled(not on)
        self.bar.setVisible(on)
        self.label.setText(text)
        self.label.setVisible(bool(text))
        if on:
            self.bar.setRange(0, 0)

    def _show_progress(self, text, done, total):
        self.label.setText(text)
        self.label.setVisible(True)
        self.bar.setRange(0, 1000)
        self.bar.setValue(int(1000 * done / max(total, 1)))

    def _on_import(self):
        folder = QFileDialog.getExistingDirectory(self, t("Choose the database folder (the one with manifest.json)"))
        if folder:
            self.import_folder(folder)

    def import_folder(self, folder):
        """Imports a database folder, after asking before replacing one."""
        if self.manager.importing:
            QMessageBox.information(self, t("Import database"), t("A database is being imported; wait until it is "
                                                                  "finished."))
            return
        try:
            build = refdb.read_json(os.path.join(folder, "manifest.json")).get("build")
        except (OSError, ValueError):
            build = None
        if build in BUILDS:
            if self.manager.proc is not None and self.manager.job["build"] == build:
                QMessageBox.warning(self, t("Import database"), t("The {build} database is being built right now. "
                                    "Stop the build first, then import.", build=build))
                return
            if self.manager.status[build]["state"] != "missing":
                reply = QMessageBox.question(self, t("Replace the database?"), t("There already is an {build} database. "
                                             "Replace it with the imported one? It is replaced only after the copy "
                                             "is complete and checked.", build=build))
                if reply != QMessageBox.Yes:
                    return
        self._busy(True, t("Copying and checking {folder}...", folder=folder))
        relay = _Relay(self)
        relay.finished.connect(lambda done, total: self._show_progress(t("Copying and checking {folder}...",
                                                                         folder=folder), done, total),
                               Qt.QueuedConnection)

        def finished(dest, error):
            self._busy(False)
            if error is not None:
                QMessageBox.warning(self, t("Import failed"), t(str(error)))
            else:
                QMessageBox.information(self, t("Database imported"), t("The database was checked and installed in {folder}.", folder=dest))

        self.manager.import_database(folder, lambda d, t: relay.finished.emit(d, t), finished)


class MissingDatabasePanel(QFrame):
    """Shown instead of the Sequence-to-cDNA tab's inputs while there is no
    usable database: what it is for, and the ways to get one."""

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        layout = QVBoxLayout(self)
        title = QLabel(t("No reference score database yet"))
        title.setStyleSheet("font-size: 13pt; font-weight: 600;")
        layout.addWidget(title)
        text = QLabel(
            t("This tab compares a sample's sequence with the reference genome using SpliceAI's predictions at every base "
              "of a gene. These predictions for the reference are stored in a database on this computer (about 1.1 GB per "
              "genome build). The database is provided on the program's USB drive (folder spliceai-refdb; Setup run from "
              "the drive copies it) or can be computed on this computer. The VCF tab does not require it.")
        )
        text.setWordWrap(True)
        layout.addWidget(text)
        self.actions = DatabaseActions(manager)
        layout.addWidget(self.actions)


class RefdbDialog(QDialog):
    """Status of both databases, the running build, import, cache and validation sets."""

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self.manager = manager
        self.setWindowTitle(t("Reference score database"))
        self.resize(760, 520)
        self.setAcceptDrops(True)
        layout = QVBoxLayout(self)

        top = QHBoxLayout()
        heading = QLabel(t("Reference score database"))
        heading.setStyleSheet("font-size: 13pt; font-weight: 600;")
        top.addWidget(heading)
        top.addStretch(1)
        top.addWidget(help_button(self, t("Reference score database"), DATABASE_HELP_TEXT))
        layout.addLayout(top)

        self.status_labels = {}
        for build in BUILDS:
            label = QLabel()
            label.setWordWrap(True)
            label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.status_labels[build] = label
            layout.addWidget(label)

        self.build_box = QFrame()
        self.build_box.setFrameShape(QFrame.StyledPanel)
        box = QVBoxLayout(self.build_box)
        self.build_title = QLabel()
        self.build_title.setStyleSheet("font-weight: 600;")
        box.addWidget(self.build_title)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1000)
        box.addWidget(self.progress_bar)
        self.progress_label = QLabel()
        self.progress_label.setWordWrap(True)
        box.addWidget(self.progress_label)
        self.stats_label = QLabel()
        self.stats_label.setStyleSheet(theme.STATUS_MUTED)
        box.addWidget(self.stats_label)
        reminder = QLabel(t("A Windows Update restart interrupts the build; it then continues at the next start."))
        reminder.setStyleSheet(theme.STATUS_MUTED)
        box.addWidget(reminder)
        controls = QHBoxLayout()
        self.pause_button = QPushButton(t("Pause"))
        self.pause_button.clicked.connect(self._on_pause)
        controls.addWidget(self.pause_button)
        self.stop_button = QPushButton(t("Stop"))
        self.stop_button.setToolTip(t("Stop this build. The genes it finished are kept and used."))
        self.stop_button.clicked.connect(self._on_stop)
        controls.addWidget(self.stop_button)
        controls.addStretch(1)
        box.addLayout(controls)
        layout.addWidget(self.build_box)

        self.error_label = QLabel()
        self.error_label.setWordWrap(True)
        self.error_label.setStyleSheet(theme.STATUS_CRIT)
        layout.addWidget(self.error_label)

        self.validation_labels = {}
        for build in BUILDS:
            row = QHBoxLayout()
            label = QLabel()
            label.setWordWrap(True)
            self.validation_labels[build] = label
            row.addWidget(label, stretch=1)
            choose = QPushButton(t("Choose..."))
            choose.setToolTip(t("Choose the {build} validation set file (validation-set-{build}.h5)", build=build))
            choose.clicked.connect(lambda _c=False, b=build: self._on_choose_validation_set(b))
            row.addWidget(choose)
            if build == BUILDS[0]:
                row.addWidget(help_button(self, t("Validation set"), VALIDATION_HELP_TEXT))
            layout.addLayout(row)

        cache_row = QHBoxLayout()
        self.cache_label = QLabel()
        cache_row.addWidget(self.cache_label, stretch=1)
        self.clear_cache_button = QPushButton(t("Clear cache"))
        self.clear_cache_button.clicked.connect(self._on_clear_cache)
        cache_row.addWidget(self.clear_cache_button)
        cache_row.addWidget(help_button(self, t("Genes computed live"), CACHE_HELP_TEXT))
        layout.addLayout(cache_row)

        layout.addStretch(1)
        self.actions = DatabaseActions(manager)
        self.import_folder = self.actions.import_folder
        layout.addWidget(self.actions)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        close = QPushButton(t("Close"))
        close.clicked.connect(self.accept)
        buttons.addWidget(close)
        layout.addLayout(buttons)

        manager.changed.connect(self.update_view)
        self.update_view()

    def update_view(self):
        m = self.manager
        for build in BUILDS:
            status = m.status[build]
            label = self.status_labels[build]
            label.setText(describe(status))
            label.setStyleSheet({"ready": theme.STATUS_OK, "refused": theme.STATUS_CRIT,
                                 "missing": theme.STATUS_MUTED}.get(status["state"], ""))
            path = m.validation_set(build)
            self.validation_labels[build].setText(t("{build} validation set: {path}", build=build,
                                                     path=path or t("none found")))
        running = m.proc is not None
        self.build_box.setVisible(running)
        if running:
            p, s = m.progress, m.stats
            genes = m.job.get("genes")
            self.build_title.setText(t("Building {build} -- {n:,} chosen gene(s)", build=m.job["build"], n=len(genes))
                                     if genes else t("Building {build} -- all genes", build=m.job["build"]))
            total, done = p.get("bases_total") or 0, p.get("bases_done") or 0
            self.progress_bar.setValue(int(1000 * done / total) if total else 0)
            self.progress_label.setText(progress_text(p))
            parts = []
            if s.get("cpu") is not None:
                parts.append(f"CPU {s['cpu']:.0f}%")
            if s.get("ram"):
                parts.append(f"RAM {s['ram'] / 1e9:.1f} GB")
            if "written" in s:
                parts.append(f"written {s['written'] / 1e6:,.0f} MB")
            if "free" in s:
                parts.append(f"free disk {s['free'] / 1e9:,.0f} GB")
            self.stats_label.setText(" · ".join(parts))
            self.pause_button.setText(t("Resume") if p.get("state") == "paused" else t("Pause"))
        self.error_label.setText("\n".join(m.errors.values()))
        self.error_label.setVisible(bool(m.errors))
        size = refdb.folder_size(m.cache_dir())
        self.cache_label.setText(t("Genes computed live (cache): {mb:,.1f} MB", mb=size / 1e6))
        self.clear_cache_button.setEnabled(size > 0)

    def _on_pause(self):
        if self.manager.progress.get("state") == "paused":
            self.manager.resume()
        else:
            self.manager.pause()

    def _on_stop(self):
        reply = QMessageBox.question(self, t("Stop the build?"), t("Stop this build? The genes it finished are kept "
                                     "and used, and building again continues from there."))
        if reply == QMessageBox.Yes:
            self.manager.cancel()

    def _on_clear_cache(self):
        self.manager.clear_cache()
        self.update_view()

    def _on_choose_validation_set(self, build):
        path, _ = QFileDialog.getOpenFileName(self, t("Choose the {build} validation set", build=build), "",
                                              t("Validation set (validation-set-{build}.h5);;HDF5 files (*.h5)",
                                                build=build))
        if path:
            self.manager.set_validation_set(build, path)
            self.update_view()

    def dragEnterEvent(self, event):
        if any(database_folder_of(u.toLocalFile()) for u in event.mimeData().urls() if u.isLocalFile()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            folder = url.isLocalFile() and database_folder_of(url.toLocalFile())
            if folder:
                event.acceptProposedAction()
                self.import_folder(os.path.normpath(folder))
                return


class _ScopePage(QWizardPage):
    def __init__(self, manager):
        super().__init__()
        self.setTitle(t("Compute the database on this computer"))
        self.genes = []
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(t("Genome build(s):")))
        self.build_boxes = {}
        for build in BUILDS:
            box = QCheckBox(build)
            status = manager.status[build]
            box.setChecked(status["state"] in ("missing", "partial") and bool(manager.fasta_for(build)))
            box.setToolTip(describe(status))
            box.toggled.connect(self.completeChanged)
            self.build_boxes[build] = box
            layout.addWidget(box)
        layout.addSpacing(8)
        self.all_radio = QRadioButton(t("All genes (recommended)"))
        self.all_radio.setChecked(True)
        self.some_radio = QRadioButton(t("Only these genes:"))
        self.all_radio.toggled.connect(self.completeChanged)
        layout.addWidget(self.all_radio)
        row = QHBoxLayout()
        row.addWidget(self.some_radio)
        self.choose_button = QPushButton(t("Choose genes..."))
        self.choose_button.clicked.connect(self._on_choose)
        row.addWidget(self.choose_button)
        self.genes_label = QLabel(t("none chosen"))
        row.addWidget(self.genes_label, stretch=1)
        layout.addLayout(row)
        note = QLabel(t("A partial build also includes the about 50 genes of the validation set, so that it can be "
                      "verified. Genes left out are computed live when they are needed."))
        note.setWordWrap(True)
        note.setStyleSheet(theme.STATUS_MUTED)
        layout.addWidget(note)
        layout.addSpacing(8)
        row = QHBoxLayout()
        row.addWidget(QLabel(t("CPU threads:")))
        self.threads = QSpinBox()
        self.threads.setRange(0, os.cpu_count() or 1)
        self.threads.setSpecialValueText(t("All"))
        row.addWidget(self.threads)
        row.addWidget(help_button(self, t("Build options"), BUILD_OPTIONS_HELP_TEXT))
        row.addStretch(1)
        layout.addLayout(row)
        self.keep_awake = QCheckBox(t("Keep this computer awake during the build"))
        self.keep_awake.setChecked(True)
        layout.addWidget(self.keep_awake)
        text = QLabel(t("It runs in the background: the program remains usable, and closing it (or a Windows Update restart) only "
                        "pauses the computation until the next start. Genes can be used as soon as they are completed. Progress is "
                        "shown under Reference DB at the top of the window."))
        text.setWordWrap(True)
        text.setStyleSheet(theme.STATUS_MUTED)
        layout.addWidget(text)
        layout.addStretch(1)

    def _on_choose(self):
        from .gene_import_dialog import GeneImportDialog

        dialog = GeneImportDialog(self)
        dialog.setWindowTitle(t("Genes to build"))
        dialog.apply_button.setText(t("Use these genes"))
        if dialog.exec() == QDialog.Accepted:
            self.genes = dialog.validated_genes()
            self.genes_label.setText(t("{n:,} gene(s):", n=len(self.genes)) + " " + ", ".join(self.genes[:8])
                                     + (", ..." if len(self.genes) > 8 else ""))
            self.some_radio.setChecked(True)
        self.completeChanged.emit()

    def builds(self):
        return [b for b, box in self.build_boxes.items() if box.isChecked()]

    def chosen_genes(self):
        return None if self.all_radio.isChecked() else list(self.genes)

    def isComplete(self):
        return bool(self.builds()) and (self.all_radio.isChecked() or bool(self.genes))


class _CheckPage(QWizardPage):
    """Prerequisites, disk space, and a short benchmark for the time estimate."""

    def __init__(self, manager, scope):
        super().__init__()
        self.manager, self.scope = manager, scope
        self.setTitle(t("Checking this computer"))
        layout = QVBoxLayout(self)
        self.label = QLabel()
        self.label.setWordWrap(True)
        self.label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.label)
        self.busy = QProgressBar()
        self.busy.setRange(0, 0)
        layout.addWidget(self.busy)
        layout.addStretch(1)
        self._ok, self._bench, self._relays = False, None, []

    def initializePage(self):
        self._ok = False
        self.busy.setVisible(True)
        builds, genes = self.scope.builds(), self.scope.chosen_genes()
        problems = [p for p in (self.manager.prerequisites(b) for b in builds) if p]
        if problems:
            self._show(problems, [])
            return
        self.label.setText(t("Reading the gene tables and measuring this computer's speed (about a minute)..."))
        bench_path = os.path.join(refdb.default_root(), "benchmark.json")
        fasta, build = self.manager.fasta_for(builds[0]), builds[0]
        mane = mane_dir(self.manager.mane_dir_for())

        def work():
            sizes = {}
            for b in builds:
                genes_b = transcripts.load_genes(b, refseq_dir(), find_cached_summary(mane))
                wanted = None if genes is None else set(genes) | self.manager.validation_genes(b)
                chosen = [g for g in genes_b if wanted is None or g.name in wanted
                          or (g.canonical and g.canonical.symbol in wanted)]
                per_chrom = {}
                for g in chosen:
                    per_chrom[g.chrom] = per_chrom.get(g.chrom, 0) + g.mask[1] - g.mask[0] + 1 + 2 * refdb.MARGIN
                sizes[b] = (len(chosen), sum(per_chrom.values()), max(per_chrom.values(), default=0))
            os.makedirs(refdb.default_root(), exist_ok=True)
            subprocess.run(builder_argv(["--build", build, "--fasta", fasta, "--benchmark", bench_path,
                                         "--refseq-dir", refseq_dir(), "--mane-dir", mane]),
                           cwd=app_root_dir(), stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL, creationflags=BUILD_PROCESS_FLAGS, check=False,
                           env=dict(os.environ, TF_CPP_MIN_LOG_LEVEL="2"))
            bench = refdb.read_json(bench_path)["bases_per_second"] if os.path.exists(bench_path) else None
            return sizes, bench

        def done(result, error):
            self.busy.setVisible(False)
            if error is not None:
                self._show([t("The check failed: {error}", error=f"{type(error).__name__}: {error}")], [])
                return
            sizes, bench = result
            self._bench = bench
            lines, need = [], 0
            for b, (count, bases, largest) in sizes.items():
                seconds, final, peak = estimate(bases, bench, largest)
                need = max(need, peak)
                lines.append(f"{b}: " + t("{n:,} genes, {bases:,} bases -> about {size}", n=count, bases=bases,
                                           size=_size(final))
                             + (t(", about {time}", time=format_elapsed(seconds)) if seconds else
                                t(", time unknown (the speed measurement failed)")))
            free = shutil.disk_usage(refdb.default_root()).free
            problems = [] if free > need else [t("Not enough free disk space: {need} needed while building, {free} "
                                                 "free on that drive.", need=_size(need), free=_size(free))]
            lines.append(t("Free disk space: {free} ({need} needed while building).", free=_size(free),
                           need=_size(need)))
            if bench:
                lines.append(t("Measured speed: {speed:,} bases/s (a whole build runs somewhat slower; the times above "
                               "allow for that).", speed=bench))
            self._show(problems, lines)

        self._relays.append(_run_in_thread(work, done))

    def _show(self, problems, lines):
        self.busy.setVisible(False)
        ok = ["✓ " + t("Reference FASTA, SpliceAI and the gene tables are in place.")] if not problems else []
        self.label.setText("\n".join(ok + lines + [f"✗ {p}" for p in problems]))
        self.label.setStyleSheet(theme.STATUS_CRIT if problems else "")
        self._ok = not problems
        self.completeChanged.emit()

    def isComplete(self):
        return self._ok


class RefdbWizard(QWizard):
    """The guided build (spec 5.8, Route 2): what to build, then the check of
    this computer (time, disk space) before it starts."""

    def __init__(self, manager, parent=None):
        super().__init__(parent)
        self.manager = manager
        self.setWindowTitle(t("Build the reference score database"))
        self.resize(640, 480)
        self.setWizardStyle(QWizard.ModernStyle)
        self.setButtonText(QWizard.FinishButton, t("Start building"))
        self.scope = _ScopePage(manager)
        for page in (self.scope, _CheckPage(manager, self.scope)):
            self.addPage(page)

    def accept(self):
        genes = self.scope.chosen_genes()
        jobs = []
        for build in self.scope.builds():
            chosen = None if genes is None else sorted(set(genes) | self.manager.validation_genes(build))
            jobs.append({"build": build, "genes": chosen, "threads": self.scope.threads.value(),
                         "keep_awake": self.scope.keep_awake.isChecked()})
        self.manager.add_jobs(jobs)
        super().accept()
