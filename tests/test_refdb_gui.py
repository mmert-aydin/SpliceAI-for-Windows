"""The reference score database in the window (spec 5.8, Phase 3 acceptance
(b) and (c) from source): the wizard queues a small build of chosen genes
(plus its check: prerequisites, disk space, the speed benchmark), the build
runs in the background and is paused, the program is closed mid-build and
started again, the build continues where it stopped and passes its
verification; then "Build this gene now", import (also by dropping the
folder onto the window), the status in the top bar and the cache. Real
window (Qt offscreen), throw-away home folder.

    .venv\\Scripts\\python tests\\test_refdb_gui.py
"""
import glob
import json
import os
import shutil
import sys
import time

import helpers
from helpers import check, skip

helpers.watchdog(1200)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg38",
               "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True,
               "mane_dir": str(helpers.MANE_DIR)}, fh)
# The build process finds SpliceAI the way this one does.
os.environ["PYTHONPATH"] = str(helpers.SPLICEAI_PKG) + os.pathsep + os.environ.get("PYTHONPATH", "")

import h5py  # noqa: E402
from PySide6.QtCore import QMimeData, QPoint, QPointF, Qt, QUrl  # noqa: E402
from PySide6.QtGui import QDragEnterEvent, QDropEvent  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

from spliceai_pipeline import refdb, transcripts  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui import config  # noqa: E402
from spliceai_gui import refdb_wizard  # noqa: E402
from spliceai_gui.main_window import MainWindow  # noqa: E402

# Every question is answered Yes (closing during a build, replacing a database).
QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Yes)
warnings = []
QMessageBox.warning = staticmethod(lambda _p, title, text, *a, **k: warnings.append((title, text)) or QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)

root = refdb.default_root()
folder = os.path.join(root, "hg38")
check("the databases live in the (throw-away) home's SpliceAI_reference_data",
      root.startswith(str(home)) and folder.endswith(os.path.join("spliceai-refdb", "hg38")))


def wait(what, until, seconds=300):
    deadline = time.time() + seconds
    while time.time() < deadline:
        app.processEvents()
        if until():
            return True
        time.sleep(0.05)
    check(f"waiting for {what}", False, f"not within {seconds} s")
    return False


def settled(window):
    return all(s["state"] != "checking" for s in window.refdb.status.values())


# The validation set: two small genes of the real hg38 one (Phase 2), if this PC built it.
real_set = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "validation", "hg38",
                        "validation-set-hg38.h5")
vs_path = os.path.join(str(home), "validation-set-hg38.h5")
validation_genes = ["C21orf140", "CAPZA3"]
if os.path.exists(real_set):
    with h5py.File(real_set, "r") as src, h5py.File(vs_path, "w") as dst:
        dst.attrs.update(dict(src.attrs))
        for name in validation_genes:
            src.copy(src[name], dst, name=name)
    real = True
else:
    skip("verification against the Phase 2 values", "this PC has no hg38 validation set; a fresh one is used")
    scores = refdb.ReferenceScores("hg38", str(helpers.HG38), cache_dir=os.path.join(str(home), "vs-cache"),
                                   refseq_dir=str(helpers.THIRDPARTY / "refseq"), mane_dir=str(helpers.MANE_DIR))
    with h5py.File(vs_path, "w") as dst:
        dst.attrs.update(build="hg38", precision="float16", model_sha256=json.dumps(scores.settings()["model_sha256"]))
        for name in validation_genes:
            values, plan, _ = scores.gene(name)
            group = dst.create_group(name)
            for i, track in enumerate(("acceptor", "donor")):
                group.create_dataset(track, data=refdb.encode(values[:, i], "float16"))
            group.attrs["plan"] = json.dumps(plan)
    real = False

genes = transcripts.load_genes("hg38", str(helpers.THIRDPARTY / "refseq"),
                               str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz"))
small = sorted((g for g in genes if g.chrom == "21" and 4_000 <= g.mask[1] - g.mask[0] <= 25_000), key=lambda g: g.start)
# Long enough on one CPU thread (about a minute) to be paused and closed mid-build.
big = sorted((g for g in genes if g.chrom == "21" and 60_000 <= g.mask[1] - g.mask[0] <= 150_000), key=lambda g: g.start)
chosen = [g.name for g in big[:4]]

# --- nothing there yet ----------------------------------------------------------------
window = MainWindow()
wait("the database check", lambda: settled(window), 60)
check("top bar: no reference database yet", window.refdb_status_button.text() == "No reference DB")
dialog = refdb_wizard.RefdbDialog(window.refdb, window)
check("dialog: a folder (the USB drive) or computing it are offered, once; nothing is downloaded",
      "folder" in dialog.actions.import_button.text() and "Compute" in dialog.actions.build_button.text()
      and not any("Download" in b.text() or "Build reference database" in b.text()
                  for b in dialog.findChildren(type(dialog.actions.build_button))))
window.refdb.set_validation_set("hg38", vs_path)
dialog.update_view()
check("dialog: the validation set chosen is shown", vs_path in dialog.validation_labels["hg38"].text())

# --- the wizard ---------------------------------------------------------------------------
wizard = refdb_wizard.RefdbWizard(window.refdb, window)
wizard.show()
scope = wizard.scope
scope.build_boxes["hg19"].setChecked(False)
scope.build_boxes["hg38"].setChecked(True)
scope.genes = list(chosen)
scope.some_radio.setChecked(True)
check("wizard step 1: the choice is complete with a build and genes", scope.isComplete())
wizard.next()
check_page = wizard.currentPage()
wait("the wizard's check (gene tables, disk space, benchmark)", check_page.isComplete, 240)
text = check_page.label.text()
check("wizard step 2: prerequisites in place, size, time and free disk space estimated from a measured speed",
      "✓ Reference FASTA, SpliceAI and the gene tables" in text and "GB" in text and "Measured speed" in text
      and "about" in text, text[:300])
check("wizard step 2: the partial build includes the validation set's genes in its estimate",
      f"hg38: {len(chosen) + len(validation_genes)} genes" in text, text[:120])
wizard.scope.threads.setValue(1)
wizard.accept()
job = (window.refdb.jobs() or [{}])[0]
check("wizard, Start building: queues the build -- chosen genes plus the validation set's -- and starts it",
      job.get("build") == "hg38" and job.get("genes") == sorted(chosen + validation_genes) and job.get("threads") == 1
      and window.refdb.proc is not None, str(job)[:200])

# --- running, paused, the program closed mid-build -----------------------------------------
work = lambda: glob.glob(os.path.join(folder, "_work", "*", "*.npz"))
wait("the first finished gene", lambda: work() and window.refdb.progress.get("bases_done"), 240)
window.refdb.pause()
wait("the build to pause", lambda: window.refdb.progress.get("state") == "paused", 60)
check("top bar: the build's progress, and that it's paused",
      window.refdb_status_button.text().startswith("Building reference DB") and "(paused)" in window.refdb_status_button.text(),
      window.refdb_status_button.text())
dialog.update_view()
check("dialog: progress per chromosome and overall, CPU, RAM, disk written and free",
      "chr" in dialog.progress_label.text() and "overall" in dialog.progress_label.text()
      and "RAM" in dialog.stats_label.text() and "free disk" in dialog.stats_label.text()
      and dialog.pause_button.text() == "Resume", f"{dialog.progress_label.text()} | {dialog.stats_label.text()}")
window.refdb.resume()
wait("the build to run again", lambda: window.refdb.progress.get("state") == "running", 60)
dialog.close()
window.close()
done_before = len(work())       # counted once the build has stopped: a gene can still finish while it stops
check("closing the program mid-build stops the build and keeps it queued",
      window.refdb.proc is None and config.load()["refdb_jobs"] == [job] and work(), f"{len(work())} genes finished")

# --- started again: the build continues and is verified ------------------------------------
window = MainWindow()
check("started again: the queued build continues by itself", window.refdb.proc is not None)
wait("the build to finish", lambda: window.refdb.proc is None, 300)
wait("the verification, and the status after it", lambda: settled(window) and window.refdb.status["hg38"].get("verified"), 120)
manifest = refdb.load_manifest(folder)
record = manifest.get("local_verification") or {}
log = open(os.path.join(folder, "build.log"), encoding="utf-8").read()
check("...continuing after the genes already finished, not redoing them",
      f"{done_before} already finished" in log or done_before == 0, [ln for ln in log.splitlines() if "to build" in ln][-1:])
check("finished: every chosen gene and the validation set's are in the database, the queue is empty",
      sorted(manifest["genes"]) == sorted(chosen + validation_genes) and not config.load()["refdb_jobs"])
check("finished: verified against the validation set" + (" (the Phase 2 values)" if real else ""),
      record.get("passed") is True and record.get("genes") == len(validation_genes), str(record)[:200])
status = window.refdb.status["hg38"]
check("status: hg38 partly built and verified; top bar says it's there",
      status["state"] == "partial" and status["verified"] and window.refdb_status_button.text() == "Reference DB: hg38",
      refdb_wizard.describe(status))
check("the build's notice was shown (not modal)", any("verified" in b.text() for b in window._refdb_notices),
      [b.text()[:80] for b in window._refdb_notices])

# --- "Build this gene now" -------------------------------------------------------------------
extra = small[6].name
window.refdb.request_gene("hg38", extra)
check("Build this gene now, with no build running: a build of just that gene starts", window.refdb.proc is not None
      and window.refdb.job["genes"] == [extra])
notices = len(window._refdb_notices)
wait("that gene, and its verification", lambda: window.refdb.proc is None and len(window._refdb_notices) > notices and settled(window), 240)
check("...and the gene is in the database", extra in (refdb.load_manifest(folder) or {}).get("genes", {}))
reader = window.refdb.reader("hg38")
check("the window's reader reads it", reader is not None and reader.values(extra).shape[1] == 2)

# --- import -------------------------------------------------------------------------------------
copy = os.path.join(str(home), "from-usb", "hg38")
shutil.copytree(folder, copy, ignore=shutil.ignore_patterns("_work", "progress.json", "control.txt", "checked.json"))
shutil.rmtree(folder)
window.refdb.refresh()
wait("the check after removing the database", lambda: settled(window) and window.refdb.status["hg38"]["state"] == "missing", 60)
opened = []
real_open = window._open_refdb_dialog
window._open_refdb_dialog = lambda import_folder=None: opened.append(import_folder)
mime = QMimeData()
mime.setUrls([QUrl.fromLocalFile(copy)])
QApplication.sendEvent(window, QDragEnterEvent(QPoint(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier))
QApplication.sendEvent(window, QDropEvent(QPointF(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier))
check("a database folder dropped onto the window is offered for import", opened == [os.path.normpath(copy)], str(opened))
window._open_refdb_dialog = real_open
dialog = refdb_wizard.RefdbDialog(window.refdb, window)
dialog.import_folder(copy)
wait("the import", lambda: dialog.actions.import_button.isEnabled() and settled(window)
     and window.refdb.status["hg38"]["state"] == "partial", 120)
check("import: checked and installed, then used", window.refdb.status["hg38"]["state"] == "partial"
      and refdb.check(folder, "hg38", window.refdb.models) is not None and not warnings, str(warnings)[:200])

# --- the cache -----------------------------------------------------------------------------------
cache = window.refdb.cache_dir()
os.makedirs(cache, exist_ok=True)
with open(os.path.join(cache, "x.npz"), "wb") as fh:
    fh.write(b"0" * 5000)
dialog.update_view()
check("dialog: the cache's size is shown", "MB" in dialog.cache_label.text() and dialog.clear_cache_button.isEnabled())
dialog.clear_cache_button.click()
check("Clear cache empties it", refdb.folder_size(cache) == 0 and not dialog.clear_cache_button.isEnabled())

dialog.close()
window.close()
helpers.finish()
