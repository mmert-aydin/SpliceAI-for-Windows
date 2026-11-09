"""Builds the reference score database (spec 5.4): SpliceAI's per-base
probabilities for every gene of its annotation on one build, predicted with
predict.py and stored with refdb.py.

    python -m spliceai_pipeline.refdb_build --build hg38 --fasta hg38.fa --out refdb\\hg38
    SpliceAI-VariantScoring.exe --build-refdb --build hg38 --fasta hg38.fa --out refdb\\hg38

It can be stopped at any point -- closed, killed, the PC restarted -- and run
again with the same arguments: it continues after the last finished gene.
Separate processes can build different chromosomes into one folder at the same
time (--chrom 1 in one, --chrom 2 in another, ...); one process per chromosome.

The window runs it as a process of its own and steers it through files in the
database folder (refdb.py): it reads progress.json, writes "pause", "run" or
"stop" into control.txt, and adds genes to priority.txt to have them built
next, whatever chromosome they are on.
"""
import argparse
import glob
import importlib.metadata
import os
import sys
import time

from . import __version__, predict, refdb, transcripts
from .control import RunCancelled
from .mane import find_cached_summary

DEFAULT_PRECISION = refdb.DEFAULT_PRECISION


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Build the SpliceAI reference score database for one genome build.")
    p.add_argument("--build", choices=("hg19", "hg38"), required=True)
    p.add_argument("--fasta", required=True, help="Reference FASTA of that build (with its .fai)")
    p.add_argument("--out", help="Database folder (created; an existing build in it is continued)")
    p.add_argument("--chrom", action="append",
                   help="Only this chromosome, e.g. 21 (repeat for more). Default: all, in annotation order.")
    p.add_argument("--genes", help="Only these genes, comma-separated (SpliceAI's or RefSeq's symbols)")
    p.add_argument("--genes-file", help="Only the genes in this file (one per line, or separated by commas)")
    p.add_argument("--ignore-unknown", action="store_true",
                   help="Leave out genes SpliceAI's gene table doesn't know (logged) instead of stopping")
    p.add_argument("--benchmark", metavar="JSON",
                   help="Build nothing: time the models on a fixed sample of genes and write the speed to JSON")
    p.add_argument("--threads", type=int, default=0, help="CPU threads for TensorFlow (default: all)")
    p.add_argument("--precision", choices=sorted(refdb.PRECISIONS), default=DEFAULT_PRECISION,
                   help=f"How probabilities are stored (default {DEFAULT_PRECISION})")
    p.add_argument("--chunk", type=int, default=predict.DEFAULT_CHUNK,
                   help=f"Output bases per model call (default {predict.DEFAULT_CHUNK:,}; results don't depend on it)")
    p.add_argument("--refseq-dir", default=None, help="Folder with <build>.ncbiRefSeqCurated.txt.gz")
    p.add_argument("--mane-dir", default=None, help="Folder with the MANE Select summary file")
    p.add_argument("--no-keep-awake", dest="keep_awake", action="store_false",
                   help="Let Windows sleep as usual during the build (default: kept awake; the display may still turn off)")
    return p.parse_args(argv)


def _bare(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001


def keep_awake(on):
    """Asks Windows not to put the PC to sleep while this process works (the
    display may still turn off), or lets go again. Returns the state before the
    call (0 off Windows). Windows ends it with the process anyway, and nothing
    stops a Windows Update restart; a build continues at the next start."""
    if sys.platform != "win32":
        return 0
    import ctypes
    call = ctypes.windll.kernel32.SetThreadExecutionState
    call.argtypes, call.restype = [ctypes.c_uint32], ctypes.c_uint32  # EXECUTION_STATE is a DWORD
    return call(ES_CONTINUOUS | (ES_SYSTEM_REQUIRED if on else 0))


class Progress:
    """progress.json, for the window: where the build is, rewritten after every
    gene. Paused time doesn't count toward seconds_this_run (the speed)."""

    def __init__(self, db, build, chroms, bases_total, bases_done):
        self.db, self._t0, self._paused, self._pause_start = db, time.time(), 0.0, None
        self.data = {"state": "running", "build": build, "pid": os.getpid(), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
                     "chroms": chroms, "chroms_done": [], "chrom": None, "chrom_genes_done": 0, "chrom_genes_total": 0,
                     "gene": None, "bases_total": bases_total, "bases_done": bases_done, "bases_this_run": 0,
                     "seconds_this_run": 0.0, "message": ""}
        self.update()

    def update(self, **fields):
        self.data.update(fields)
        now = time.time()
        paused = self._paused + (now - self._pause_start if self._pause_start else 0.0)
        self.data["seconds_this_run"] = round(now - self._t0 - paused, 1)
        self.data["updated"] = time.strftime("%Y-%m-%d %H:%M:%S")
        refdb.write_json(self.db.path("progress.json"), self.data)

    def add(self, bases, new=False):
        """A gene finished; new: it wasn't counted in bases_total."""
        self.data["bases_total"] += bases if new else 0
        self.data["bases_done"] += bases
        self.data["bases_this_run"] += bases

    def state(self, state):
        if state == "paused":
            self._pause_start = time.time()
        elif self._pause_start is not None:
            self._paused += time.time() - self._pause_start
            self._pause_start = None
        self.update(state=state)


class FileControl:
    """The window's Pause and Stop for a build in its own process, with
    RunControl's check() (spliceai_pipeline.control), called between genes and
    model calls: blocks while control.txt says "pause" -- letting the PC sleep
    meanwhile if the build keeps it awake -- and raises RunCancelled once it
    says "stop"."""

    def __init__(self, path, keeping_awake, on_state):
        self.path, self.keeping_awake, self.on_state = path, keeping_awake, on_state

    def read(self):
        try:
            with open(self.path, encoding="utf-8") as fh:
                return fh.read().strip()
        except OSError:
            return "run"

    def check(self):
        state = self.read()
        if state == "pause":
            self.on_state("paused")
            if self.keeping_awake:
                keep_awake(False)
            while state == "pause":
                time.sleep(0.5)
                state = self.read()
            if self.keeping_awake:
                keep_awake(True)
            self.on_state("running")
        if state == "stop":
            raise RunCancelled("Stopped on request.")


def take_priority(db):
    """The genes the window asked to have built next (priority.txt), which is emptied."""
    path = db.path("priority.txt")
    if not os.path.exists(path):
        return []
    try:
        refdb.replace(path, path + ".taking")
    except OSError:
        return []
    with open(path + ".taking", encoding="utf-8") as fh:
        names = [line.strip() for line in fh if line.strip()]
    os.remove(path + ".taking")
    return names


def _by_symbol(genes):
    """Genes by SpliceAI's symbol, and by their canonical transcript's RefSeq one."""
    by = {g.name: g for g in genes}
    for g in genes:
        if g.canonical:
            by.setdefault(g.canonical.symbol, g)
    return by


def benchmark(args, genes):
    """Times the models on a fixed sample of chromosome 21's genes (about
    600,000 bases of all sizes up to 150 kb; the first gene only warms the
    models up) and writes the speed to args.benchmark. The window turns it into
    an estimate for the whole build (refdb_wizard.py)."""
    from .normalize import resolve_contig

    t0 = time.time()
    predictor = predict.Predictor(args.fasta, args.build)
    load = time.time() - t0
    fasta = predictor.annotator.ref_fasta
    record = fasta[resolve_contig("21", set(fasta.keys()))]
    length = len(record)
    sample, bases = [], 0
    for g in sorted((g for g in genes if g.chrom == "21"), key=lambda g: g.start)[::6]:
        plan = refdb.plan_gene(g, length)
        n = plan["end"] - plan["start"] + 1
        if n <= 150_000:
            sample.append((g, plan))
            bases += n if len(sample) > 1 else 0
            if bases >= 600_000:
                break

    def run(g, plan):
        predictor.predict(lambda a, b: record[a - 1:b].seq, length, g.strand, *plan["mask"], plan["start"],
                          plan["end"], chunk=args.chunk)

    run(*sample[0])
    t0 = time.time()
    for g, plan in sample[1:]:
        run(g, plan)
    seconds = time.time() - t0
    refdb.write_json(args.benchmark, {"build": args.build, "genes": [g.name for g, _ in sample[1:]], "bases": bases,
                                      "seconds": round(seconds, 2), "bases_per_second": round(bases / seconds),
                                      "model_load_seconds": round(load, 1), "threads": args.threads})
    return 0


def main(argv=None):
    args = parse_args(argv)
    if not args.out and not args.benchmark:
        sys.exit("--out is required (the database folder).")
    if args.threads:
        # Has to happen before TensorFlow starts its thread pools.
        import tensorflow as tf
        tf.config.threading.set_intra_op_parallelism_threads(args.threads)
        tf.config.threading.set_inter_op_parallelism_threads(args.threads)

    model_sha256 = predict.verify_models()
    mane_path = find_cached_summary(args.mane_dir)
    table_path = transcripts.table_path(args.build, args.refseq_dir)
    annotation_path = transcripts.annotation_path(args.build)
    for path, what in ((mane_path, "MANE Select summary"), (table_path, "curated RefSeq table")):
        if not path or not os.path.exists(path):
            sys.exit(f"No {what} found ({path or args.mane_dir}) -- see installer/fetch-thirdparty.ps1.")
    genes = transcripts.load_genes(args.build, args.refseq_dir, mane_path)
    if args.benchmark:
        return benchmark(args, genes)

    import keras
    import tensorflow as tf

    settings = {
        "build": args.build,
        "precision": args.precision,
        "format": refdb.FORMAT_VERSION,
        "software_version": __version__,
        "spliceai_version": importlib.metadata.version("spliceai"),
        "tensorflow_version": tf.__version__,
        "keras_version": keras.__version__,
        "model_sha256": model_sha256,
        "annotation": {os.path.basename(p): refdb.sha256_of(p) for p in (annotation_path, table_path, mane_path)},
        "fasta_fai_sha256": refdb.sha256_of(args.fasta + ".fai"),
        "built_locally": True,
        "validation": None,
    }
    db = refdb.Database(args.out)
    try:
        db.start(settings)
    except refdb.DatabaseMismatchError as exc:
        sys.exit(str(exc))

    wanted = None
    if args.genes or args.genes_file:
        text = args.genes or ""
        if args.genes_file:
            with open(args.genes_file, encoding="utf-8-sig") as fh:
                text += "," + fh.read()
        wanted = {s.strip() for s in text.replace("\n", ",").split(",") if s.strip()}
        unknown = sorted(wanted - set(_by_symbol(genes)))
        if unknown and not args.ignore_unknown:
            sys.exit(f"Not in SpliceAI's {args.build} gene table, under either symbol: {', '.join(unknown)}")
        if unknown:
            db.log(f"Not in SpliceAI's {args.build} gene table, left out: {', '.join(unknown)}")
            wanted -= set(unknown)
    chroms = [_bare(c) for c in args.chrom] if args.chrom else list(dict.fromkeys(g.chrom for g in genes))

    db.log(f"Building {args.build} into {args.out}: chromosomes {', '.join(chroms)}"
           + (f", {len(wanted)} gene(s)" + (f" ({', '.join(sorted(wanted))})" if len(wanted) <= 20 else "") if wanted is not None else "")
           + f", {args.precision}")
    if args.keep_awake:
        keep_awake(True)
        db.log("Keeping the PC awake until the build ends (the display may still turn off)")
    try:
        return _build(args, db, genes, wanted, chroms)
    finally:
        if args.keep_awake:
            keep_awake(False)


def _build(args, db, genes, wanted, chroms):
    t_model = time.time()
    predictor = predict.Predictor(args.fasta, args.build)
    fasta = predictor.annotator.ref_fasta
    db.log(f"Models loaded in {time.time() - t_model:.1f} s")
    from .normalize import resolve_contig
    contigs = set(fasta.keys())
    all_names = [g.name for g in genes]
    by_symbol = _by_symbol(genes)
    contig = {c: resolve_contig(c, contigs) for c in dict.fromkeys(g.chrom for g in genes)}
    length = {c: len(fasta[name]) for c, name in contig.items() if name}
    packed = {}

    def packed_genes(chrom):
        if chrom not in packed:
            packed[chrom] = (db.packed(chrom) or {"genes": {}})["genes"]
        return packed[chrom]

    def size(g):
        plan = refdb.plan_gene(g, length[g.chrom])
        return plan["end"] - plan["start"] + 1

    scope = [g for g in genes if g.chrom in chroms
             and (wanted is None or g.name in wanted or (g.canonical and g.canonical.symbol in wanted))]
    buildable = [g for g in scope if g.chrom in length]
    buildable_names = {g.name for g in buildable}
    progress = Progress(db, args.build, chroms, sum(map(size, buildable)),
                        sum(size(g) for g in buildable if g.name in packed_genes(g.chrom) or db.has_work(g.chrom, g.name)))
    control = FileControl(db.path("control.txt"), args.keep_awake, progress.state)
    asked_chroms = set()

    def build_one(g, label):
        plan, t0 = refdb.plan_gene(g, length[g.chrom]), time.time()
        record = fasta[contig[g.chrom]]
        probs = predictor.predict(lambda a, b: record[a - 1:b].seq, length[g.chrom], g.strand, plan["mask"][0],
                                  plan["mask"][1], plan["start"], plan["end"], chunk=args.chunk, checkpoint=control.check)
        db.write_work(g.chrom, g.name, probs, plan)
        n, dt = plan["end"] - plan["start"] + 1, time.time() - t0
        db.log(f"{label} {g.name}: {n:,} bases in {dt:.1f} s ({n / dt:,.0f} bases/s)")
        return n

    def build_asked():
        """Genes the window asked for ("Build this gene now"), before the next one in line."""
        names = take_priority(db)
        for i, name in enumerate(names):
            g = by_symbol.get(name)
            if g is None or g.chrom not in length:
                db.log(f"Asked to build {name} first, but SpliceAI's {args.build} gene table or the FASTA lacks it")
                continue
            if g.name in packed_genes(g.chrom) or db.has_work(g.chrom, g.name):
                continue
            progress.update(gene=g.name, message=f"{g.name}, on request")
            try:
                n = build_one(g, f"chr{g.chrom} on request")
            except RunCancelled:
                # Still wanted next time.
                with open(db.path("priority.txt"), "a", encoding="utf-8") as fh:
                    fh.write("".join(f"{x}\n" for x in names[i:]))
                raise
            progress.add(n, new=g.name not in buildable_names)
            progress.update(message="")
            asked_chroms.add(g.chrom)

    try:
        for chrom in chroms:
            todo = [g for g in scope if g.chrom == chrom]
            if not todo:
                continue
            if chrom not in length:
                for g in todo:
                    db.skip(g.name, chrom, "chromosome not in the reference FASTA")
                db.log(f"chr{chrom}: not in the FASTA -- {len(todo)} gene(s) skipped")
                continue
            missing = [g for g in todo if g.name not in packed_genes(chrom)]
            if not missing:
                db.log(f"chr{chrom}: already built ({len(todo)} gene(s))")
                continue
            plans = {g.name: refdb.plan_gene(g, length[chrom]) for g in missing}
            resumed = sum(db.has_work(chrom, g.name) for g in missing)
            db.log(f"chr{chrom}: {len(missing)} gene(s) to build" + (f", {resumed} already finished" if resumed else ""))
            progress.update(chrom=chrom, chrom_genes_total=len(missing), chrom_genes_done=resumed)
            t_chrom, bases = time.time(), 0
            for n, g in enumerate(missing, 1):
                control.check()
                build_asked()
                if db.has_work(chrom, g.name):
                    continue
                progress.update(gene=g.name)
                size_ = build_one(g, f"chr{chrom} {n}/{len(missing)}")
                bases += size_
                progress.add(size_)
                progress.update(chrom_genes_done=n)
            dt = time.time() - t_chrom
            t_pack = time.time()
            db.pack(chrom, plans)
            packed.pop(chrom, None)
            db.write_manifest(all_names)
            db.log(f"chr{chrom} done: {bases:,} bases predicted in {dt:.0f} s"
                   + (f" ({bases / dt:,.0f} bases/s)" if bases and dt else "") + f", packed in {time.time() - t_pack:.1f} s")
            progress.update(chroms_done=progress.data["chroms_done"] + [chrom], gene=None)
        build_asked()
        # Genes built on request, now or by a run that was stopped, on a chromosome
        # that is already packed or not part of this build. A run limited to some
        # chromosomes leaves the others alone: another process may be building them.
        leftover = asked_chroms if args.chrom else {os.path.basename(d)[3:] for d in glob.glob(db.path("_work", "chr*"))}
        for chrom in sorted(leftover):
            work = glob.glob(db.path("_work", f"chr{chrom}", "*.npz"))
            if work:
                db.pack(chrom, {os.path.basename(p)[:-4]: None for p in work})
    except RunCancelled:
        db.log("Stopped on request; run the same command again to continue")
        progress.state("stopped")
        return 0
    except BaseException as exc:
        progress.update(state="failed", message=f"{type(exc).__name__}: {exc}")
        raise
    manifest = db.write_manifest(all_names)
    db.log(f"Finished: {manifest['gene_count']} gene(s) packed, scope {manifest['scope']}")
    progress.update(state="finished", gene=None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
