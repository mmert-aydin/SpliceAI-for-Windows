"""Reference score database: SpliceAI's per-base acceptor and donor
probabilities of the reference genome, per gene (spec 5.2, 5.4, 5.5).

One folder per build:

    build.json        settings the build started with; a resumed build must match
    manifest.json     what the folder holds (spec 5.5), rewritten after each chromosome
    chr21.h5          one group per gene, datasets "acceptor" and "donor"
    chr21.json        that file's size, SHA-256 and genes; manifest.json is assembled
                      from these, so builds of different chromosomes can run at once
    skipped.tsv       genes not built, with the reason
    build.log         progress and timings
    _work/chr21/      finished genes not packed yet, one .npz each (atomic rename)

Each gene is stored from 50 bases before its mask to 50 after it, so that every
101-base delta-score window of a variant inside the gene is covered. Within the
mask the bases are real sequence; outside it they are N, the way
get_delta_scores masks a gene (predict.py). Where the mask is SpliceAI's own
interval the values are the VCF tab's; within 5,000 bases inside a boundary the
mask moved they are not, and neither are positions whose 5,000 bases of
context run off the chromosome (N there; the VCF tab skips such variants).
Each gene records the stretch that is comparable with the VCF tab.

While a build runs, the window talks to it through three more files:

    progress.json     where the build is, rewritten after each gene
    control.txt       "pause", "run" or "stop", written by the window
    priority.txt      genes to build next, one per line ("Build this gene now")

and checked.json remembers which packed files were already found intact.

This module writes the database (Database) and reads it for the app (spec
5.7): check() decides whether a database may be used, Reader reads it,
ReferenceScores falls back to a live prediction for genes it lacks, and
verify_against_set() / import_database() are spec 5.8's verification and
import.
"""
import glob
import hashlib
import json
import os
import shutil
import time

import numpy as np

from .predict import CONTEXT
from .texts import tx

FORMAT_VERSION = 1
MARGIN = 50     # score.DIST_VAR: half a delta-score window

# Measured on hg38 chr21 (2026-10-26) against the lossless float32 build:
# float16 is 0.91 bytes per base (about 1.1 GB per genome build), off by at most
# 0.00024, moves no site across the 0.05 transcript threshold, and changes a
# delta score by at most 0.00022 (DP never). Fixed 0.001 steps ("q001") would
# be 7x smaller but move 66 of chr21's sites across that threshold.
DEFAULT_PRECISION = "float16"

# A database built on another PC against a validation set (both float16): each
# side is rounded by up to 2^-12, and the two PCs' float32 results may differ
# a little more (the note in every validation set says the same).
VERIFY_TOLERANCE = 2 * 2 ** -12 + 1e-5

MASK_RULE = (
    "Per gene, input bases outside the mask are N (as get_delta_scores does with SpliceAI's gene "
    "interval). Mask = SpliceAI's gene interval joined with the canonical transcript and every curated "
    "NM_ of the same symbol overlapping either (spliceai_pipeline/transcripts.py); an NM_ sharing exons "
    "with another gene's MANE Select transcript is left out. Bases beyond a chromosome end are N."
)

VCF_COMPARABLE_RULE = (
    "Per gene, vcf_comparable [start, end] is where the stored values come from exactly the input "
    "get_delta_scores builds (equal to the VCF tab's up to float rounding): start = SpliceAI's start + "
    "5,000 if the mask begins before it, else the stored start; end = SpliceAI's end - 5,000 if the mask "
    "ends after it, else the stored end; then limited to 5,001 .. chromosome length - 5,000; null if "
    "empty. Conservative: the model's outermost few hundred bases of context change its output by less "
    "than float32 resolves."
)

# Stored as: numpy dtype, and the step for fixed-point storage (None = floating point).
PRECISIONS = {
    "float32": (np.float32, None),
    "float16": (np.float16, None),
    "q001": (np.uint16, 0.001),
}


class DatabaseMismatchError(RuntimeError):
    """A folder already holds a database built with other settings."""


class DatabaseRefused(RuntimeError):
    """A database the app won't use; the message says why, for the user."""


def default_root():
    """Where the app keeps its databases, one folder per build, next to the
    reference genomes: ~/SpliceAI_reference_data/spliceai-refdb."""
    return os.path.join(os.path.expanduser("~"), "SpliceAI_reference_data", "spliceai-refdb")


def folder_size(path):
    """Bytes in a folder and its subfolders (0 if it doesn't exist)."""
    total = 0
    for d, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(d, f))
            except OSError:     # removed meanwhile, e.g. by a running build
                pass
    return total


def plan_gene(gene, chrom_len):
    """What gets stored for a transcripts.Gene on a chromosome chrom_len long:
    a dict that goes into the manifest and the gene's HDF5 attributes."""
    mask_start, mask_end = gene.mask
    start, end = max(1, mask_start - MARGIN), min(chrom_len, mask_end + MARGIN)
    comp_start = gene.start + CONTEXT if mask_start < gene.start else start
    comp_end = gene.end - CONTEXT if mask_end > gene.end else end
    comp_start, comp_end = max(comp_start, CONTEXT + 1), min(comp_end, chrom_len - CONTEXT)
    chrom_end = [[a, b] for a, b in ((start, min(end, CONTEXT)), (max(start, chrom_len - CONTEXT + 1), end)) if a <= b]
    return {
        "gene": gene.name, "chrom": gene.chrom, "strand": gene.strand,
        "start": start, "end": end,
        "mask": [mask_start, mask_end],
        "spliceai": [gene.start, gene.end],
        "vcf_comparable": [comp_start, comp_end] if comp_start <= comp_end else None,
        "chrom_end_padded": chrom_end,
        "canonical": gene.canonical.accession if gene.canonical else None,
        "canonical_symbol": gene.canonical.symbol if gene.canonical else None,
        "rule": gene.rule,
        "transcripts": [t.accession for t in gene.transcripts],
        "excluded": [list(x) for x in gene.excluded],
        "other_mane": list(gene.other_mane),
    }


def encode(probs, precision):
    dtype, step = PRECISIONS[precision]
    return np.rint(probs / step).astype(dtype) if step else probs.astype(dtype)


def decode(data, precision):
    _dtype, step = PRECISIONS[precision]
    return data.astype(np.float32) * np.float32(step) if step else data.astype(np.float32)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def replace(src, dst, tries=50):
    """os.replace, retried for a few seconds: on Windows a reader holding dst
    open (or a virus scanner looking at src) makes it fail for a moment."""
    for attempt in range(tries):
        try:
            return os.replace(src, dst)
        except PermissionError:
            if attempt == tries - 1:
                raise
            time.sleep(0.1)


def write_json(path, data):
    with open(path + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1)
    replace(path + ".tmp", path)


def read_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class Database:
    """A database folder being built (or read back, one gene at a time)."""

    def __init__(self, folder):
        self.folder = folder

    def path(self, *parts):
        return os.path.join(self.folder, *parts)

    def start(self, settings):
        """Begins or resumes a build. settings (build, precision, model and
        annotation hashes, ...) must equal those the folder was started with."""
        os.makedirs(self.path("_work"), exist_ok=True)
        if os.path.exists(self.path("build.json")):
            old = read_json(self.path("build.json"))
            changed = sorted(k for k in set(old) | set(settings) if old.get(k) != settings.get(k))
            if changed:
                raise DatabaseMismatchError(
                    f"{self.folder} holds a database built with other settings ({', '.join(changed)}). "
                    "Use an empty folder, or the same settings to continue that build."
                )
        else:
            write_json(self.path("build.json"), settings)
        for stale in glob.glob(self.path("_work", "*", "*.tmp")):
            os.remove(stale)
        self.settings = settings

    def log(self, line):
        print(line, flush=True)
        with open(self.path("build.log"), "a", encoding="utf-8") as fh:
            fh.write(time.strftime("%Y-%m-%d %H:%M:%S ") + line + "\n")

    def skip(self, gene, chrom, reason):
        new = not os.path.exists(self.path("skipped.tsv"))
        with open(self.path("skipped.tsv"), "a", encoding="utf-8", newline="\n") as fh:
            if new:
                fh.write("gene\tchrom\treason\n")
            fh.write(f"{gene}\t{chrom}\t{reason}\n")

    # --- finished genes, before their chromosome is packed ---------------------
    def work_path(self, chrom, gene):
        return self.path("_work", f"chr{chrom}", f"{gene}.npz")

    def has_work(self, chrom, gene):
        return os.path.exists(self.work_path(chrom, gene))

    def write_work(self, chrom, gene, probs, plan):
        path = self.work_path(chrom, gene)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path + ".tmp", "wb") as fh:
            np.savez_compressed(fh, probs=probs, plan=json.dumps(plan))
        replace(path + ".tmp", path)

    def read_work(self, chrom, gene):
        with np.load(self.work_path(chrom, gene)) as z:
            return z["probs"], json.loads(str(z["plan"]))

    # --- packed chromosomes ---------------------------------------------------------
    def packed(self, chrom):
        """The chromosome's sidecar (chr21.json), or None if it isn't packed."""
        path = self.path(f"chr{chrom}.json")
        return read_json(path) if os.path.exists(path) else None

    def pack(self, chrom, plans):
        """Writes chr<chrom>.h5 from the finished genes of `plans` (work files,
        plus any gene already in an earlier chr<chrom>.h5), then its sidecar,
        then deletes the work files."""
        import h5py

        precision = self.settings["precision"]
        h5_path = self.path(f"chr{chrom}.h5")
        old = h5py.File(h5_path, "r") if os.path.exists(h5_path) else None
        genes = {}
        try:
            with h5py.File(h5_path + ".tmp", "w") as out:
                out.attrs["build"] = self.settings["build"]
                out.attrs["precision"] = precision
                out.attrs["format"] = FORMAT_VERSION
                for name in sorted(set(plans) | set(old.keys() if old else ())):
                    if name in plans and self.has_work(chrom, name):
                        probs, plan = self.read_work(chrom, name)
                        group = out.create_group(name)
                        for i, track in enumerate(("acceptor", "donor")):
                            data = encode(probs[:, i], precision)
                            group.create_dataset(track, data=data, chunks=(min(len(data), 1 << 16),),
                                                 compression="gzip", compression_opts=4, shuffle=True)
                        group.attrs["plan"] = json.dumps(plan)
                    elif old is not None and name in old:
                        old.copy(old[name], out, name=name)
                        plan = json.loads(old[name].attrs["plan"])
                    else:
                        continue
                    genes[name] = plan
        finally:
            if old is not None:
                old.close()
        replace(h5_path + ".tmp", h5_path)
        write_json(self.path(f"chr{chrom}.json"), {
            "file": os.path.basename(h5_path), "size": os.path.getsize(h5_path),
            "sha256": sha256_of(h5_path), "precision": precision, "genes": genes,
        })
        for name in plans:
            if self.has_work(chrom, name):
                os.remove(self.work_path(chrom, name))

    def write_manifest(self, all_genes, extra=None):
        """manifest.json from build.json, the sidecars and the work folder.
        all_genes: names of every gene in the annotation, to state the scope."""
        files, genes = [], {}
        for side in sorted(glob.glob(self.path("chr*.json"))):
            sidecar = read_json(side)
            files.append({k: sidecar[k] for k in ("file", "size", "sha256")})
            for name, plan in sidecar["genes"].items():
                genes[name] = dict(plan, file=sidecar["file"])
        unpacked = sorted(os.path.basename(p)[:-4] for p in glob.glob(self.path("_work", "*", "*.npz")))
        skipped = []
        if os.path.exists(self.path("skipped.tsv")):
            with open(self.path("skipped.tsv"), encoding="utf-8") as fh:
                skipped = sorted({line.split("\t")[0] for line in fh.readlines()[1:]})
        done = set(genes) | set(unpacked) | set(skipped)
        # A validation or verification result stays as long as the packed files
        # it was made on do.
        old = read_json(self.path("manifest.json")) if os.path.exists(self.path("manifest.json")) else {}
        kept = {k: old[k] for k in ("validation", "local_verification") if old.get(k) and old.get("files") == files}
        manifest = dict(self.settings, **kept)
        manifest.update({
            "format": FORMAT_VERSION,
            "mask_rule": MASK_RULE,
            "vcf_comparable_rule": VCF_COMPARABLE_RULE,
            "stored_region": f"each gene's mask +/- {MARGIN} bases (one delta-score window)",
            "precision_note": "q001: integer thousandths (p = value * 0.001)" if self.settings["precision"] == "q001" else "",
            "scope": "all genes" if done >= set(all_genes) else f"partial: {len(done & set(all_genes))} of {len(all_genes)} genes",
            "gene_count": len(genes),
            "files": files,
            "genes": genes,
            "unpacked_genes": unpacked,
            "skipped_genes": skipped,
            "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
        })
        manifest.update(extra or {})
        write_json(self.path("manifest.json"), manifest)
        return manifest

    def read_gene(self, chrom, gene):
        """(float32 probabilities (n, 2), plan) of one gene, packed or not."""
        if self.has_work(chrom, gene):
            return self.read_work(chrom, gene)
        import h5py

        with h5py.File(self.path(f"chr{chrom}.h5"), "r") as f:
            group, precision = f[gene], f.attrs["precision"]
            probs = np.stack([decode(group[t][:], precision) for t in ("acceptor", "donor")], axis=1)
            return probs, json.loads(group.attrs["plan"])


# --- the app's side (spec 5.7, 5.8) -------------------------------------------------

def load_manifest(folder):
    """The folder's manifest.json; while a build is still on its first
    chromosome (no manifest yet), its build.json with no packed genes; None
    if the folder holds neither."""
    if os.path.exists(os.path.join(folder, "manifest.json")):
        return read_json(os.path.join(folder, "manifest.json"))
    if os.path.exists(os.path.join(folder, "build.json")):
        return dict(read_json(os.path.join(folder, "build.json")), genes={}, files=[], unpacked_genes=[],
                    skipped_genes=[], scope="partial: building", gene_count=0)
    return None


def check(folder, build, model_sha256, files=True):
    """The manifest of the database in `folder` if the app may use it for
    `build` with the installed models (model_sha256: predict.verify_models()),
    None if there is no database there. Otherwise raises DatabaseRefused with
    the reason. files: also check every packed file's size and SHA-256 against
    the manifest -- in full once; checked.json then remembers each file found
    intact by its size and modification time."""
    manifest = load_manifest(folder)
    if manifest is None:
        return None
    name = tx("The {build} reference score database in {folder}", build=manifest.get("build", ""), folder=folder)
    if manifest.get("format") != FORMAT_VERSION:
        raise DatabaseRefused(tx("{name} has format {format}, and this program reads format {ours}. It was made by "
                                 "another version of the program; build it again.", name=name,
                                 format=manifest.get("format"), ours=FORMAT_VERSION))
    if manifest.get("build") != build:
        raise DatabaseRefused(tx("{name} is for {other}, not {build}.", name=name, other=manifest.get("build"),
                                 build=build))
    if model_sha256 is not None and manifest.get("model_sha256") != model_sha256:
        differ = sorted(k for k in model_sha256 if (manifest.get("model_sha256") or {}).get(k) != model_sha256[k])
        raise DatabaseRefused(tx("{name} was computed with other SpliceAI model files ({files} differ from the "
                                 "installed ones), so its values wouldn't match this program's own scoring.",
                                 name=name, files=", ".join(differ)))
    for key in ("validation", "local_verification"):
        record = manifest.get(key) or {}
        if record.get("passed") is False:
            problem = record["problems"][0] if record.get("problems") else ""
            what = (tx("{name} failed its validation on {date}", name=name, date=record.get("date"))
                    if key == "validation" else
                    tx("{name} failed its verification on {date}", name=name, date=record.get("date")))
            raise DatabaseRefused(tx("{what} ({problem}). Build it again, or import a database that passed.",
                                     what=what, problem=problem) if problem else
                                  tx("{what}. Build it again, or import a database that passed.", what=what))
    if files:
        _check_files(folder, manifest, name)
    return manifest


def _check_files(folder, manifest, name):
    stamp_path = os.path.join(folder, "checked.json")
    try:
        stamps = read_json(stamp_path)
    except (OSError, ValueError):
        stamps = {}
    changed = False
    for entry in manifest["files"]:
        path = os.path.join(folder, entry["file"])
        if not os.path.exists(path):
            raise DatabaseRefused(tx("{name} is incomplete: {file} is missing.", name=name, file=entry["file"]))
        st = os.stat(path)
        stamp = [st.st_size, st.st_mtime_ns, entry["sha256"]]
        if stamps.get(entry["file"]) == stamp:
            continue
        if st.st_size != entry["size"] or sha256_of(path) != entry["sha256"]:
            raise DatabaseRefused(tx("{name} is damaged: {file} isn't the file its manifest describes (size or "
                                     "SHA-256 differ). Build it again, or import an intact copy.", name=name,
                                     file=entry["file"]))
        stamps[entry["file"]] = stamp
        changed = True
    if changed:
        try:
            write_json(stamp_path, stamps)
        except OSError:     # a read-only copy: it is simply checked in full next time
            pass


class Reader:
    """One database folder: its packed chromosomes, and genes a running build
    has finished but not packed yet. A file is open only for the moment of a
    read, so a build can replace it meanwhile; refresh() rereads the manifest."""

    def __init__(self, folder, manifest=None):
        self.folder = folder
        self.manifest = manifest or load_manifest(folder)
        self.build = self.manifest["build"]
        self.precision = self.manifest["precision"]

    def refresh(self):
        self.manifest = load_manifest(self.folder) or self.manifest

    def _loose(self, gene):
        hits = glob.glob(os.path.join(glob.escape(self.folder), "_work", "chr*", glob.escape(gene) + ".npz"))
        return hits[0] if hits else None

    def plan(self, gene):
        """What the database stores for `gene` (plan_gene()'s dict), or None."""
        if gene in self.manifest["genes"]:
            return self.manifest["genes"][gene]
        path = self._loose(gene)
        if path is None:
            return None
        with np.load(path) as z:
            return json.loads(str(z["plan"]))

    def values(self, gene, start=None, end=None):
        """float32 array (n, 2), acceptor and donor, of `gene` at start..end
        (1-based, inclusive, genomic order; default: everything it stores).
        Raises KeyError for a gene the database doesn't have."""
        plan = self.plan(gene)
        if plan is None:
            raise KeyError(gene)
        start = plan["start"] if start is None else start
        end = plan["end"] if end is None else end
        if not plan["start"] <= start <= end <= plan["end"]:
            raise ValueError(f"{gene} is stored at {plan['start']:,}..{plan['end']:,}, not {start:,}..{end:,}")
        a, b = start - plan["start"], end - plan["start"] + 1
        if gene not in self.manifest["genes"]:
            with np.load(self._loose(gene)) as z:
                # Rounded the way packing rounds them, so a gene reads the same before and after.
                return decode(encode(z["probs"][a:b], self.precision), self.precision)
        import h5py

        with h5py.File(os.path.join(self.folder, plan["file"]), "r") as f:
            group = f[gene]
            return np.stack([decode(group[t][a:b], self.precision) for t in ("acceptor", "donor")], axis=1)

    def genes_at(self, chrom, start, end, strand=None):
        """Names of the genes whose stored region overlaps chrom:start..end
        ("21" or "chr21"), only those on `strand` if it is given."""
        chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
        plans = [p for p in self.manifest["genes"].values() if p["chrom"] == chrom]
        for path in glob.glob(os.path.join(glob.escape(self.folder), "_work", f"chr{chrom}", "*.npz")):
            if os.path.basename(path)[:-4] not in self.manifest["genes"]:
                with np.load(path) as z:
                    plans.append(json.loads(str(z["plan"])))
        return sorted(p["gene"] for p in plans
                      if p["start"] <= end and p["end"] >= start and strand in (None, p["strand"]))


def open_reader(folder, build, model_sha256, files=True):
    """A Reader for the database in `folder` (check()'s rules), or None."""
    manifest = check(folder, build, model_sha256, files=files)
    return Reader(folder, manifest) if manifest is not None else None


class ReferenceScores:
    """Reference probabilities of whole genes for the app (spec 5.7): from the
    database where it has the gene; otherwise predicted live the database's
    way (predict.py, the same mask and stored region) and kept in a cache
    folder, so a gene is computed at most once. TensorFlow and the gene
    tables are loaded only when a gene has to be computed."""

    def __init__(self, build, fasta_path, reader=None, cache_dir=None, refseq_dir=None, mane_dir=None):
        self.build, self.fasta_path, self.reader = build, fasta_path, reader
        self.cache_dir = cache_dir or os.path.join(default_root(), "live-cache", build)
        self.refseq_dir, self.mane_dir = refseq_dir, mane_dir
        self.precision = reader.precision if reader else DEFAULT_PRECISION
        self._genes = self._predictor = self._settings = None

    def gene(self, name, checkpoint=None, on_progress=None):
        """(values, plan, source) of a gene by SpliceAI's or RefSeq's symbol:
        values as Reader.values(); source "database" or "computed live" (now,
        or earlier and read from the cache). checkpoint() is called between
        model calls (pause / cancel), on_progress(fraction) after each while it
        is computed live. KeyError for a gene SpliceAI doesn't know."""
        if self.reader is not None:
            plan = self.reader.plan(name)
            if plan is not None:
                return self.reader.values(name), plan, "database"
        path = os.path.join(self.cache_dir, f"{name}.npz")
        if os.path.exists(path):
            with np.load(path) as z:
                if json.loads(str(z["settings"])) == self.settings():
                    precision = str(z["precision"])
                    return decode(z["probs"], precision), json.loads(str(z["plan"])), "computed live"
        values, plan = self._compute(name, checkpoint, on_progress)
        return values, plan, "computed live"

    def settings(self):
        """What a cached gene depends on besides the genome: models, annotation, precision."""
        if self._settings is None:
            from . import predict, transcripts
            from .mane import find_cached_summary

            paths = (transcripts.annotation_path(self.build), transcripts.table_path(self.build, self.refseq_dir),
                     find_cached_summary(self.mane_dir))
            self._settings = {"model_sha256": predict.verify_models(),
                              "annotation": {os.path.basename(p): sha256_of(p) for p in paths if p},
                              "precision": self.precision}
        return self._settings

    def resolve(self, name):
        """The transcripts.Gene a symbol (SpliceAI's or its RefSeq one) names; KeyError if none."""
        if self._genes is None:
            from . import transcripts
            from .mane import find_cached_summary

            genes = transcripts.load_genes(self.build, self.refseq_dir, find_cached_summary(self.mane_dir))
            self._genes = {g.name: g for g in genes}
            for g in genes:
                if g.canonical:
                    self._genes.setdefault(g.canonical.symbol, g)
        if name not in self._genes:
            raise KeyError(f"{name} is not in SpliceAI's {self.build} gene table, under either symbol")
        return self._genes[name]

    def predictor(self):
        """The models (predict.Predictor), loaded once; the patient's prediction uses them too."""
        if self._predictor is None:
            from . import predict

            self._predictor = predict.Predictor(self.fasta_path, self.build)
        return self._predictor

    def _compute(self, name, checkpoint, on_progress=None):
        from .normalize import resolve_contig
        from .predict import DEFAULT_CHUNK

        gene = self.resolve(name)
        fasta = self.predictor().annotator.ref_fasta
        record = fasta[resolve_contig(gene.chrom, set(fasta.keys()))]
        plan = plan_gene(gene, len(record))
        chunks, done = -(-(plan["end"] - plan["start"] + 1) // DEFAULT_CHUNK), [0]

        def step():
            done[0] += 1
            if on_progress:
                on_progress(done[0] / chunks)
            if checkpoint:
                checkpoint()

        probs = self._predictor.predict(lambda a, b: record[a - 1:b].seq, len(record), gene.strand, *plan["mask"],
                                        plan["start"], plan["end"], checkpoint=step)
        stored = encode(probs, self.precision)
        os.makedirs(self.cache_dir, exist_ok=True)
        path = os.path.join(self.cache_dir, f"{name}.npz")
        with open(path + ".tmp", "wb") as fh:
            np.savez_compressed(fh, probs=stored, plan=json.dumps(plan), precision=self.precision,
                                settings=json.dumps(self.settings()))
        replace(path + ".tmp", path)
        return decode(stored, self.precision), plan


def verify_against_set(folder, set_path, tolerance=VERIFY_TOLERANCE):
    """Checks a locally built database against a validation set (spec 5.8):
    every gene in both has to store the same region and the same values within
    `tolerance`. Returns the record kept in the manifest as
    "local_verification"; its "passed" is None when no gene is in both."""
    import h5py

    reader = Reader(folder)
    problems, worst, worst_gene, genes = [], 0.0, None, 0
    with h5py.File(set_path, "r") as vs:
        if vs.attrs["build"] != reader.build:
            problems.append(f"the validation set is for {vs.attrs['build']}, the database for {reader.build}")
        elif json.loads(vs.attrs["model_sha256"]) != reader.manifest.get("model_sha256"):
            problems.append("the validation set was made with other SpliceAI model files")
        else:
            for name in sorted(vs.keys()):
                plan = reader.plan(name)
                if plan is None:
                    continue
                expected = json.loads(vs[name].attrs["plan"])
                if (plan["start"], plan["end"], plan["mask"]) != (expected["start"], expected["end"], expected["mask"]):
                    problems.append(f"{name} is stored at {plan['start']:,}..{plan['end']:,} with mask {plan['mask']}; "
                                    f"the validation set has {expected['start']:,}..{expected['end']:,} with mask "
                                    f"{expected['mask']} (another gene annotation?)")
                    continue
                ref = np.stack([decode(vs[name][t][:], vs.attrs["precision"]) for t in ("acceptor", "donor")], axis=1)
                diff = float(np.abs(reader.values(name) - ref).max())
                genes += 1
                if diff > worst:
                    worst, worst_gene = diff, name
                if diff > tolerance:
                    problems.append(f"{name}'s values differ by up to {diff:.6f} (allowed: {tolerance:.6f})")
    return {"passed": False if problems else (True if genes else None),
            "date": time.strftime("%Y-%m-%d %H:%M:%S"), "genes": genes, "max_abs_diff": worst,
            "worst_gene": worst_gene, "tolerance": tolerance, "problems": problems[:20],
            "validation_set": {"file": os.path.basename(set_path), "sha256": sha256_of(set_path)}}


def record_verification(folder, record):
    """Writes a verify_against_set() result into the folder's manifest."""
    path = os.path.join(folder, "manifest.json")
    write_json(path, dict(read_json(path), local_verification=record))


def import_database(src, model_sha256, root=None, on_progress=None):
    """Copies a database obtained separately (a folder, e.g. on a USB drive)
    to root/<build> (default_root()): checks it as check() does, copies its
    manifest and packed files to <build>.importing, checks the copies'
    SHA-256, and only then replaces the database of that build, if there is
    one. on_progress(bytes_done, bytes_total). Returns the new folder; raises
    DatabaseRefused if `src` can't be used."""
    if not os.path.exists(os.path.join(src, "manifest.json")):
        raise DatabaseRefused(tx("{folder} holds no finished reference score database (no manifest.json).",
                                 folder=src))
    manifest = read_json(os.path.join(src, "manifest.json"))
    build = manifest.get("build")
    if build not in ("hg19", "hg38"):
        raise DatabaseRefused(tx("{folder} holds a database for an unknown genome build ({build}).", folder=src,
                                 build=build))
    check(src, build, model_sha256, files=False)
    root = root or default_root()
    dest = os.path.join(root, build)
    if os.path.normcase(os.path.abspath(src)) == os.path.normcase(os.path.abspath(dest)):
        raise DatabaseRefused(tx("{folder} is already the {build} database this program uses.", folder=src,
                                 build=build))
    names = ["manifest.json"] + [n for n in ("build.json", "skipped.tsv") if os.path.exists(os.path.join(src, n))]
    for entry in manifest["files"]:
        names += [entry["file"], entry["file"][:-3] + ".json"]
    missing = [n for n in names if not os.path.exists(os.path.join(src, n))]
    if missing:
        raise DatabaseRefused(tx("{folder} is incomplete: {files} missing.", folder=src, files=", ".join(missing)))
    total, done = sum(os.path.getsize(os.path.join(src, n)) for n in names), 0
    tmp = dest + ".importing"
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    try:
        for n in names:
            with open(os.path.join(src, n), "rb") as fin, open(os.path.join(tmp, n), "wb") as fout:
                for block in iter(lambda: fin.read(1 << 22), b""):
                    fout.write(block)
                    done += len(block)
                    if on_progress:
                        on_progress(done, total)
        check(tmp, build, model_sha256, files=True)
        old = dest + ".replaced"
        if os.path.exists(dest):
            shutil.rmtree(old, ignore_errors=True)
            os.replace(dest, old)
        os.replace(tmp, dest)
        shutil.rmtree(old, ignore_errors=True)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return dest
