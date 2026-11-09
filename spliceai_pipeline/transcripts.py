"""Curated RefSeq transcripts, and which of them belong to each SpliceAI gene
-- the canonical one and the per-gene mask of the reference score database
(spec 5.2, 5.3).

Transcripts come from UCSC's ncbiRefSeqCurated table of each build (fetched by
installer/fetch-thirdparty.ps1, pinned by SHA-256): exons, CDS, every curated
NM_/NR_ transcript and never a predicted XM_. MANE Select (mane.py) names the
canonical one.

SpliceAI's gene table uses older symbols than RefSeq does today (its FAM132A
is C1QTNF12), so a SpliceAI gene is matched to transcripts by locus -- same
chromosome and strand, overlapping -- rather than through an alias list:

1. MANE Select: the MANE Select transcript at the locus whose MANE symbol is
   the SpliceAI gene's; else the one sharing the most exonic bases with
   SpliceAI's own exons for the gene. When several MANE genes share exonic
   bases with it (read-through and nested genes), the gene is flagged
   ambiguous and the others are listed.
2. Otherwise the longest-CDS NM_ at the locus, preferring the gene's own
   symbol -- the VCF tab's curated-NM_ preference, limited to transcripts
   overlapping the gene on its strand. Under another symbol, a transcript
   must share exonic bases with SpliceAI's gene, as in 1.
3. Otherwise none: the gene stays in the database with SpliceAI's interval,
   but no reference cDNA can be built for it.

The mask (the stretch whose bases are real sequence, not N) is SpliceAI's
interval joined with the canonical transcript and every NM_ of the same
symbol overlapping either. An NM_ sharing exons with another gene's MANE
Select transcript (a read-through, or a gene overlapping another) is left out
and listed instead.
"""
import bisect
import csv
import gzip
import importlib.util
import os
import sys
from dataclasses import dataclass, field

from .app_paths import app_root_dir
from .texts import tx

# score.BUILD_TO_ANNOTATION, repeated so this module loads without TensorFlow.
BUILD_TO_ANNOTATION = {"hg19": "grch37", "hg38": "grch38"}

DEFAULT_REFSEQ_DIR = os.path.join(app_root_dir(), "refseq")
TABLE_NAME = "{build}.ncbiRefSeqCurated.txt.gz"


def default_refseq_dir():
    """The folder with the curated RefSeq tables of both builds: "refseq" next
    to the program (in a source checkout: the repository's own), inside the
    frozen bundle (the .spec puts them there), or installer/thirdparty/refseq
    (fetch-thirdparty.ps1).
    DEFAULT_REFSEQ_DIR when none has them, so the error names it."""
    root = app_root_dir()
    for folder in (DEFAULT_REFSEQ_DIR, os.path.join(getattr(sys, "_MEIPASS", root), "refseq"),
                   os.path.join(root, "installer", "thirdparty", "refseq")):
        if all(os.path.isfile(os.path.join(folder, TABLE_NAME.format(build=b))) for b in BUILD_TO_ANNOTATION):
            return folder
    return DEFAULT_REFSEQ_DIR


def missing_table(build, refseq_dir=None):
    """None, or what to tell the user when the build's RefSeq table isn't there."""
    path = table_path(build, refseq_dir)
    if os.path.isfile(path):
        return None
    return tx("The curated RefSeq transcript table for {build} ({name}) is missing from {folder}. The Sequence to "
              "cDNA tab needs it to find genes and their transcripts. Reinstall the program; when running from "
              "source, they are in the repository's refseq folder, and installer\\fetch-thirdparty.ps1 fetches them "
              "too.", build=build, name=os.path.basename(path), folder=os.path.dirname(path))


PRIMARY_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}

RULE_MANE = "MANE Select"
RULE_MANE_VERSION = "MANE Select, other version"
RULE_LONGEST_CDS = "longest CDS (no MANE Select at this locus)"
RULE_NONE = "none (no NM_ at this locus)"


@dataclass(frozen=True)
class Transcript:
    accession: str          # with version, e.g. NM_000059.4
    symbol: str             # RefSeq gene symbol (UCSC name2)
    chrom: str              # UCSC name, e.g. chr13
    strand: str
    start: int              # 1-based, inclusive: first base of the lowest exon
    end: int                # last base of the highest exon
    cds: tuple              # (first, last) CDS base, 1-based, genomic order; None if non-coding
    exons: tuple            # ((start, end), ...) 1-based inclusive, ascending

    @property
    def base_accession(self):
        return self.accession.split(".")[0]

    @property
    def cds_length(self):
        if self.cds is None:
            return 0
        return sum(max(0, min(e, self.cds[1]) - max(s, self.cds[0]) + 1) for s, e in self.exons)


@dataclass
class Gene:
    """One gene of SpliceAI's annotation, with its transcripts and mask."""
    name: str
    chrom: str              # as in SpliceAI's table, e.g. "21"
    strand: str
    start: int              # SpliceAI's interval: TX_START+1 .. TX_END (1-based)
    end: int
    exons: tuple            # SpliceAI's exons, 1-based inclusive
    canonical: Transcript = None
    rule: str = RULE_NONE
    transcripts: list = field(default_factory=list)     # canonical first, then alternatives
    excluded: list = field(default_factory=list)        # (accession, why)
    other_mane: list = field(default_factory=list)      # other MANE genes at an ambiguous locus

    @property
    def mask(self):
        """(start, end) of the bases that are real sequence for this gene."""
        return (min([self.start] + [t.start for t in self.transcripts]),
                max([self.end] + [t.end for t in self.transcripts]))


def table_path(build, refseq_dir=None):
    return os.path.join(refseq_dir or default_refseq_dir(), TABLE_NAME.format(build=build))


def spliceai_package_dir():
    """The installed spliceai package's folder, found without importing it
    (its __init__ calls signal.signal(), which fails outside the main thread;
    see score.py)."""
    spec = importlib.util.find_spec("spliceai")
    if spec is None or not spec.submodule_search_locations:
        raise ModuleNotFoundError("SpliceAI isn't installed.", name="spliceai")
    return spec.submodule_search_locations[0]


def annotation_path(build):
    """SpliceAI's own gene table for a build (inside the spliceai package)."""
    return os.path.join(spliceai_package_dir(), "annotations", f"{BUILD_TO_ANNOTATION[build]}.txt")


def read_table(path):
    """Curated NM_/NR_ transcripts on chr1-22, X, Y from a UCSC genePred table
    (0-based half-open, converted to 1-based inclusive here)."""
    transcripts = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for f in csv.reader(fh, delimiter="\t"):
            if f[2] not in PRIMARY_CHROMS or not f[1].startswith(("NM_", "NR_")):
                continue
            cds_start, cds_end = int(f[6]), int(f[7])
            starts = [int(x) for x in f[9].rstrip(",").split(",")]
            ends = [int(x) for x in f[10].rstrip(",").split(",")]
            transcripts.append(Transcript(
                accession=f[1], symbol=f[12], chrom=f[2], strand=f[3],
                start=int(f[4]) + 1, end=int(f[5]),
                cds=(cds_start + 1, cds_end) if cds_end > cds_start else None,
                exons=tuple((s + 1, e) for s, e in zip(starts, ends)),
            ))
    return transcripts


def read_spliceai_genes(path):
    """SpliceAI's gene table (TX_START and EXON_START are 0-based there)."""
    genes = []
    with open(path, encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            starts = [int(x) + 1 for x in r["EXON_START"].rstrip(",").split(",")]
            ends = [int(x) for x in r["EXON_END"].rstrip(",").split(",")]
            genes.append(Gene(name=r["#NAME"], chrom=r["CHROM"], strand=r["STRAND"],
                              start=int(r["TX_START"]) + 1, end=int(r["TX_END"]),
                              exons=tuple(zip(starts, ends))))
    return genes


def bare_chrom(chrom):
    """"chr7" and "7" alike, as SpliceAI's Annotator compares them."""
    c = str(chrom)
    return c[3:] if c.lower().startswith("chr") else c


class GeneRegions:
    """SpliceAI's intervals of a chosen set of genes in one build: the stretches
    where its Annotator reports one of them for a variant (by the variant's
    position, TX_START+1 .. TX_END). Lets a VCF run score only those genes."""

    def __init__(self, build, genes):
        by_upper = {}
        wanted = {g.strip().upper() for g in genes if g and g.strip()}
        spans = {}
        for g in read_spliceai_genes(annotation_path(build)):
            if g.name.upper() in wanted:
                by_upper[g.name.upper()] = g.name
                spans.setdefault(bare_chrom(g.chrom), []).append((g.start, g.end))
        self.genes = set(by_upper.values())          # as SpliceAI spells them
        self.unknown = sorted(wanted - set(by_upper))  # not in this build's table
        self._starts, self._ends = {}, {}
        for chrom, items in spans.items():
            merged = []
            for s, e in sorted(items):
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self._starts[chrom] = [s for s, _ in merged]
            self._ends[chrom] = [e for _, e in merged]

    def contains(self, chrom, pos, margin=0):
        """Whether pos lies in one of the genes' intervals, widened by margin."""
        starts = self._starts.get(bare_chrom(chrom))
        if not starts:
            return False
        i = bisect.bisect_right(starts, pos + margin) - 1
        return i >= 0 and self._ends[bare_chrom(chrom)][i] >= pos - margin


def exonic_overlap(a, b):
    """Bases two ascending exon lists have in common."""
    i = j = shared = 0
    while i < len(a) and j < len(b):
        shared += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]) + 1)
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return shared


class _Index:
    """Transcripts by (chrom, strand), for overlap queries."""

    def __init__(self, transcripts):
        self._by_key = {}
        for t in sorted(transcripts, key=lambda t: t.start):
            self._by_key.setdefault((t.chrom, t.strand), []).append(t)
        self._starts = {k: [t.start for t in v] for k, v in self._by_key.items()}
        self._longest = {k: max(t.end - t.start for t in v) for k, v in self._by_key.items()}

    def overlapping(self, chrom, strand, a, b):
        key = (chrom, strand)
        if key not in self._by_key:
            return []
        starts, items = self._starts[key], self._by_key[key]
        lo = bisect.bisect_left(starts, a - self._longest[key])
        hi = bisect.bisect_right(starts, b)
        return [t for t in items[lo:hi] if t.end >= a]


def assign(genes, transcripts, mane_select):
    """Fills in each Gene's canonical transcript, rule, transcripts and
    exclusions (see the module docstring). mane_select is
    mane.load_mane_select()'s {base accession: (accession, symbol)}."""
    index = _Index(transcripts)

    def mane_symbol(t):
        hit = mane_select.get(t.base_accession)
        return hit[1] if hit else None

    def shares_exons_with_other_mane(t):
        for m in index.overlapping(t.chrom, t.strand, t.start, t.end):
            if mane_symbol(m) is not None and m.symbol != t.symbol and exonic_overlap(m.exons, t.exons):
                return f"shares exons with {m.symbol}'s MANE Select transcript {m.accession}"
        return None

    for g in genes:
        chrom = "chr" + g.chrom
        here = index.overlapping(chrom, g.strand, g.start, g.end)
        mane_here = [t for t in here if mane_symbol(t) is not None]
        canonical = None
        named = [t for t in mane_here if mane_symbol(t) == g.name]
        if named:
            canonical = named[0]
        else:
            shared = [(exonic_overlap(t.exons, g.exons), t) for t in mane_here]
            shared = [(n, t) for n, t in shared if n > 0]
            if shared:
                canonical = max(shared, key=lambda p: p[0])[1]
                g.other_mane = sorted({mane_symbol(t) for _, t in shared} - {mane_symbol(canonical)})
        if canonical is not None:
            g.rule = RULE_MANE if canonical.accession == mane_select[canonical.base_accession][0] else RULE_MANE_VERSION
        else:
            # Under another symbol only when it shares exons with SpliceAI's
            # gene: a host gene's transcript isn't that of a gene nested in its
            # intron.
            nm = [t for t in here if t.accession.startswith("NM_")
                  and (t.symbol == g.name or exonic_overlap(t.exons, g.exons))]
            nm = [t for t in nm if not shares_exons_with_other_mane(t)]
            nm = [t for t in nm if t.symbol == g.name] or nm
            if nm:
                canonical = max(nm, key=lambda t: (t.cds_length, t.end - t.start, t.accession))
                g.rule = RULE_LONGEST_CDS
        if canonical is None:
            continue
        g.canonical = canonical
        g.transcripts = [canonical]
        span = (min(g.start, canonical.start), max(g.end, canonical.end))
        for t in index.overlapping(chrom, g.strand, *span):
            if t is canonical or not t.accession.startswith("NM_") or t.symbol != canonical.symbol:
                continue
            why = shares_exons_with_other_mane(t)
            if why:
                g.excluded.append((t.accession, why))
            else:
                g.transcripts.append(t)
    return genes


def load_genes(build, refseq_dir=None, mane_path=None):
    """SpliceAI's genes for a build, each assigned its transcripts and mask."""
    from .mane import find_cached_summary, load_mane_select

    mane_path = mane_path or find_cached_summary()
    if not mane_path:
        raise FileNotFoundError("No MANE Select summary file found (see spliceai_pipeline/mane.py).")
    if missing_table(build, refseq_dir):
        raise FileNotFoundError(missing_table(build, refseq_dir))
    return assign(read_spliceai_genes(annotation_path(build)), read_table(table_path(build, refseq_dir)),
                  load_mane_select(mane_path))
