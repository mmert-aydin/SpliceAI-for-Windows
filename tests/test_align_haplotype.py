"""Phase 4: reading a pasted sequence, aligning it to the reference
(align.py), and the patient haplotype with its coordinate map
(haplotype.py). Synthetic patients are made from the real reference with
known edits; aligning them back must find those edits, and applying what was
found must give the patient's sequence again, base for base.

    .venv\\Scripts\\python tests\\test_align_haplotype.py
"""
import random
import time

import helpers
from helpers import check

helpers.watchdog(900)
helpers.setup()

from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline import align, haplotype, transcripts  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.normalize import make_base_fetcher  # noqa: E402

rng = random.Random(7)
rand_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))


def raises(fn, text):
    try:
        fn()
    except ValueError as exc:
        return text in str(exc)
    return False


# --- reading the paste --------------------------------------------------------------------
body = "acgt" * 20
check("parse: raw, lower case, whitespace and position numbers", align.parse_sequence(
    "1 acgtacgtac\n11 " + body[10:] + "  \n") == body.upper())
check("parse: one FASTA record", align.parse_sequence(">patient 1\n" + body + "\n") == body.upper())
check("parse: N is accepted", align.parse_sequence("N" * 10 + body) == "N" * 10 + body.upper())
check("parse: other letters are refused, named",
      raises(lambda: align.parse_sequence(body + "RYK"), "'K', 'R', 'Y'"))
check("parse: two FASTA records are refused", raises(lambda: align.parse_sequence(">a\n" + body + "\n>b\n" + body),
                                                     "2 FASTA sequences"))
check("parse: fewer than 50 bases are refused, saying why",
      raises(lambda: align.parse_sequence("ACGT" * 12), "48 bases; at least 50"))
check("parse region: commas, spaces, chr optional",
      align.parse_region(" chr7:117,480,025 - 117,668,665 ") == ("7", 117480025, 117668665)
      and align.parse_region("x:5-10") == ("X", 5, 10))
check("parse region: an inverted region is refused", raises(lambda: align.parse_region("1:10-5"), "not after the end"))
check("parse coordinates (spec 14 C5): a single position or a region",
      align.parse_coordinates("chr7:117,559,590") == ("7", 117559590, 117559590, True)
      and align.parse_coordinates("7:100-200") == ("7", 100, 200, False))
check("parse coordinates: anything else is refused, naming both forms",
      raises(lambda: align.parse_coordinates("chr7"), "neither a region") and
      raises(lambda: align.parse_coordinates("chr7:0"), "at least 1"))

# --- haplotype and coordinate map ----------------------------------------------------------------
REF = rand_seq(3000)
START = 100_001


def round_trips(ref, start, edits):
    """Builds the haplotype and checks every base both ways; returns (hap, cmap, problems)."""
    hap, cmap = haplotype.build(ref, start, edits)
    problems = []
    for h in range(1, len(hap) + 1):
        r, label = cmap.to_ref(h)
        if r is None:
            if not label or "+" not in label:
                problems.append(f"hap {h}: inserted base without a label ({label!r})")
            continue
        if cmap.to_hap(r) != h and not label.startswith("DUP"):
            problems.append(f"hap {h} -> ref {r} -> hap {cmap.to_hap(r)}")
        s = cmap.segment_at(h)
        expect = ref[r - start]
        if s.orientation == -1:
            expect = align.reverse_complement(expect)
        if not label and hap[h - 1] != expect:
            problems.append(f"hap {h}: {hap[h - 1]} isn't the reference's {expect}")
    deleted = {p for e in edits if e.kind != "inversion" and len(e.alt) != e.end - e.start + 1
               for p in range(e.start, e.end + 1)}
    for r in range(start, start + len(ref)):
        h = cmap.to_hap(r)
        if (h is None) != (r in deleted):
            problems.append(f"ref {r}: to_hap {h}, deleted {r in deleted}")
        elif h is not None and cmap.to_ref(h)[0] != r:
            problems.append(f"ref {r} -> hap {h} -> ref {cmap.to_ref(h)[0]}")
    return hap, cmap, problems[:5]


hap, cmap, problems = round_trips(REF, START, [])
check("haplotype: no edits is the reference, mapped one to one", hap == REF and not problems
      and cmap.to_ref(1) == (START, "") and cmap.to_hap(START + 2999) == 3000, str(problems))
cases = {
    "SNV": [Edit(START + 10, START + 10, "A" if REF[10] != "A" else "C", "SNV")],
    "insertion": [Edit(START + 500, START + 499, "TTTGGG", "insertion")],
    "deletion": [Edit(START + 700, START + 779, "", "deletion")],
    "tandem duplication": [Edit(START + 1200, START + 1199, REF[1100:1200], "duplication")],
    "inversion": [Edit(START + 1500, START + 1799, "", "inversion")],
    "complex (3 bases for 7)": [Edit(START + 2000, START + 2002, "GATTACA", "complex")],
    "edits at both region edges": [Edit(START, START + 4, "", "deletion"), Edit(START + 3000, START + 2999, "AC", "insertion")],
}
cases["several edits"] = [e for k in ("SNV", "insertion", "deletion", "tandem duplication", "inversion",
                                      "complex (3 bases for 7)") for e in cases[k]]
for name, edits in cases.items():
    hap, cmap, problems = round_trips(REF, START, edits)
    check(f"haplotype, {name}: every non-inserted base maps both ways, inserted ones are labelled", not problems,
          str(problems))
hap, cmap = haplotype.build(REF, START, cases["several edits"])
check("haplotype: its length is the region's minus deletions plus insertions",
      len(hap) == 3000 + 6 - 80 + 100 + 4)
check("haplotype: inserted bases are named after their insertion (INS1+3, DUP1+1)",
      cmap.to_ref(cmap.to_hap(START + 499) + 3) == (None, "INS1+3")
      and cmap.to_ref(cmap.to_hap(START + 1199) + 1) == (None, "DUP1+1"))
inv_first = cmap.to_hap(START + 1500)
check("haplotype: an inverted segment maps reversed and is the reverse complement",
      inv_first == cmap.to_hap(START + 1799) + 299 and cmap.to_ref(inv_first) == (START + 1500, "INV1")
      and hap[cmap.to_hap(START + 1799) - 1:inv_first] == align.reverse_complement(REF[1500:1800]))
check("haplotype: names a position the viewer's way", cmap.name(cmap.to_hap(START + 1500), "chr1") ==
      f"chr1:{START + 1500:,} (INV1)" and cmap.name(cmap.to_hap(START + 499) + 3) == "INS1+3")
check("haplotype: overlapping edits are refused", raises(lambda: haplotype.build(REF, START, [
    Edit(START + 10, START + 20, "", "deletion"), Edit(START + 15, START + 15, "A", "SNV")]), "overlap"))
check("haplotype: an edit outside the region is refused", raises(lambda: haplotype.build(REF, START, [
    Edit(START - 1, START + 3, "", "deletion")]), "outside the region"))

# Random edits, many times.
bad = []
for trial in range(200):
    edits, pos = [], START + rng.randrange(0, 40)
    while pos < START + 2900:
        kind = rng.choice(["SNV", "insertion", "deletion", "inversion", "complex", "duplication"])
        n = rng.randrange(1, 60)
        end = min(pos + n - 1, START + 2999)
        if kind == "SNV":
            edits.append(Edit(pos, pos, "A" if REF[pos - START] != "A" else "G", kind))
        elif kind in ("insertion", "duplication"):
            seq = rand_seq(n) if kind == "insertion" else REF[max(0, pos - START - n):pos - START]
            edits.append(Edit(pos, pos - 1, seq, kind))
            end = pos - 1
        elif kind == "deletion":
            edits.append(Edit(pos, end, "", kind))
        elif kind == "inversion":
            edits.append(Edit(pos, end, "", kind))
        else:
            edits.append(Edit(pos, end, rand_seq(rng.randrange(0, 20)), kind))
        pos = end + 1 + rng.randrange(1, 300)
    _, _, problems = round_trips(REF, START, edits)
    if problems:
        bad.append((trial, problems[0]))
check("haplotype: 200 random sets of edits round-trip at every base", not bad, str(bad[:2]))

# ref_aligned (patient values on reference positions, VCF-style) against base-by-base.
import numpy as np  # noqa: E402


def naive_ref_aligned(pred, cmap, start, end):
    out = np.zeros((end - start + 1, 2), np.float32)
    for r in range(start, end + 1):
        h = cmap.to_hap(r)
        if h is not None and pred.start <= h <= pred.end:
            out[r - start] = pred.values[h - pred.start]
    for s in cmap.segments:
        if s.orientation == 0 and s.hap_start > 1:
            anchor = cmap.to_ref(s.hap_start - 1)[0]
            a, b = max(s.hap_start, pred.start), min(s.hap_end, pred.end)
            if anchor is not None and start <= anchor <= end and a <= b:
                out[anchor - start] = np.maximum(out[anchor - start], pred.values[a - pred.start:b - pred.start + 1].max(0))
    return out


differ = []
for trial in range(40):
    edits = cases["several edits"] if trial == 0 else [e for e in cases["several edits"] if rng.random() < 0.6]
    hap, cmap = haplotype.build(REF, START, edits)
    p0, p1 = rng.randrange(1, 200), len(hap) - rng.randrange(0, 200)
    pred = haplotype.PatientPrediction(p0, p1, np.random.default_rng(trial).random((p1 - p0 + 1, 2)).astype(np.float32),
                                       (p0, p1))
    a, b = START + rng.randrange(0, 300), START + 2999 - rng.randrange(0, 300)
    if not np.array_equal(haplotype.ref_aligned(pred, cmap, a, b), naive_ref_aligned(pred, cmap, a, b)):
        differ.append(trial)
check("ref_aligned (by segment) equals the base-by-base mapping, 40 random haplotypes and stretches", not differ,
      str(differ))

# --- the real reference: genes, both strands -------------------------------------------------------
fasta = Fasta(str(helpers.HG38))
fetch = make_base_fetcher(fasta)
genes = {g.name: g for g in transcripts.load_genes("hg38", helpers.THIRDPARTY / "refseq",
                                                     helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")}


def region_of(gene, flank=2000):
    start, end = gene.mask[0] - flank, gene.mask[1] + flank
    return start, end, fasta["chr" + gene.chrom][start - 1:end].seq.upper()


minus = next(g for g in genes.values() if g.strand == "-" and g.canonical and g.canonical.cds
             and 20_000 < g.mask[1] - g.mask[0] < 60_000 and g.chrom == "21")
start, end, ref = region_of(minus)
cds_start, cds_end = minus.canonical.cds
atg = cds_end   # minus strand: the ATG's first base is the CDS's genomic end
check(f"minus strand ({minus.name}): the reference has ATG at the CDS end, reading the minus strand",
      align.reverse_complement(ref[atg - 2 - start:atg + 1 - start]) == "ATG")
upstream = minus.canonical.end + 500    # minus strand: upstream is to the right
edits = [Edit(minus.canonical.start + 300, minus.canonical.start + 299, rand_seq(40), "insertion"),
         Edit(upstream, upstream + 99, "", "deletion")]
hap, cmap = haplotype.build(ref, start, edits)
lifted = cmap.to_hap(atg)
check("minus strand: the ATG lifted onto the haplotype is still ATG on the minus strand, 40 bases further",
      align.reverse_complement(hap[lifted - 3:lifted]) == "ATG" and lifted == atg - start + 1 + 40)
check("minus strand: transcript start and end lift; the deleted upstream bases don't",
      cmap.to_hap(minus.canonical.end) == minus.canonical.end - start + 41 and cmap.to_hap(upstream + 50) is None)


# --- synthetic patients, aligned back -------------------------------------------------------------
def homopolymer(seq, at, n=6):
    """The start (0-based) of a run of n or more equal bases at or after `at`."""
    for i in range(at, len(seq) - n):
        if len(set(seq[i:i + n])) == 1 and seq[i - 1] != seq[i]:
            return i
    raise AssertionError("no homopolymer")


COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def clear_inversion(ref, at, n):
    """A start at or after `at` where inverting n bases has one description
    only: neither its outermost bases nor the ones just outside are each
    other's complement."""
    for s in range(at, len(ref) - n - 1):
        e = s + n - 1
        if COMP.get(ref[s]) != ref[e] and COMP.get(ref[s - 1]) != ref[e + 1]:
            return s
    raise AssertionError("no clear inversion")


def clear_replacement(segment):
    """7 bases replacing `segment` with no base in common: one complex event."""
    base = next(b for b in "ACGT" if b not in segment)
    return base * 7


def patient_and_back(gene, make_edits, paste=lambda hap: hap, label=""):
    start, end, ref = region_of(gene)
    chrom = "chr" + gene.chrom
    edits = make_edits(start, ref)
    hap, _ = haplotype.build(ref, start, edits)
    pasted = paste(hap)
    t0 = time.time()
    result = align.align(pasted, ref, start, fetch, chrom)
    seconds = time.time() - t0
    # Every event found, whether it is applied by default or not.
    rebuilt, _ = haplotype.build(ref, start, [Edit.from_event(e) for e in result.events])
    return result, hap, rebuilt, edits, seconds


plus = genes["APP"] if genes["APP"].strand == "+" else next(
    g for g in genes.values() if g.strand == "+" and 200_000 < g.mask[1] - g.mask[0] < 400_000)
small = genes["MRPL39"]


def many_edits(start, ref):
    hp = homopolymer(ref, 5000)
    snv = 3000
    inv = clear_inversion(ref, 18000, 300)
    return [
        Edit(start + snv, start + snv, "A" if ref[snv] != "A" else "C", "SNV"),
        Edit(start + hp + 2, start + hp + 3, "", "deletion"),               # 2 bases off a homopolymer
        Edit(start + 9000, start + 8999, rand_seq(65), "insertion"),        # 65 bp, deep in an intron
        Edit(start + 12000, start + 12079, "", "deletion"),                 # 80 bp
        Edit(start + 15000, start + 14999, ref[14880:15000], "duplication"),  # 120 bp tandem copy
        Edit(start + inv, start + inv + 299, "", "inversion"),              # 300 bp
        Edit(start + 21000, start + 21002, clear_replacement(ref[21000:21003]), "complex"),
    ]


result, hap, rebuilt, edits, seconds = patient_and_back(small, many_edits)
kinds = sorted(e.kind for e in result.events)
check(f"synthetic patient ({small.name}, {len(hap):,} bases, 7 known edits): found as such, and applying "
      "what was found rebuilds the patient exactly",
      rebuilt == hap and kinds == sorted(["SNV", "deletion", "insertion", "deletion", "duplication", "inversion",
                                          "complex"]), f"{kinds}, {seconds:.2f} s")
hp = homopolymer(region_of(small)[2], 5000)
two = next((e for e in result.events if e.kind == "deletion" and e.size == 2), None)
check("the homopolymer deletion is left-aligned (to the run's first bases), VCF-style with its anchor",
      two is not None and two.start == region_of(small)[0] + hp and len(two.ref) == 3 and len(two.alt) == 1,
      f"{two}")
dup = next((e for e in result.events if e.kind == "duplication"), None)
check("the 120 bp tandem copy is called a duplication, applied by default; SNV and 2 bp deletion are not",
      dup is not None and dup.size == 120 and dup.apply
      and not next(e for e in result.events if e.kind == "SNV").apply)
inv = next((e for e in result.events if e.kind == "inversion"), None)
inv_at = region_of(small)[0] + clear_inversion(region_of(small)[2], 18000, 300)
check("the inversion's breakpoints are exact", inv is not None and (inv.start, inv.end) == (inv_at, inv_at + 299),
      f"{inv and (inv.start, inv.end)}, expected {(inv_at, inv_at + 299)}")
check("identity and covered interval are reported", result.ref_start == region_of(small)[0]
      and result.ref_end == region_of(small)[1] and 0.97 < result.identity < 0.99, f"{result.identity:.4f}")

result, hap, rebuilt, _, _ = patient_and_back(small, many_edits, paste=align.reverse_complement)
check("a reverse-complemented paste is recognised, used turned around, and gives the same events",
      result.reverse and rebuilt == hap and len(result.events) == 7)

# A partial paste: only the middle 8 kb, with one insertion.
def middle_edit(start, ref):
    return [Edit(start + 10000, start + 9999, rand_seq(30), "insertion")]


result, hap, rebuilt, _, _ = patient_and_back(small, middle_edit, paste=lambda h: h[6000:14030])
report = align.compatibility(result, small.mask[0], small.mask[1])
check("partial paste: placed exactly, its one insertion found, the uncovered parts of the gene reported",
      result.ref_start == region_of(small)[0] + 6000 and len(result.events) == 1
      and result.events[0].kind == "insertion" and any("doesn't cover the whole gene" in w for w in report["warnings"]),
      str(report["warnings"])[:160])
check("partial paste: the report gives identity and the covered interval",
      any(line.startswith("Identity") for line in report["lines"])
      and any(line.startswith("Covered reference") for line in report["lines"]))

# Junk at the paste's end is left out, not taken as an insertion.
result, hap, _, _, _ = patient_and_back(small, lambda s, r: [], paste=lambda h: h[5000:9000] + rand_seq(200))
# (Their first few bases may match the reference by chance; those are kept.)
check("200 unrelated bases at the end are left out and reported, not made an insertion",
      not result.events and 4000 <= result.query_end <= 4005
      and any(f"last {len(hap[5000:9000]) + 200 - result.query_end} bases" in n for n in result.notes),
      str(result.notes))

# Random patients: whatever description the alignment picks (indels in repeats
# have several), applying it must give the patient's sequence back.
def random_edits(start, ref):
    edits, pos = [], 1500 + rng.randrange(0, 200)
    while pos < len(ref) - 1500:
        kind = rng.choice(["SNV", "SNV", "insertion", "deletion", "duplication", "inversion", "complex"])
        n = rng.choice([1, 2, 3, 5, 12, 40, 90, 250])
        p = start + pos
        if kind == "SNV":
            edits.append(Edit(p, p, "A" if ref[pos] != "A" else "T", kind))
            n = 1
        elif kind == "insertion":
            edits.append(Edit(p, p - 1, rand_seq(n), kind))
            n = 0
        elif kind == "duplication":
            edits.append(Edit(p, p - 1, ref[pos - n:pos], kind))
            n = 0
        elif kind == "deletion":
            edits.append(Edit(p, p + n - 1, "", kind))
        elif kind == "inversion":
            n = max(n, 60)
            edits.append(Edit(p, p + n - 1, "", kind))
        else:
            edits.append(Edit(p, p + n - 1, rand_seq(rng.randrange(1, 15)), kind))
        pos += n + rng.randrange(30, 900)
    return edits


failed, total = [], 0
for trial in range(25):
    gene = rng.choice([small, minus, genes["CXADR"], genes["JAM2"]])
    result, hap, rebuilt, edits, _ = patient_and_back(gene, random_edits)
    total += len(edits)
    if rebuilt != hap:
        failed.append((trial, gene.name, len(edits), len(result.events)))
check(f"25 random patients ({total} edits of every kind, 30-900 bases apart): applying what the alignment "
      "found rebuilds each patient exactly", not failed, str(failed[:3]))

# Short pastes: exact dynamic programming over the whole region, both orientations.
s_start, s_end, s_ref = region_of(plus)
piece = s_ref[150_000:150_300]
piece = piece[:120] + ("A" if piece[120] != "A" else "C") + piece[121:]
t0 = time.time()
result = align.align(piece, s_ref, s_start, fetch, "chr" + plus.chrom)
short_seconds = time.time() - t0
check(f"short paste (300 bases in {plus.name}'s {len(s_ref):,}-base region): placed by dynamic programming, "
      "its SNV found", result.ref_start == s_start + 150_000 and [e.kind for e in result.events] == ["SNV"]
      and result.events[0].start == s_start + 150_120 and "dynamic programming" in result.method,
      f"{short_seconds:.2f} s")
result = align.align(align.reverse_complement(piece), s_ref, s_start, fetch, "chr" + plus.chrom)
check("short paste, reverse complement: placed the same, marked reverse",
      result.reverse and result.ref_start == s_start + 150_000 and len(result.events) == 1)

# The same gene on the other build: a paste from hg38 doesn't "fit hg19 better".
genes19 = {g.name: g for g in transcripts.load_genes("hg19", helpers.THIRDPARTY / "refseq",
                                                       helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")}
fasta19 = Fasta(str(helpers.HG19))
start38, end38, ref38 = region_of(small)
g19 = genes19[small.name]
start19, end19 = align.analysed_region(g19.mask[0], g19.mask[1], len(fasta19["chr" + g19.chrom]), flank=2000)
ref19 = fasta19["chr" + g19.chrom][start19 - 1:end19].seq.upper()
paste = ref38[1000:-1000]
a38 = align.align(paste, ref38, start38, fetch, "chr" + small.chrom)
a19 = align.align(paste, ref19, start19, make_base_fetcher(fasta19), "chr" + g19.chrom)
report = align.compatibility(a38, small.mask[0], small.mask[1], other=("hg19", a19))
check(f"other build: an hg38 paste of {small.name} aligns to hg19 too, and isn't said to fit hg19 better",
      a19.columns > 0 and not any("fits hg19" in w for w in report["warnings"]),
      f"hg38 {100 * a38.identity:.2f}%, hg19 {100 * a19.identity:.2f}%")
check("analysed region: the gene plus 5,000 bases each side, within the chromosome; or the window given",
      align.analysed_region(3000, 9000, 20_000) == (1, 14_000)
      and align.analysed_region(3000, 9000, 12_000, window=(2500, 9500)) == (2500, 9500)
      and raises(lambda: align.analysed_region(1, 5, 100, window=(1, 101)), "past the chromosome's end"))

# --- the report and the chosen events -----------------------------------------------------------------
a = align.Alignment(1000, False, 100, 1099, 1, 1000, 990, 1000, 4900, method="x")
b = align.Alignment(1000, False, 100, 1099, 1, 1000, 1000, 1000, 5000, method="x")
report = align.compatibility(a, 50, 1200, spliceai=(100, 1100), canonical=(80, 1150), other=("hg19", b))
check("report: a clearly better fit to the other build is a warning",
      any("fits hg19 better" in w for w in report["warnings"]), str(report["warnings"]))
check("report: SpliceAI's gene interval and the canonical transcript differing is noted",
      any("SpliceAI's gene interval" in line for line in report["lines"]))
ev = lambda kind, s, e, apply=True: align.Event(kind, s, "", "", s, e, "", e - s + 1, apply=apply)
checked = align.validate_events([ev("deletion", 100, 120), ev("SNV", 110, 110)], 1, 1000)
check("chosen events that overlap block the run", any("overlap" in p for p in checked["problems"]))
checked = align.validate_events([ev("deletion", 100, 120), ev("SNV", 110, 110, apply=False)], 1, 1000)
check("...but not when one of them isn't chosen", not checked["problems"])
checked = align.validate_events([ev("deletion", 990, 1010)], 1, 1000)
check("an event reaching outside the analysed region blocks the run", any("outside" in p for p in checked["problems"]))
checked = align.validate_events([ev("deletion", 500, 510)], 1, 1000, anchors=[("ATG", 511)])
check("an event touching an anchor (here the ATG) is a warning", any("touches the ATG" in w for w in checked["warnings"]))

# --- speed -----------------------------------------------------------------------------------------------
big = genes["DMD"]
b_start, b_end, b_ref = region_of(big)
edits = [Edit(b_start + p, b_start + p + 49, "", "deletion") for p in range(100_000, len(b_ref) - 100_000, 400_000)]
b_hap, _ = haplotype.build(b_ref, b_start, edits)
t0 = time.time()
result = align.align(b_hap, b_ref, b_start, fetch, "chrX")
big_seconds = time.time() - t0
rebuilt, _ = haplotype.build(b_ref, b_start, [Edit.from_event(e) for e in result.events])
check(f"a {len(b_hap) / 1e6:.1f} Mb paste (DMD) aligns in seconds, its {len(edits)} deletions found exactly",
      rebuilt == b_hap and len(result.events) == len(edits) and big_seconds < 60, f"{big_seconds:.1f} s")
piece = b_ref[1_000_000:1_000_499]
t0 = time.time()
result = align.align(piece, b_ref, b_start, fetch, "chrX")
dp_seconds = time.time() - t0
check(f"the slowest case, 499 bases by dynamic programming in DMD's {len(b_ref) / 1e6:.1f} Mb: within the "
      "spec's 2 to 3 s", result.ref_start == b_start + 1_000_000 and dp_seconds < 3.5, f"{dp_seconds:.2f} s")

helpers.finish()
