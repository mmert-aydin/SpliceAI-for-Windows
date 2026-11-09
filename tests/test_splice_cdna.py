"""Phase 6: transcripts from SpliceAI's probabilities (splice_paths.py), and
their cDNA, protein and comparisons (cdna.py). Synthetic probability arrays
(canonical, forced skip, forced retention, pseudoexon, no skipping bias, the
k best against brute force), then real reference transcripts on both
strands.

    .venv\\Scripts\\python tests\\test_splice_cdna.py
"""
import math
import random

import helpers
from helpers import check

helpers.watchdog(600)
helpers.setup()

import numpy as np  # noqa: E402
from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline import cdna, haplotype, splice_paths as sp, transcripts  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402

rng = random.Random(5)
OFFSET = 1001


def arrays(sites, n=1500):
    """acceptor, donor arrays for positions OFFSET.. with sites {("A"|"D", pos): p}."""
    acc, don = np.zeros(n), np.zeros(n)
    for (kind, pos), p in sites.items():
        (acc if kind == "A" else don)[pos - OFFSET] = p
    return acc, don


EXONS = [(1001, 1100), (1301, 1400), (1601, 1700), (2101, 2200)]
CANON = {("D", 1100): .95, ("A", 1301): .95, ("D", 1400): .95, ("A", 1601): .95, ("D", 1700): .95, ("A", 2101): .95}


def best(sites, strand="+", **kw):
    acc, don = arrays(sites)
    return sp.top_paths(acc, don, OFFSET, strand, 1001, 2200, **kw)


paths = best(CANON)
check("synthetic: the canonical transcript is the top path", paths[0].exons == EXONS, str(paths[0].exons))
check("synthetic: relative scores are a softmax over the paths reported (sum 1, best first)",
      abs(sum(p.relative for p in paths) - 1) < 1e-12 and paths[0].relative == max(p.relative for p in paths)
      and len(paths) <= 5)
skip = {**CANON, ("A", 1601): 0.0}
check("synthetic: no acceptor for exon 3 forces it to be skipped",
      best(skip)[0].exons == [EXONS[0], EXONS[1], EXONS[3]], str(best(skip)[0].exons))
retain = {**CANON, ("D", 1100): 0.0, ("A", 1301): 0.0}
check("synthetic: no donor and acceptor around intron 1 keeps it (intron retention, no special rule)",
      best(retain)[0].exons == [(1001, 1400), EXONS[2], EXONS[3]], str(best(retain)[0].exons))
pseudo = {**CANON, ("A", 1800): .9, ("D", 1849): .9}
check("synthetic: a strong acceptor and donor deep in an intron become a pseudoexon",
      best(pseudo)[0].exons == [EXONS[0], EXONS[1], EXONS[2], (1800, 1849), EXONS[3]], str(best(pseudo)[0].exons))
moderate = {**CANON, ("A", 1601): .6, ("D", 1700): .6}
check("synthetic: no skipping bias -- an exon with two moderate sites (0.6) is kept, though scoring used sites "
      "alone would skip it (log 0.6 + log 0.6 < 0)", best(moderate)[0].exons == EXONS
      and 2 * math.log(.6) < 0, str(best(moderate)[0].exons))
check("synthetic: a weak site (0.1) isn't used just because it's a candidate",
      best({**CANON, ("A", 1850): .1, ("D", 1900): .1})[0].exons == EXONS)
short = {**CANON, ("A", 1800): .99, ("D", 1810): .99}
check("synthetic: an exon shorter than the minimum isn't allowed; at the default of 3 this one is",
      all((1800, 1810) not in p.exons for p in best(short, min_exon=25))
      and any((1800, 1810) in p.exons for p in best(short)))
check("synthetic: minimum exon and intron lengths are configurable",
      any((1800, 1810) in p.exons for p in best(short, min_exon=5)))

# The same on the minus strand: mirror everything.
mirror = lambda x: 3201 - x
macc, mdon = arrays(CANON)
m_sites = {("A" if k == "A" else "D", mirror(p)): v for (k, p), v in CANON.items()}
m_acc, m_don = arrays(m_sites)
mp = sp.top_paths(m_acc, m_don, OFFSET, "-", mirror(1001), mirror(2200))
check("synthetic, minus strand: the mirrored transcript, in transcript order",
      mp[0].exons == [tuple(sorted((mirror(a), mirror(b)))) for a, b in EXONS] and mp[0].exons[0][1] == mirror(1001),
      str(mp[0].exons))


# Exact k-best: against every valid path, on random small cases.
def brute_force(acc, don, length, min_exon=sp.MIN_EXON, min_intron=sp.MIN_INTRON, threshold=sp.SITE_THRESHOLD):
    """Every path's score, by enumeration (transcript coordinates 0..length, plus strand)."""
    donors = [t for t in range(length) if don[t] >= threshold]
    acceptors = [t for t in range(1, length + 1) if acc[t] >= threshold]
    out = []

    def walk(exon_start, used):
        out.append(sp.path_score(exons_of(used, length), acc, don, OFFSET_T, "+"))
        for d in donors:
            if d < exon_start or (used and d - exon_start + 1 < min_exon):
                continue
            for a in acceptors:
                if a - d - 1 >= min_intron:
                    walk(a, used + [d, a])

    walk(0, [])
    return sorted(out, reverse=True)


OFFSET_T = 0


def exons_of(used, length):
    bounds = [0] + used + [length]
    return list(zip(bounds[0::2], bounds[1::2]))


mismatch = []
for case in range(60):
    length = 420
    acc, don = np.zeros(length + 1), np.zeros(length + 1)
    for _ in range(rng.randrange(4, 9)):
        acc[rng.randrange(1, length + 1)] = rng.choice([rng.uniform(0.05, 1), 1.0, 0.0])
        don[rng.randrange(0, length)] = rng.choice([rng.uniform(0.05, 1), 1.0, 0.0])
    scores = brute_force(acc, don, length)[:5]
    got = [p.score for p in sp.top_paths(acc, don, OFFSET_T, "+", 0, length)]
    if len(got) != len(scores) or any(abs(a - b) > 1e-9 for a, b in zip(got, scores)):
        mismatch.append((case, got[:3], scores[:3]))
check("the k-best dynamic programme gives the brute force's 5 best scores (60 random cases)", not mismatch,
      str(mismatch[:2]))
check("path_score() gives the score top_paths() reports for its path",
      abs(sp.path_score(paths[0].exons, *arrays(CANON), OFFSET, "+") - paths[0].score) < 1e-9)

# --- real reference transcripts ---------------------------------------------------------------
fasta = Fasta(str(helpers.HG38))
genes = {g.name: g for g in transcripts.load_genes("hg38", helpers.THIRDPARTY / "refseq",
                                                     helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")}
known = {"BRCA2": ("+", 3418), "APP": ("-", 770)}
refs = {}
for name, (strand, aa) in known.items():
    t = genes[name].canonical
    record = fasta["chr" + t.chrom.replace("chr", "")]
    fetch = lambda a, b, record=record: record[a - 1:b].seq
    ref = cdna.reference(fetch, t)
    refs[name] = (t, fetch, ref)
    stop_base = t.cds[1] if t.strand == "+" else t.cds[0]
    regions = ref.regions()
    check(f"{name} ({strand}, {t.accession}): its cDNA translates from the ATG to the annotated stop, {aa} aa",
          t.strand == strand and ref.protein.startswith("M") and len(ref.protein) == aa and "*" not in ref.protein
          and ref.sequence[ref.cds_start:ref.cds_start + 3] == "ATG" and ref.index(stop_base) == ref.cds_end - 1,
          f"{len(ref.protein)} aa")
    check(f"{name}: 5'UTR, CDS and 3'UTR cover the cDNA exactly, CDS a whole number of codons",
          [r[0] for r in regions] == ["5'UTR", "CDS", "3'UTR"] and regions[0][1] == 0
          and regions[-1][2] == len(ref.sequence) and (regions[1][2] - regions[1][1]) % 3 == 0
          and sum(b - a + 1 for a, b in t.exons) == len(ref.sequence))
    check(f"{name}: summary of the reference itself: no change, no PTC", (s := cdna.summarize(ref, ref)).frame ==
          "no change" and not s.ptc and cdna.hgvs_r(ref, ref) == "r.(=)" and cdna.hgvs_p(ref, ref) == "p.(=)")

t, fetch, ref = refs["APP"]
exons = list(reversed(t.exons))                 # transcript order (minus strand)
coding = [k for k in range(1, len(exons) - 1)
          if ref.cds_start < ref.exon_starts[k] and ref.exon_starts[k + 1] < ref.cds_end - 3]
in_frame = next(k for k in coding if (exons[k][1] - exons[k][0] + 1) % 3 == 0)
shifting = next(k for k in coding if (exons[k][1] - exons[k][0] + 1) % 3 != 0)
for k, expect in ((in_frame, "in frame"), (shifting, "frameshift")):
    skipped = cdna.build(fetch, "-", exons[:k] + exons[k + 1:], t.cds[1])
    s = cdna.summarize(skipped, ref)
    first = cdna.c_position(ref, ref.exon_starts[k])
    last = cdna.c_position(ref, ref.exon_starts[k + 1] - 1)
    r_desc = cdna.hgvs_r(ref, skipped)
    check(f"APP skipping exon {k + 1} ({exons[k][1] - exons[k][0] + 1} nt): {expect}, r.({first}_{last}del) "
          "(3' rule aside, same length)", s.frame == expect and r_desc.endswith("del)")
          and len(skipped.sequence) == len(ref.sequence) - (exons[k][1] - exons[k][0] + 1), f"{s.frame}, {r_desc}")
    p_desc = cdna.hgvs_p(ref, skipped)
    check(f"APP skipping exon {k + 1}: p. is {'an in-frame deletion' if expect == 'in frame' else 'a frameshift'}",
          ("del" in p_desc and "fs" not in p_desc) if expect == "in frame" else ("fs" in p_desc or "Ter)" in p_desc),
          p_desc)
    check(f"APP skipping exon {k + 1}: " + ("the normal stop is kept, no premature stop" if expect == "in frame"
                                            else "a premature stop codon, NMD assessed"),
          (not s.ptc and not s.nmd) if expect == "in frame" else (s.ptc and bool(s.nmd)), f"{s.ptc} {s.nmd}")
atg_exon = next(k for k in range(len(exons)) if ref.exon_starts[k] <= ref.cds_start
                < (ref.exon_starts[k + 1] if k + 1 < len(exons) else len(ref.sequence)))
lost = cdna.build(fetch, "-", exons[:atg_exon] + exons[atg_exon + 1:], t.cds[1])
check("APP without the exon of its ATG: ATG lost, p.(Met1?), no start invented",
      lost.cds_start is None and cdna.summarize(lost, ref).frame == "ATG lost" and cdna.hgvs_p(ref, lost) == "p.(Met1?)"
      and not lost.protein)
cut = cdna.build(fetch, "-", exons[:-1], t.cds[1])
check("APP ending before its stop codon: 'no stop codon' noted", cut.cds_end is None and any(
    "no stop codon" in n for n in cut.notes))

# "Change vs WT" reads in the reference's c. numbering: the WT prediction lacks
# exon a, the patient lacks a and b -> the same r. as b lost from the reference.
a, b = in_frame, next(k for k in coding if abs(k - in_frame) > 1)
without = lambda *ks: cdna.build(fetch, "-", [e for k, e in enumerate(exons) if k not in ks], t.cds[1])
wt_like = without(a)
via_wt, direct = cdna.hgvs_r(wt_like, without(a, b), cdna.number_on(ref, wt_like)), cdna.hgvs_r(ref, without(b))
check(f"APP vs a WT prediction without exon {a + 1}: exon {b + 1} lost too is in the reference's numbering ({direct})",
      via_wt == direct, f"{via_wt} / {direct}")

# NMD, on a synthetic transcript: 4 exons, stop codon placed at chosen distances.
def synthetic(stop_at, exon_lengths=(60, 120, 120, 200)):
    """A plus-strand transcript ATG at cDNA 10, stop codon at cDNA index stop_at."""
    seq = list(("ACGTCAGCCC" + "GCC" * 400)[:sum(exon_lengths)])
    seq[10:13] = "ATG"
    seq[stop_at:stop_at + 3] = "TAA"
    genome = "".join(seq)
    exons, pos = [], 1
    for n in exon_lengths:
        exons.append((pos, pos + n - 1))
        pos += n + 1000
    flat = {}
    for (a, b), off in zip(exons, [sum(exon_lengths[:i]) for i in range(len(exons))]):
        for x in range(a, b + 1):
            flat[x] = genome[off + x - a]
    return cdna.build(lambda a, b: "".join(flat.get(x, "N") for x in range(a, b + 1)), "+", exons, 11)


far = synthetic(178)        # last junction at 300: stop ends 119 nt before it
check("NMD: a stop codon 119 nt before the last exon-exon junction -> expected",
      cdna.nmd(far).startswith("expected") and "119 nt" in cdna.nmd(far), cdna.nmd(far))
near = synthetic(250)       # ends 47 nt before the junction
check("NMD: a stop codon 47 nt before the last junction escapes (50-55 nt rule)",
      cdna.nmd(near).startswith("escapes") and "47 nt" in cdna.nmd(near), cdna.nmd(near))
last = synthetic(319)        # in frame (ATG at 10), inside the last exon (from 300)
check("NMD: a stop codon in the last exon escapes", cdna.nmd(last) == "escapes: stop codon in the last exon", cdna.nmd(last))
check("NMD: without a stop codon it is not assessed", cdna.nmd(synthetic(320)) == "not assessed: no stop codon")
early = synthetic(40, exon_lengths=(60, 300, 120, 200))
check("NMD: a start-proximal stop is flagged as possibly escaping", "start-proximal" in cdna.nmd(early), cdna.nmd(early))

# --- exon concordance ------------------------------------------------------------------------------
ref_ex = [(100, 200), (500, 600), (900, 1000), (1300, 1400)]
other = [(100, 200), (500, 637), (1150, 1200), (1300, 1400)]
rows, summary = cdna.concordance(ref_ex, other, "+")
status = [r["status"] for r in rows]
check("concordance: identical, end shifted +37 (a cryptic donor), missing, extra (novel)",
      status == ["identical", "end shifted +37", "missing", "identical", "extra (novel)"]
      and summary == "2 identical, 1 shifted, 1 missing, 1 extra (novel)", f"{status} | {summary}")
rows, _ = cdna.concordance(ref_ex, [(100, 600), (900, 1000), (1300, 1400)], "+")
check("concordance: a retained intron shows as both exons in one, with their shifts",
      rows[0]["status"].startswith("end shifted +400") and "reference exon 2" in rows[0]["status"]
      and rows[1]["status"].startswith("start shifted -400"), str([r["status"] for r in rows[:2]]))
m_ref = [(1300, 1400), (900, 1000)]
rows, _ = cdna.concordance(m_ref, [(1300, 1400), (880, 1000)], "-")
check("concordance, minus strand: + is still further 3' (an exon ending 20 nt further 3')",
      rows[1]["status"] == "end shifted +20", rows[1]["status"])
ref_seq = "".join(rng.choice("ACGT") for _ in range(2000))
hap, cmap = haplotype.build(ref_seq, 1, [Edit(1101, 1100, "".join(rng.choice("ACGT") for _ in range(300)), "insertion")])
mapped = cdna.to_reference_exons([(100, 200), (1150, 1250), (1500, 1600)], cmap)
check("patient exons onto the reference: one inside inserted sequence becomes extra (None), others map back",
      mapped == [(100, 200), None, (1200, 1300)], str(mapped))

# --- HGVS details and the pairwise view ------------------------------------------------------------
t, fetch, ref = refs["BRCA2"]
seq = ref.sequence


def altered(new_sequence):
    """ref's transcript with another sequence, translated from the same ATG."""
    return cdna.translate(cdna.CDNA(ref.exons, "+", new_sequence, ref.exon_starts), ref.cds_start)


i = ref.cds_start + 99                      # c.100, the first base of codon 34
new_base = next(b for b in "ACGT" if b != seq[i] and seq[i:i + 3].replace(seq[i], b, 1) not in cdna.STOPS)
mut = altered(seq[:i] + new_base + seq[i + 1:])
r_desc, p_desc = cdna.hgvs_r(ref, mut), cdna.hgvs_p(ref, mut)
rna = lambda b: b.lower().replace("t", "u")
check("HGVS: one base changed at c.100 is r.(100x>y), and p. names codon 34 (or p.(=) if silent)",
      r_desc == f"r.(100{rna(seq[i])}>{rna(new_base)})"
      and (p_desc == "p.(=)" or "34" in p_desc), f"{r_desc} {p_desc}")
ins = altered(seq[:ref.cds_start + 30] + "GGG" + seq[ref.cds_start + 30:])
check("HGVS: a 3-base insertion in the CDS: r. ins of 3 bases, in frame, p. an insertion or delins",
      "ins" in cdna.hgvs_r(ref, ins) and cdna.summarize(ins, ref).frame == "in frame"
      and "fs" not in cdna.hgvs_p(ref, ins) and "ins" in cdna.hgvs_p(ref, ins),
      f"{cdna.hgvs_r(ref, ins)} {cdna.hgvs_p(ref, ins)}")
a, m, b = cdna.pairwise("ACGTACGTTTGCA", "ACGTTTGCA")
check("pairwise: aligned strings of equal length, the deletion as gaps", len(a) == len(m) == len(b)
      and b.count("-") == 4 and a.replace("-", "") == "ACGTACGTTTGCA", f"{a} / {b}")
pa, _, pb = cdna.pairwise("MKTAYIAKQR", "MKTAYQR", protein=True)
check("pairwise: proteins too", pb.count("-") == 3 and pa == "MKTAYIAKQR", f"{pa} / {pb}")

# --- the wild type from the database's values ----------------------------------------------------
import os  # noqa: E402

from spliceai_pipeline import refdb  # noqa: E402

real_db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
if os.path.exists(os.path.join(real_db, "manifest.json")):
    reader = refdb.Reader(real_db)
    for name in ("BRCA2", "SAMSN1", "MRPL39"):
        g = genes[name]
        c, plan, v = g.canonical, reader.plan(name), reader.values(name)
        tss, tes = (c.start, c.end) if c.strand == "+" else (c.end, c.start)
        top = sp.top_paths(v[:, 0], v[:, 1], plan["start"], c.strand, tss, tes)[0]
        ex = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
        _, summary = cdna.concordance(ex, top.exons, c.strand)
        check(f"wild type, {name} ({c.strand}): the top path from the database's values is the canonical transcript",
              top.exons == ex, summary)
else:
    helpers.skip("wild type from the database's values", "the Phase 2 database isn't on this PC")

helpers.finish()
