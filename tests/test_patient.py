"""Phase 5: the patient's prediction on top of the reference score database
(haplotype.predict_patient, spec 7.2). Mandatory: predicting only near the
edits gives what predicting the whole haplotype gives, and the deltas at a
variant's +/-50 window are SpliceAI's own get_delta_scores' raw DS/DP (spec
5.6 rule) where the conventions agree.

    .venv\\Scripts\\python tests\\test_patient.py
"""
import os
import random
import tempfile

import helpers
from helpers import check

helpers.watchdog(1200)
helpers.setup()

import numpy as np  # noqa: E402

from spliceai_pipeline import align, haplotype, predict, refdb, transcripts  # noqa: E402
from spliceai_pipeline import refdb_validate as rv  # noqa: E402
from spliceai_pipeline.control import RunCancelled  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402

rng = random.Random(11)
rand_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
REFSEQ, MANE = helpers.THIRDPARTY / "refseq", helpers.MANE_DIR

# The reference values come from the database this PC built (Phase 2) if it is
# there, else they are computed live (the fallback) into a throw-away cache.
real_db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
reader = refdb.Reader(real_db) if os.path.exists(os.path.join(real_db, "manifest.json")) else None
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=reader, cache_dir=tempfile.mkdtemp(),
                               refseq_dir=str(REFSEQ), mane_dir=str(MANE))
if reader is None:
    helpers.skip("reference values from the Phase 2 database", "not built on this PC; computed live instead")
predictor = scores.predictor()
fasta = predictor.annotator.ref_fasta
genes = {g.name: g for g in transcripts.load_genes("hg38", REFSEQ, MANE / "MANE.GRCh38.v1.5.summary.txt.gz")}


def setup_gene(name):
    gene = genes[name]
    record = fasta["chr" + gene.chrom]
    start, end = align.analysed_region(gene.mask[0], gene.mask[1], len(record))
    ref = record[start - 1:end].seq.upper()
    wt, plan, source = scores.gene(name)
    return gene, record, start, ref, wt, plan, source


def fresh_reference(record, plan):
    return predictor.predict(lambda a, b: record[a - 1:b].seq, len(record), plan["strand"], *plan["mask"],
                             plan["start"], plan["end"])


def full_patient(hap, cmap, plan):
    start, end, mask = haplotype.patient_region(cmap, plan)
    return predictor.predict(lambda a, b: hap[a - 1:b], len(hap), plan["strand"], mask[0], mask[1], start, end)


def edits_in(gene, start, ref):
    """Seven kinds of edit spread over the gene's mask, one close to its start."""
    m0 = gene.mask[0] - start
    span = gene.mask[1] - gene.mask[0]
    at = lambda f: m0 + int(span * f)
    p = [at(0.002), at(0.15), at(0.3), at(0.45), at(0.6), at(0.75), at(0.9)]
    return [
        Edit(start + p[0], start + p[0], "A" if ref[p[0]] != "A" else "C", "SNV"),
        Edit(start + p[1], start + p[1] - 1, rand_seq(65), "insertion"),
        Edit(start + p[2], start + p[2] + 79, "", "deletion"),
        Edit(start + p[3], start + p[3] - 1, ref[p[3] - 120:p[3]], "duplication"),
        Edit(start + p[4], start + p[4] + 299, "", "inversion"),
        Edit(start + p[5], start + p[5] + 2, "GGGGGGG", "complex"),
        Edit(start + p[6], start + p[6] + 1, "", "deletion"),
    ]


# --- local recomputation = full recomputation (mandatory) ---------------------------------
# Large enough that most of each gene is reused from the database, not predicted again.
big = lambda strand: next(g for g in sorted(genes.values(), key=lambda g: (g.chrom != "21", g.start))
                          if g.strand == strand and 150_000 <= g.mask[1] - g.mask[0] <= 300_000)
plus, minus = big("+"), big("-")
check("the genes: one on each strand, 150-300 kb", plus.strand == "+" and minus.strand == "-",
      f"{plus.name} {plus.mask[1] - plus.mask[0]:,}, {minus.name} {minus.mask[1] - minus.mask[0]:,}")
for gene in (plus, minus):
    gene, record, start, ref, wt, plan, source = setup_gene(gene.name)
    exact = fresh_reference(record, plan)
    check(f"{gene.name} ({gene.strand}): reference values from the {source} are the fresh prediction (float16 at most)",
          float(np.abs(wt - exact).max()) <= 2 ** -12 + 1e-6)
    hap, cmap = haplotype.build(ref, start, edits_in(gene, start, ref))
    full = full_patient(hap, cmap, plan)
    seen = []
    local = haplotype.predict_patient(predictor, hap, cmap, plan, exact, on_progress=lambda d, t: seen.append((d, t)))
    recomputed = sum(b - a + 1 for a, b in local.windows)
    diff = float(np.abs(local.values - full).max())
    check(f"{gene.name}: predicting only within 5,000 bases of the 7 edits gives the whole haplotype's prediction "
          "(exact reference values), and most of the gene is reused", diff <= 1e-6 and recomputed < 0.6 * len(full),
          f"max |diff| {diff:.2g}; {recomputed:,} of {len(full):,} bases predicted again, {local.seconds:.1f} s")
    local16 = haplotype.predict_patient(predictor, hap, cmap, plan, wt)
    diff16 = float(np.abs(local16.values - full).max())
    check(f"{gene.name}: ...and with the database's values (float16) within its rounding",
          diff16 <= 2 ** -12 + 1e-6, f"max |diff| {diff16:.2g}")
    check(f"{gene.name}: progress is reported after each window, ending at the total",
          seen and seen[-1][0] == seen[-1][1] == recomputed and len(seen) == len(local.windows))
    check(f"{gene.name}: the estimate before the run is in the right range",
          0.2 < haplotype.estimate_seconds(local.windows) / max(local.seconds, 0.1) < 5,
          f"estimated {haplotype.estimate_seconds(local.windows):.1f} s, took {local.seconds:.1f} s")

gene, record, start, ref, wt, plan, source = setup_gene(plus.name)
hap, cmap = haplotype.build(ref, start, [])
same = haplotype.predict_patient(predictor, hap, cmap, plan, wt)
check("no edits: nothing is predicted again, and the patient is the reference", not same.windows
      and np.array_equal(same.values, wt))


def stop():
    raise RunCancelled("stop")


hap, cmap = haplotype.build(ref, start, edits_in(gene, start, ref))
try:
    haplotype.predict_patient(predictor, hap, cmap, plan, wt, checkpoint=stop)
    check("a cancelled prediction stops", False)
except RunCancelled:
    check("a cancelled prediction stops (RunCancelled from the checkpoint)", True)
everything = [Edit(gene.mask[0] - 100, gene.mask[1] + 100, "", "deletion")]
hap, cmap = haplotype.build(ref, start, everything)
try:
    haplotype.predict_patient(predictor, hap, cmap, plan, wt)
    check("deleting the whole gene is refused", False)
except ValueError as exc:
    check("deleting the whole gene is refused, saying so", "delete the whole gene" in str(exc))

# --- the same deltas as the VCF pipeline (mandatory) ---------------------------------------------
annotator = predictor.annotator
problems, checked, worst, strong = [], 0, 0.0, 0
for gene in (plus, minus):
    gene, record, start, ref, wt, plan, source = setup_gene(gene.name)
    lo, hi = plan["vcf_comparable"]
    chrom = "chr" + gene.chrom
    # Random positions, then the canonical transcript's splice sites: the base
    # just outside each internal exon end (donor or acceptor, by strand), where
    # DS is large and DP must agree too.
    sites = [p for a, b in gene.canonical.exons[1:-1][:4] for p in (a - 1, b + 1)]
    for k in range(12 + len(sites)):
        pos = rng.randrange(lo + 60, hi - 60) if k < 12 else sites[k - 12]
        if not lo + 60 <= pos <= hi - 60:
            continue
        base = ref[pos - start:pos - start + 5]
        if "N" in base:
            continue
        kind = ["SNV", "insertion", "deletion"][k % 3]
        if kind == "SNV":
            vref, valt = base[0], "A" if base[0] != "A" else "G"
            edit = Edit(pos, pos, valt, "SNV")
        elif kind == "insertion":
            vref, valt = base[0], base[0] + rand_seq(1)
            edit = Edit(pos + 1, pos, valt[1:], "insertion")
        else:
            n = 1 + k % 4
            vref, valt = base[:n + 1], base[0]
            edit = Edit(pos + 1, pos + n, "", "deletion")
        _, live = rv.delta_scores(annotator, chrom, pos, vref, valt, masked=False)
        if gene.name not in live:
            continue
        hap, cmap = haplotype.build(ref, start, [edit])
        pred = haplotype.predict_patient(predictor, hap, cmap, plan, wt)
        y_ref = np.zeros((1, rv.COV, 3), np.float32)
        y_ref[0, :, 1:] = wt[pos - rv.HALF - plan["start"]:pos + rv.HALF + 1 - plan["start"]]
        y_alt = np.zeros((1, rv.COV, 3), np.float32)
        y_alt[0, :, 1:] = haplotype.ref_aligned(pred, cmap, pos - rv.HALF, pos + rv.HALF)
        ds, dp, _ = rv.scores_from(y_ref, y_alt, live[gene.name]["dist_exon"], masked=False)
        found, _flips = rv.compare(live[gene.name]["ds"], live[gene.name]["dp"], live[gene.name]["series"], ds, dp)
        worst = max(worst, max(abs(a - b) for a, b in zip(live[gene.name]["ds"], ds)))
        checked += 1
        strong += max(live[gene.name]["ds"]) >= 0.5
        if found:
            problems.append(f"{gene.name} {chrom}:{pos} {vref}>{valt}: {found}")
check(f"the engine's raw DS/DP at the variant's +/-50 window are get_delta_scores' (SNVs, 1-base insertions, "
      f"1-4 base deletions; random and at splice sites; both strands; spec 5.6 rule)",
      checked >= 25 and strong >= 6 and not problems,
      f"{checked} variants ({strong} with DS >= 0.5), max |dDS| {worst:.2g}; {problems[:2]}")

# --- end to end: a donor-site SNV, through the patient's transcripts ------------------------------
from spliceai_pipeline import cdna, splice_paths  # noqa: E402

for name, k, expect in (("SAMSN1", 3, ("in frame", False)), ("BRCA2", 5, ("frameshift", True))):
    gene, record, start, ref, wt, plan, source = setup_gene(name)
    c = gene.canonical
    tss, tes = (c.start, c.end) if c.strand == "+" else (c.end, c.start)
    ex = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
    a, b = ex[k - 1]
    site = b + 1 if c.strand == "+" else a - 1      # the intron's first base after exon k (donor +1)
    base = ref[site - start]
    hap, cmap = haplotype.build(ref, start, [Edit(site, site, "A" if base != "A" else "C", "SNV")])
    pred = haplotype.predict_patient(predictor, hap, cmap, plan, wt)
    patient = splice_paths.top_paths(pred.values[:, 0], pred.values[:, 1], pred.start, c.strand,
                                     cmap.to_hap(tss), cmap.to_hap(tes))
    wild = splice_paths.top_paths(wt[:, 0], wt[:, 1], plan["start"], c.strand, tss, tes)
    rows, summary = cdna.concordance(wild[0].exons, cdna.to_reference_exons(patient[0].exons, cmap), c.strand)
    ref_cdna = cdna.reference(lambda x, y: record[x - 1:y].seq, c)
    pc = cdna.build(lambda x, y: hap[x - 1:y], c.strand, patient[0].exons,
                    cmap.to_hap(c.cds[0] if c.strand == "+" else c.cds[1]))
    s = cdna.summarize(pc, ref_cdna)
    check(f"end to end, {name} ({c.strand}) donor +1 of exon {k}: the wild type is the canonical transcript, the "
          f"patient's top path skips exon {k} ({expect[0]}{', premature stop' if expect[1] else ''})",
          wild[0].exons == ex and [r["status"] for r in rows].count("missing") == 1
          and rows[k - 1]["status"] == "missing" and (s.frame, s.ptc) == expect,
          f"{summary}; {s.frame}, PTC {s.ptc}, {s.nmd[:30]}; {cdna.hgvs_r(ref_cdna, pc)} {cdna.hgvs_p(ref_cdna, pc)}")

helpers.finish()
