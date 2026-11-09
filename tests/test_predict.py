"""predict.py against SpliceAI's own get_delta_scores (spec 5.1, Phase 1).

The reference prediction get_delta_scores makes for a variant is taken from
inside it (the five models' outputs are recorded while it runs), not re-derived
here, and compared with predict.py for the same 101 positions: both strands,
windows crossing a gene's start or end (N-masking), a position in two
overlapping genes, and the ANKRD26 regression variant. Then the database's own
(extended) mask, where it gives the VCF tab's input, and chunk-size
independence.

    .venv\\Scripts\\python tests\\test_predict.py
"""
import helpers
from helpers import check

helpers.watchdog(900)
helpers.setup()

import numpy as np  # noqa: E402

from spliceai_pipeline import predict, refdb, transcripts  # noqa: E402
from spliceai_pipeline.score import _Record, get_delta_scores  # noqa: E402

# Float rounding only: the tf.function path vs Keras predict() measured 7e-9.
TOL = 1e-6

# --- masked_bases: where N goes ------------------------------------------------------
seq = "ACGTACGTAC"  # positions 1..10


def fetch(a, b):
    return seq[a - 1:b]


check("masked_bases: inside mask and sequence -> the bases", predict.masked_bases(fetch, 10, 1, 10, 3, 6) == "GTAC")
check("masked_bases: N before position 1 and after the sequence end",
      predict.masked_bases(fetch, 10, 1, 10, -1, 12) == "NN" + seq + "NN")
check("masked_bases: N outside the mask, bases inside it",
      predict.masked_bases(fetch, 10, 4, 6, 2, 8) == "NN" + "TAC" + "NN")
check("masked_bases: a window entirely outside the mask is all N",
      predict.masked_bases(fetch, 10, 8, 10, 1, 5) == "NNNNN")
check("masked_bases: length always b - a + 1",
      all(len(predict.masked_bases(fetch, 10, m1, m2, a, a + w)) == w + 1
          for m1 in (1, 4) for m2 in (6, 10) for a in range(-3, 12) for w in (0, 3, 20)))

check("the installed models are the pinned SpliceAI 1.3.1 files", predict.verify_models() == predict.MODEL_SHA256)


class Recorder:
    """Stands in for one model inside get_delta_scores and keeps its outputs."""

    def __init__(self, model):
        self.model, self.outputs = model, []

    def predict(self, x, *args, **kwargs):
        y = self.model.predict(x, *args, **kwargs)
        self.outputs.append(y)
        return y


def reference_inside_get_delta_scores(pred, chrom, pos):
    """{gene: (strand, tx_start, tx_end, y_ref)} for an SNV at chrom:pos, y_ref
    being the (101, 2) acceptor/donor reference prediction get_delta_scores
    made for that gene, in genomic order."""
    ann = pred.annotator
    ref = ann.ref_fasta[chrom][pos - 1:pos].seq.upper()
    alt = next(b for b in "ACGT" if b != ref)
    models, ann.models = ann.models, [Recorder(m) for m in ann.models]
    try:
        scores = get_delta_scores(_Record(chrom, pos, ref, alt), ann, 50, 0)
        recorders = ann.models
    finally:
        ann.models = models
    genes, strands, idxs = ann.get_name_and_strand(chrom, pos)
    assert len(scores) == len(idxs), scores
    out = {}
    for k, (gene, strand, i) in enumerate(zip(genes, strands, idxs)):
        # Per gene, each model is called for the reference, then the alternative.
        y = np.mean([r.outputs[2 * k] for r in recorders], axis=0)
        if strand == "-":
            y = y[:, ::-1]
        out[gene] = (strand, int(ann.tx_starts[i]), int(ann.tx_ends[i]), y[0, :, 1:3])
    return out


def compare(pred, chrom, pos, label, expect_genes=1):
    ref = pred.annotator.ref_fasta[chrom]
    n = len(ref)
    got = reference_inside_get_delta_scores(pred, chrom, pos)
    worst = 0.0
    for gene, (strand, tx_s, tx_e, y_ref) in got.items():
        y = pred.predict(lambda a, b: ref[a - 1:b].seq, n, strand, tx_s, tx_e, pos - 50, pos + 50)
        worst = max(worst, float(np.abs(y - y_ref).max()))
    check(f"{label}: predict() == get_delta_scores' reference prediction ({', '.join(got)})",
          len(got) == expect_genes and worst <= TOL, f"max |diff| {worst:.2g}")
    return worst


hg38 = predict.Predictor(str(helpers.HG38), "hg38")
ann = hg38.annotator
genes = list(ann.genes)


def gene_span(name):
    i = genes.index(name)
    return "chr" + str(ann.chroms[i]), int(ann.tx_starts[i]), int(ann.tx_ends[i]), ann.strands[i]


compare(hg38, "chr22", 28725242, "minus strand, inside CHEK2 (c.444+1)")
compare(hg38, "chr13", 32370955, "plus strand, inside BRCA2 (c.8488-1)")
for name in ("BRCA2", "CHEK2"):
    chrom, s, e, strand = gene_span(name)
    compare(hg38, chrom, s + 10, f"{strand} strand, window across {name}'s start (N before it)")
    compare(hg38, chrom, e - 10, f"{strand} strand, window across {name}'s end (N after it)")

# FLCN and the lncRNA RP11-45M22.4 overlap (the same pair test_pipeline.py uses).
c1, s1, e1, _ = gene_span("FLCN")
_, s2, e2, _ = gene_span("RP11-45M22.4")
compare(hg38, c1, (max(s1, s2) + min(e1, e2)) // 2, "a position in two overlapping genes, each with its own mask",
        expect_genes=2)

hg19 = predict.Predictor(str(helpers.HG19), "hg19")
compare(hg19, "chr10", 27326999, "ANKRD26 chr10:27326999 (hg19 regression variant)")

# --- the database's own mask, where it gives the VCF tab's input -------------------
db_genes = {g.name: g for g in transcripts.load_genes(
    "hg38", helpers.THIRDPARTY / "refseq", helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")}
for name, pos in (("CHEK2", 28725242), ("BRCA2", 32370955)):
    chrom, *_ = gene_span(name)
    ref = ann.ref_fasta[chrom]
    plan = refdb.plan_gene(db_genes[name], len(ref))
    lo, hi = plan["vcf_comparable"]
    inside = lo <= pos - 50 and pos + 50 <= hi
    (_strand, _s, _e, y_ref), = reference_inside_get_delta_scores(hg38, chrom, pos).values()
    y = hg38.predict(lambda a, b: ref[a - 1:b].seq, len(ref), plan["strand"], *plan["mask"], pos - 50, pos + 50)
    worst = float(np.abs(y - y_ref).max())
    check(f"{name}: with the database's extended mask {plan['mask']}, positions in its VCF-comparable stretch "
          f"equal get_delta_scores", inside and worst <= TOL, f"comparable {plan['vcf_comparable']}, max |diff| {worst:.2g}")

# Inside the 5,000 bases where CHEK2's mask moved, the values must differ --
# otherwise the mask isn't being used. Not right at the edge: the outermost few
# hundred bases of the 5,000 change the output by less than float32 resolves
# (measured: identical up to 100 bases in, 3e-14 at 500, 7e-9 at 1,000).
chrom, s, e, strand = gene_span("CHEK2")
ref = ann.ref_fasta[chrom]
plan = refdb.plan_gene(db_genes["CHEK2"], len(ref))
p = plan["vcf_comparable"][1] + 2000
a = hg38.predict(lambda x, y: ref[x - 1:y].seq, len(ref), strand, s, e, p - 50, p + 50)
b = hg38.predict(lambda x, y: ref[x - 1:y].seq, len(ref), strand, *plan["mask"], p - 50, p + 50)
check("...and outside it (the moved mask within reach) they differ", float(np.abs(a - b).max()) > 0,
      f"max |diff| {float(np.abs(a - b).max()):.2g}")

# --- chunk size doesn't change the result ----------------------------------------------
chrom, s, e, strand = gene_span("BRCA2")
ref = ann.ref_fasta[chrom]
region = (s + 20_000, s + 50_000)
whole = hg38.predict(lambda a, b: ref[a - 1:b].seq, len(ref), strand, s, e, *region, chunk=30_001)
for chunk in (7_777, 5_000, 1_000):
    part = hg38.predict(lambda a, b: ref[a - 1:b].seq, len(ref), strand, s, e, *region, chunk=chunk)
    check(f"chunks of {chunk:,} give the same 30 kb as one call", float(np.abs(part - whole).max()) <= TOL,
          f"max |diff| {float(np.abs(part - whole).max()):.2g}")

helpers.finish()
