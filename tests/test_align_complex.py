"""Spec 14 C8: complex rearrangements in the aligner. Synthetic samples made
from the real hg38 reference -- at least 20 per type -- are aligned back:

- delins (a deletion plus novel bases), an inversion with flanking
  deletions, an inverted duplication, an insertion of an inverted copy of
  nearby sequence, a deletion with a templated insertion from elsewhere in
  the gene, a deletion with an inverted duplication next to it, a mosaic of
  two copies (a complex replacement), deletions with 1 to 10 bp of
  microhomology, and deletions with both breakpoints inside Alu repeats.

The events found must rebuild the sample exactly, every time (the haplotype
is always the paste; the classification is only a label). The share labelled
as made is reported, per type; a microhomology must be reported with its
length and left-aligned, an Alu breakpoint's range must hold the true one.

    .venv\\Scripts\\python tests\\test_align_complex.py [report.md]
"""
import random
import sys

import helpers
from helpers import check

helpers.setup()

from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline import align, cdna_analysis as ca, haplotype  # noqa: E402
from spliceai_pipeline.align import reverse_complement as rc  # noqa: E402

REFSEQ, MANE = str(helpers.THIRDPARTY / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
table = ca.GeneTable("hg38", REFSEQ, MANE)
fasta = Fasta(str(helpers.HG38))
rng = random.Random(14)
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
PER_TYPE = 20

genes = []
for name in ("SAMSN1", "BRCA2", "USP25", "APP"):
    g = table.find(name)[0]
    record = fasta["chr" + g.chrom]
    region = (g.mask[0] - 5000, g.mask[1] + 5000)
    genes.append((g, record, region, record[region[0] - 1:region[1]].seq.upper()))


def place(ref, span):
    """A random position with `span` bases of clean sequence after it and 12 kb from both region ends."""
    while True:
        p = rng.randint(12_000, len(ref) - 12_000 - span)
        if "N" not in ref[p - 6000:p + span + 6000]:
            return p


def make(kind, ref):
    """(sample, (lo, hi) of the change in ref, expected label, extra) for one case."""
    if kind == "delins":
        n, m = rng.randint(20, 500), rng.randint(10, 300)
        p = place(ref, n)
        return ref[:p] + rand(m) + ref[p + n:], (p, p + n), "delins (novel sequence)", None
    if kind == "del-inv-del":
        d1, i, d2 = rng.randint(10, 200), rng.randint(100, 2000), rng.randint(10, 200)
        p = place(ref, d1 + i + d2)
        return (ref[:p] + rc(ref[p + d1:p + d1 + i]) + ref[p + d1 + i + d2:], (p, p + d1 + i + d2),
                "deletion-inversion-deletion", None)
    if kind == "dup-inv":
        n = rng.randint(50, 1000)
        p = place(ref, 0)
        return ref[:p] + rc(ref[p - n:p]) + ref[p:], (p - n, p), "inverted duplication", None
    if kind == "ins-inv":
        n, d = rng.randint(50, 500), rng.randint(500, 5000)
        p = place(ref, 0)
        src = ref[p + d:p + d + n] if rng.random() < 0.5 else ref[p - d - n:p - d]
        return ref[:p] + rc(src) + ref[p:], (p, p), "inverted templated insertion", None
    if kind == "templated":
        n, m, d = rng.randint(10, 300), rng.randint(50, 500), rng.randint(1000, 5000)
        p = place(ref, n)
        src = ref[p + n + d:p + n + d + m] if rng.random() < 0.5 else ref[p - d - m:p - d]
        return ref[:p] + src + ref[p + n:], (p, p + n), "deletion with templated insertion", None
    if kind == "del-dup-inv":
        n, m = rng.randint(10, 300), rng.randint(50, 500)
        p = place(ref, n)
        return ref[:p] + rc(ref[p - m:p]) + ref[p + n:], (p - m, p + n), "deletion with inverted duplication", None
    if kind == "mosaic":
        n, a, b = rng.randint(10, 200), rng.randint(100, 300), rng.randint(100, 300)
        p = place(ref, n)
        s1, s2 = ref[p - 3000:p - 3000 + a], rc(ref[p + n + 2000:p + n + 2000 + b])
        return ref[:p] + s1 + s2 + ref[p + n:], (p, p + n), "complex replacement", 2
    if kind == "microhomology":
        k = 1 + len(made[kind]) % 10
        while True:
            s = place(ref, 3000)
            hit = ref.find(ref[s:s + k], s + 50, s + 3000)
            n = hit - s
            if hit > 0 and ref[s + k] != ref[s + n + k] and ref[s - 1] != ref[s + n - 1]:
                return ref[:s] + ref[s + n:], (s, s + n), "deletion", (k, s, n)
    if kind == "alu":
        # Two Alus in the same orientation, found by the core of the Alu left arm; the sample keeps the first
        # one's start and the second one's end (as an Alu-Alu recombination does), cut at the same place in both.
        pairs = []
        for motif in ("CCTGTAATCCCAGC", rc("CCTGTAATCCCAGC")):
            hits, at = [], ref.find(motif)
            while at >= 0:
                hits.append(at)
                at = ref.find(motif, at + 1)
            pairs += [(a, b, len(motif)) for a in hits for b in hits
                      if 1000 <= b - a <= 15_000 and 12_000 < a and b < len(ref) - 12_000]
        a, b, m = rng.choice(pairs)
        o = rng.randint(0, m - 1)
        return ref[:a + o] + ref[b + o:], (a + o, b + o), "deletion", (a + o, b - a)
    raise ValueError(kind)


TYPES = ["delins", "del-inv-del", "dup-inv", "ins-inv", "templated", "del-dup-inv", "mosaic", "microhomology", "alu"]
made = {t: [] for t in TYPES}
rows, failures = [], []
for kind in TYPES:
    exact = labelled = 0
    for i in range(PER_TYPE):
        g, record, region, ref = genes[i % len(genes)]
        sample, (lo, hi), expected, extra = make(kind, ref)
        made[kind].append(1)
        shift = len(sample) - len(ref)
        paste = sample[lo - 4000:hi + shift + 4000] if hi + shift > lo else sample[lo - 4000:lo + 4000]
        a = align.align(paste, ref, region[0], lambda c, x, y, r=record: r[x - 1:y].seq.upper(), "chr" + g.chrom)
        for e in a.events:
            e.apply = True
        try:
            rebuilt, _ = haplotype.apply_events(ref, region[0], a.events)
        except ValueError as exc:
            rebuilt = f"refused: {exc}"
        ok = rebuilt == sample
        exact += ok
        labels = sorted({e.label for e in a.events})
        good = labels == [expected]
        if good and kind == "microhomology":
            k, s, n = extra
            e = a.events[0]
            good = (len(a.events) == 1 and e.microhomology == ref[s:s + k] and e.start == region[0] + s
                    and e.end == region[0] + s + n - 1)
        elif good and kind == "alu":
            s, n = extra
            e = a.events[0]
            true = region[0] + s
            low, high = (e.breakpoints[0] if e.breakpoints else [e.start, e.start])
            good = len(a.events) == 1 and e.end - e.start + 1 == n and low <= true <= high
        elif good and kind == "mosaic":
            good = len(a.events) == 1 and len(a.events[0].blocks) == extra
        labelled += good
        if not ok or not good:
            failures.append(f"{kind} #{i + 1} ({g.name}): rebuilt {'exactly' if ok else 'WRONG'}; labels {labels}, "
                            f"expected {expected!r}; events {[(e.kind, e.start, e.end, len(e.inserted)) for e in a.events]}")
    rows.append((kind, PER_TYPE, exact, labelled))

total = sum(r[1] for r in rows)
check(f"every synthetic sample ({total}, {PER_TYPE} per type) is rebuilt exactly from the events found",
      all(r[2] == r[1] for r in rows), "; ".join(f for f in failures if "WRONG" in f)[:600])
report = ["| Type | Cases | Rebuilt exactly | Labelled as made | Accuracy |", "|---|---|---|---|---|"]
report += [f"| {k} | {n} | {x} | {lab} | {lab / n:.0%} |" for k, n, x, lab in rows]
print("\n".join(report))
for f in failures:
    print("  -", f)
check("labelled as made: at least 90% of each type", all(r[3] >= 0.9 * r[1] for r in rows),
      ", ".join(f"{k} {lab}/{n}" for k, n, _x, lab in rows))
if len(sys.argv) > 1:
    with open(sys.argv[1], "w", encoding="utf-8") as fh:
        fh.write("\n".join(report) + "\n\n" + "\n".join(f"- {f}" for f in failures) + "\n")

# Spec 17 item 6: in repeats the alignment finds no inverted seeds and gives an inversion as a replacement.
dice = random.Random(5)
stretch = "".join(dice.choice("ACGT") for _ in range(400))
piece = rc(stretch[100:300])
exact = align.Event("complex", 1100, stretch[100:300], piece, 1100, 1299, piece, 200, apply=True)
nearly = "A" if piece[50] != "A" else "C"
changed = align.Event("complex", 1100, stretch[100:300], piece, 1100, 1299, piece[:50] + nearly + piece[51:], 200, apply=True)
align.classify([exact], stretch, 1000)
align.classify([changed], stretch, 1000)
check("an exact inversion given as a replacement is an inversion; with a changed base it stays a replacement",
      (exact.kind, exact.label) == ("inversion", "inversion") and (changed.kind, changed.label) == ("complex", "inversion"),
      f"{exact.kind}/{exact.label}, {changed.kind}/{changed.label}")
helpers.finish()
