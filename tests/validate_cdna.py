"""Phase 8 validation of the Sequence-to-cDNA engine (spec 11, Phase 8 items
2-5), with the default thresholds -- nothing is tuned to make a control
pass. Not a pass/fail check: it writes a report.

1. Wild-type concordance on 30 MANE genes (hg38) of different sizes and
   exon counts: exon-level identity of the WT top transcript with the
   canonical one, per gene and overall, and the time per gene.
2. Literature positive controls with RNA-proven outcomes, entered as the
   patient's sequence: the HGVS c. position is resolved with the local
   annotation and the reference base checked against the FASTA (a mismatch
   stops that case); whether the documented outcome appears, its rank and
   relative score.
3. Synthetic controls on several genes: a random 65 bp insertion deep in an
   intron, an exon with its flanks inserted into an intron, an 80 bp deletion
   removing an annotated acceptor, an inversion containing an exon.
4. (Local only) the patient variants of docs/local_regression.md, when
   SPLICEAI_PATIENT_TSV names their file(s) (gene, chrom, pos, ref, alt and
   build per line) and SPLICEAI_PATIENT_REPORT where to write; never into the
   repository's tracked files.

Needs the Phase 2 databases in ~/SpliceAI_reference_data/spliceai-refdb.

    .venv\\Scripts\\python tests\\validate_cdna.py docs\\phase_reports\\phase_8_validation.md
"""
import os
import random
import sys
import time

import helpers

helpers.setup()

from spliceai_pipeline import cdna, cdna_analysis as ca, haplotype, refdb  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.score import _Record, get_delta_scores  # noqa: E402

REFSEQ, MANE = str(helpers.THIRDPARTY / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
ROOT = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb")
rng = random.Random(8)
out = []
say = lambda line="": (out.append(line), print(line, flush=True))


class Build:
    def __init__(self, build, fasta):
        self.build = build
        self.reader = refdb.Reader(os.path.join(ROOT, build))
        self.scores = refdb.ReferenceScores(build, str(fasta), reader=self.reader, refseq_dir=REFSEQ,
                                            mane_dir=str(helpers.MANE_DIR))
        self.fasta = self.scores.predictor().annotator.ref_fasta
        self.table = ca.GeneTable(build, REFSEQ, MANE)

    def record(self, gene):
        rec = self.fasta["chr" + gene.chrom]
        return rec, (lambda a, b, rec=rec: rec[a - 1:b].seq)

    def analyse(self, gene, edits, apply_all=True):
        rec, fetch = self.record(gene)
        region = (max(1, gene.mask[0] - 5000), min(len(rec), gene.mask[1] + 5000))
        ref = fetch(*region).upper()
        patient, _ = haplotype.build(ref, region[0], edits)
        compat = ca.check_compatibility(patient, gene, fetch, len(rec))
        for e in compat.alignment.events:
            e.apply = e.apply or apply_all
        return ca.run_analysis(self.build, gene, compat, self.scores, fetch, self.table)


def transcripts(result, kind):
    return [t for t in result.transcripts if t.kind == kind]


def ref_exons_of(g):
    c = g.canonical
    return list(c.exons) if c.strand == "+" else list(reversed(c.exons))


def extra_exons(t, base):
    """Exons of t (reference coordinates) overlapping no exon of base: [(start, end)]."""
    return [e for e in t.ref_exons if e and not any(e[0] <= b[1] and b[0] <= e[1] for b in base)]


def spliceai_ds(gene, pos, ref, alt):
    """SpliceAI's own get_delta_scores for the variant (VCF form, unmasked): 'AG/AL/DG/DL'."""
    out = get_delta_scores(_Record("chr" + gene.chrom, pos, ref, alt), hg38.scores.predictor().annotator, 50, 0)
    return next(("/".join(f.split("|")[2:6]) for f in out if f.split("|")[1] == gene.name), "--")


hg38 = Build("hg38", helpers.HG38)

# --- 1. wild-type concordance on 30 MANE genes ----------------------------------------------------
say("## 1. Wild-type concordance, 30 MANE genes (hg38)")
say()
mane = [g for g in hg38.table.genes if g.rule == "MANE Select" and g.canonical and g.canonical.cds
        and len(g.canonical.exons) >= 2 and g.chrom not in ("X", "Y") and hg38.reader.plan(g.name)]
bins = [(0, 10e3), (10e3, 30e3), (30e3, 100e3), (100e3, 300e3), (300e3, 1e6), (1e6, 3e6)]
chosen = []
for lo, hi in bins:
    pool = sorted((g for g in mane if lo <= g.mask[1] - g.mask[0] < hi), key=lambda g: (len(g.canonical.exons), g.name))
    picks = [pool[int((len(pool) - 1) * f)] for f in (0.1, 0.3, 0.5, 0.7, 0.9)] if len(pool) >= 5 else pool
    chosen += picks
say("| Gene | Strand | Size (kb) | Exons | Identical | Rate | Relative score of #1 | Time (s) |")
say("|---|---|---|---|---|---|---|---|")
total_same = total = 0
from spliceai_pipeline import splice_paths  # noqa: E402

for g in chosen:
    t0 = time.time()
    c = g.canonical
    wt, plan, _ = hg38.scores.gene(g.name)
    tss, tes = (c.start, c.end) if c.strand == "+" else (c.end, c.start)
    top = splice_paths.top_paths(wt[:, 0], wt[:, 1], plan["start"], c.strand, tss, tes)
    ex = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
    rows, _ = cdna.concordance(ex, top[0].exons, c.strand)
    same = sum(r["status"] == "identical" for r in rows if r["ref"])
    total_same, total = total_same + same, total + len(ex)
    seconds = time.time() - t0
    say(f"| {g.name} | {g.strand} | {(g.mask[1] - g.mask[0]) / 1e3:,.1f} | {len(ex)} | {same} | {same / len(ex):.0%} | "
        f"{top[0].relative:.3f} | {seconds:.2f} |")
say()
say(f"Overall: {total_same} of {total} canonical exons identical in the WT top transcript ({total_same / total:.1%}).")
say()

# --- 2. literature controls ------------------------------------------------------------------------
say("## 2. Literature positive controls (hg38, entered as the patient's sequence)")
say()
# (gene, HGVS c. on the MANE transcript, ref, alt, outcome kind, size or exon, source).
# The outcomes as reported in the literature cited (checked against ClinVar's
# summaries of it, 2026-10-27); the coordinates come from the local
# annotation, and the reference base is checked below.
LITERATURE = [
    ("CFTR", "3718-2477", "C", "T", "pseudoexon", 84,
     "Highsmith et al. 1994, N Engl J Med (3849+10kbC>T): an 84 bp cryptic exon"),
    ("CFTR", "1680-886", "A", "G", "pseudoexon", 49,
     "Chillon et al. 1995, Am J Hum Genet 56:623 (1811+1.6kbA>G): a 49 bp exon"),
    ("ATM", "5763-1050", "A", "G", "pseudoexon", 137,
     "McConville et al. 1996 (IVS40-1050A>G, reported as 5762ins137): a 137 nt insertion"),
    ("CFTR", "2657+5", "G", "A", "skipping", "2657",
     "Highsmith et al. 1997, Hum Mutat (2789+5G>A): skipping of exon 16 (legacy 14b)"),
]
say("| Variant | Reference base | Documented outcome | SpliceAI AG/AL/DG/DL | Found in patient top 5 | Rank | "
    "Relative score | Patient #1 vs WT |")
say("|---|---|---|---|---|---|---|---|")
for name, c_pos, ref_base, alt, kind, what, source in LITERATURE:
    g = hg38.table.find(name)[0]
    rec, fetch = hg38.record(g)
    ref_cdna = cdna.reference(lambda a, b: fetch(a, b).upper(), g.canonical)
    pos = cdna.c_to_genomic(ref_cdna, c_pos)
    base = fetch(pos, pos).upper()
    genomic_ref = ref_base if g.strand == "+" else {"A": "T", "C": "G", "G": "C", "T": "A"}[ref_base]
    genomic_alt = alt if g.strand == "+" else {"A": "T", "C": "G", "G": "C", "T": "A"}[alt]
    label = f"{name} {g.canonical.accession}:c.{c_pos}{ref_base}>{alt} (chr{g.chrom}:{pos:,})"
    if base != genomic_ref:
        say(f"| {label} | **mismatch: {base} at the position, expected {genomic_ref}; case stopped** | {kind} | | | | | |")
        continue
    result = hg38.analyse(g, [Edit(pos, pos, genomic_alt, "SNV")])
    wt1 = transcripts(result, "WT predicted")[0]
    found = None
    for t in transcripts(result, "Patient"):
        if kind == "pseudoexon":
            extra = extra_exons(t, [e for e in wt1.ref_exons if e])
            if any(abs((e[1] - e[0] + 1) - what) <= 5 for e in extra):
                found = t
                break
        else:
            # Skipping of the canonical exon whose last base is c.<what>.
            donor = cdna.c_to_genomic(ref_cdna, what)
            exon = next(k for k, e in enumerate(ref_exons_of(g), 1) if min(e) <= donor <= max(e))
            if any(r["ref"] == exon and r["status"] == "missing" for r in t.vs_ref[1]):
                found = t
                break
    p1 = transcripts(result, "Patient")[0]
    ds = spliceai_ds(g, pos, genomic_ref, genomic_alt)
    say(f"| {label} | {base} (verified) | {source} | {ds} | yes | {found.rank} | {found.relative:.3f} | "
        f"{p1.change_vs_wt} |" if found else
        f"| {label} | {base} (verified) | {source} | {ds} | no | | | {p1.change_vs_wt} |")
say()

# --- 3. synthetic controls --------------------------------------------------------------------------
say("## 3. Synthetic controls (hg38)")
say()
say("| Gene | Control | Expected | Patient #1 vs WT | Outcome | Rank of the expected | Relative score |")
say("|---|---|---|---|---|---|---|")
for name in ("BRCA2", "SAMSN1", "MRPL39", "USP25", "APP"):
    g = hg38.table.find(name)[0]
    rec, fetch = hg38.record(g)
    c = g.canonical
    ex = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
    introns = sorted(((ex[i], ex[i + 1]) for i in range(len(ex) - 1)), key=lambda p: -abs(p[1][0] - p[0][1]))
    a, b = sorted((introns[0][0][1] if c.strand == "+" else introns[0][0][0],
                   introns[0][1][0] if c.strand == "+" else introns[0][1][1]))
    deep = (a + b) // 2                                   # the middle of the largest intron
    k = len(ex) // 2                                      # an internal exon
    e_lo, e_hi = sorted(ex[k])
    flank = fetch(e_lo - 150, e_hi + 150).upper()
    acceptor = ex[k][0] if c.strand == "+" else ex[k][1]  # the exon's first base, in transcript order
    deletion = (acceptor - 60, acceptor + 19) if c.strand == "+" else (acceptor - 19, acceptor + 60)
    seq = "".join(rng.choice("ACGT") for _ in range(65))
    anchor = fetch(deep - 1, deep - 1).upper()
    controls = [
        (f"65 bp random insertion deep in an intron (SpliceAI AG/AL/DG/DL "
         f"{spliceai_ds(g, deep - 1, anchor, anchor + seq)})", [Edit(deep, deep - 1, seq, "insertion")],
         "no change in the top transcript"),
        (f"exon {k + 1} with 150 bp flanks inserted into an intron", [Edit(deep, deep - 1, flank, "insertion")],
         "a pseudoexon candidate"),
        (f"80 bp deletion removing exon {k + 1}'s acceptor", [Edit(*deletion, "", "deletion")],
         f"exon {k + 1} skipped or a cryptic acceptor"),
        (f"inversion containing exon {k + 1}", [Edit(e_lo - 100, e_hi + 100, "", "inversion")], "(behaviour reported)"),
    ]
    for label, edits, expected in controls:
        result = hg38.analyse(g, edits)
        wt1 = transcripts(result, "WT predicted")[0]
        patients = transcripts(result, "Patient")
        hit = None
        for t in patients:
            vs_wt = t.vs_wt[1]          # against the WT prediction
            vs_ref = t.vs_ref[1]        # against the canonical transcript: exon k + 1 is canonical numbering
            if expected.startswith("no change") and all(r["status"] == "identical" for r in vs_wt):
                hit = t
            elif expected.startswith("a pseudoexon") and any(r["status"] == "extra (novel)" for r in vs_wt):
                hit = t
            elif "skipped or" in expected and any(r["ref"] == k + 1 and (r["status"] == "missing" or "start shifted"
                                                                          in r["status"]) for r in vs_ref):
                hit = t
            if hit:
                break
        if expected.startswith("("):
            outcome = "--"
        elif hit is patients[0]:
            outcome = "as expected, top transcript"
        elif hit:
            outcome = f"as expected at rank {hit.rank}"
        else:
            outcome = "not among the 5"
        say(f"| {name} ({c.strand}) | {label} | {expected} | {patients[0].change_vs_wt} | {outcome} | "
            f"{hit.rank if hit else ''} | {f'{hit.relative:.3f}' if hit else ''} |")
say()

# --- 4. the patient variants (local only) ---------------------------------------------------------------
tsvs, report_path = os.environ.get("SPLICEAI_PATIENT_TSV"), os.environ.get("SPLICEAI_PATIENT_REPORT")
if tsvs and report_path:
    builds = {"hg38": hg38}
    lines = ["## Phase 8 item 4: the patient variants as patient sequences (local, never committed)", ""]
    for path in tsvs.split(os.pathsep):
        for row in open(path, encoding="utf-8"):
            if not row.strip() or row.startswith("#"):
                continue
            name, chrom, pos, ref_allele, alt, build = row.split()[:6]
            if build not in builds:
                builds[build] = Build(build, helpers.HG19 if build == "hg19" else helpers.HG38)
            bd = builds[build]
            g = bd.table.find(name)[0]
            rec, fetch = bd.record(g)
            pos = int(pos)
            if fetch(pos, pos + len(ref_allele) - 1).upper() != ref_allele.upper():
                lines.append(f"- {name}: the reference allele doesn't match {build} at {chrom}:{pos}; stopped.")
                continue
            e = Edit(pos, pos, alt, "SNV") if len(ref_allele) == len(alt) == 1 else \
                Edit(pos + 1, pos + len(ref_allele) - 1, alt[1:], "deletion" if len(alt) == 1 else "complex")
            result = bd.analyse(g, [e])
            lines.append(f"- {name} {build} {chrom}:{pos} {ref_allele}>{alt}:")
            for t in result.transcripts:
                if t.kind == "Patient" or (t.kind == "WT predicted" and t.rank == 1):
                    lines.append(f"  - {t.kind} {t.rank} (relative {t.relative:.3f}): vs reference: {t.change_vs_ref}; "
                                 f"vs WT: {t.change_vs_wt}")
    with open(report_path, "a", encoding="utf-8") as fh:
        fh.write("\n" + "\n".join(lines) + "\n")
    say("## 4. The patient variants")
    say()
    say(f"Run; the predicted transcripts are in the local, git-ignored file given by SPLICEAI_PATIENT_REPORT.")

if len(sys.argv) > 1:
    with open(sys.argv[1], "w", encoding="utf-8") as fh:
        fh.write("# Phase 8 validation of the Sequence-to-cDNA engine\n\n"
                 "Written by `tests/validate_cdna.py` (default thresholds, nothing tuned).\n\n" + "\n".join(out) + "\n")
