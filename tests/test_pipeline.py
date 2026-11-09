"""Pipeline checks: the six fixes of 2026-09-11, the ANKRD26 known value and the
ClinVar demo VCF against its frozen output.

    .venv\\Scripts\\python tests\\test_pipeline.py
"""
import csv
import gzip
import os
import random
import sys

import helpers
from helpers import check, skip

helpers.isolate_home()   # the settings (e.g. the language) are never the real ones
helpers.watchdog(900)
helpers.setup()

from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline.cli import ReferenceMismatchError, run_pipeline_core  # noqa: E402
from spliceai_pipeline.merge import merge_variant  # noqa: E402
from spliceai_pipeline.normalize import normalize_variant, resolve_contig  # noqa: E402
from spliceai_pipeline.snpeff import _gene_annotations  # noqa: E402
from spliceai_pipeline.vcfio import count_vcf_records, is_gzip, read_vcf_records  # noqa: E402
from spliceai_pipeline.writer import FIELDNAMES, write_rows  # noqa: E402


def quiet(*_a, **_k):
    pass


hg19 = Fasta(str(helpers.HG19), rebuild=False)
hg38 = Fasta(str(helpers.HG38), rebuild=False)

# --- complex indels: normalization keeps the same change ---------------------------
random.seed(7)
bad = []
for _ in range(3000):
    alphabet = "ACGT" if random.random() < 0.5 else "AT"  # "AT": lots of repeats
    seq = "".join(random.choice(alphabet) for _ in range(60))

    def fetch(_c, s, e, seq=seq):
        return seq[s - 1:e]

    p = random.randint(10, 50)
    ref = seq[p - 1:p - 1 + random.randint(1, 5)]
    alt = "".join(random.choice("ACGT") for _ in range(random.randint(1, 5)))
    if ref == alt:
        continue
    q, r, a = normalize_variant(fetch, "c", p, ref, alt)

    def hap(P, R, A, seq=seq):
        return seq[:P - 1] + A + seq[P - 1 + len(R):]

    ok = seq[q - 1:q - 1 + len(r)] == r and hap(p, ref, alt) == hap(q, r, a)
    ok = ok and not (len(r) > 1 and len(a) > 1 and (r[0] == a[0] or r[-1] == a[-1]))  # minimal
    ok = ok and not (len(r) != len(a) and q > 1 and r[-1] == a[-1])                   # left-aligned
    if not ok:
        bad.append((seq, p, ref, alt, q, r, a))
check("normalization: 3000 random variants keep their sequence, minimal and left-aligned",
      not bad, f"{len(bad)} bad, e.g. {bad[:1]}")
seq = "N" * 99 + "XAT" + "GGGG"
check("normalization: complex indel AT>GCA stays itself",
      normalize_variant(lambda c, s, e: seq[s - 1:e], "c", 101, "AT", "GCA") == (101, "AT", "GCA"))

# --- blank lines --------------------------------------------------------------------
p = helpers.write_vcf(["chr10\t27326999\t.\tT\tC\t.\t.\t.\n", "\n", "   \n"])
try:
    n = len(list(read_vcf_records(p)))
    check("blank / whitespace-only lines are ignored", n == 1, f"{n} record(s)")
except Exception as exc:
    check("blank / whitespace-only lines are ignored", False, repr(exc))

# --- a VCF saved with a byte-order mark (some Windows editors) ----------------------
p = helpers.write_vcf(["chr10\t27326999\t.\tT\tC\t.\t.\t.\n"])
with open(p, "rb") as fh:
    raw = fh.read()
with open(p, "wb") as fh:
    fh.write(b"\xef\xbb\xbf" + raw)
try:
    n = len(list(read_vcf_records(p)))
    check("a VCF saved with a byte-order mark reads normally", n == 1, f"{n} record(s)")
except Exception as exc:
    check("a VCF saved with a byte-order mark reads normally", False, repr(exc))

# --- gzip is recognised by content, not by name ---------------------------------------
plain = helpers.write_vcf(["chr10\t27326999\t.\tT\tC\t.\t.\t.\n", "chr10\t27327036\t.\tA\tG\t.\t.\t.\n"])
with open(plain, "rb") as fh:
    raw = fh.read()
for name, data in (("bgzipped, named .vcf", gzip.compress(raw)), ("plain text, named .vcf.gz", raw)):
    path = helpers.write_vcf([]).replace(".vcf", ".vcf.gz" if "named .vcf.gz" in name else "_bgz.vcf")
    with open(path, "wb") as fh:
        fh.write(data)
    try:
        n = len(list(read_vcf_records(path)))
        check(f"a VCF {name} reads normally (the CLI used to crash on it)", n == 2 and count_vcf_records(path) == 2,
              f"{n} record(s)")
    except Exception as exc:
        check(f"a VCF {name} reads normally (the CLI used to crash on it)", False, repr(exc)[:120])

# --- chromosome names ---------------------------------------------------------------
contigs = set(hg19.keys())
check("MT / chrMT / M resolve to hg19's chrM",
      all(resolve_contig(c, contigs) == "chrM" for c in ("MT", "chrMT", "M", "chrM")))
check("'10' resolves to chr10; an unknown contig to None",
      resolve_contig("10", contigs) == "chr10" and resolve_contig("chrFake", contigs) is None)

# --- one run: known score, unknown chromosome, too-long deletion, MT indel ---------
long_ref = hg19["chr10"][27326950 - 1:27326950 - 1 + 151].seq.upper()
mt_ref = hg19["chrM"][99:101].seq.upper()
vcf = helpers.write_vcf([
    "chr10\t27326999\t.\tT\tC\t.\t.\t.\n",
    "chrFake\t100\t.\tA\tG\t.\t.\t.\n",
    f"chr10\t27326950\t.\t{long_ref}\t{long_ref[0]}\t.\t.\t.\n",
    f"MT\t100\t.\t{mt_ref}\t{mt_ref[0]}\t.\t.\t.\n",
    "\n",
])
rows = run_pipeline_core(vcf, "hg19", "masked", str(helpers.HG19), on_progress=quiet)
ank = [r for r in rows if r.pos == 27326999 and r.gene == "ANKRD26"]
check("ANKRD26 chr10:27326999 T>C scores DS_AG 0.29 (hg19)",
      ank and ank[0].DS_AG == 0.29 and ank[0].DP_AG == 12 and ank[0].source == "live",
      str([(r.gene, r.DS_AG, r.DP_AG, r.source) for r in rows if r.pos == 27326999]))
check("unknown chromosome: the run continues, the row says why",
      any(r.chrom == "chrFake" and r.source == "skipped: chromosome not in reference FASTA" for r in rows))
longdel = [r for r in rows if len(r.ref) > 100]
check("SpliceAI's own skip (REF > 100 bp) is shown, not a silent blank",
      longdel and longdel[0].source == "skipped: REF allele too long for SpliceAI",
      str([(r.gene, r.source) for r in longdel]))
check("an MT indel doesn't crash the run (mapped to chrM)", any(r.chrom == "chrM" for r in rows))

# --- wrong build: hg38 variants scored as hg19 stop before scoring -----------------
random.seed(11)
hg38_only = helpers.snv_lines(hg38, "chr1", sorted(random.sample(range(10_000_000, 240_000_000), 80)))
try:
    run_pipeline_core(helpers.write_vcf(hg38_only), "hg19", "masked", str(helpers.HG19), on_progress=quiet)
    check("wrong build stops with a clear message", False, "no error raised")
except ReferenceMismatchError as exc:
    check("wrong build stops with a clear message", "wrong genome build" in str(exc), str(exc)[:120])

# --- a FASTA whose chromosome names match none of the VCF's -----------------------
fa_path = helpers.write_vcf([]).replace(".vcf", ".fa")
with open(fa_path, "w") as fh:
    fh.write(">NC_000010.10\n" + "ACGT" * 50 + "\n")
Fasta(fa_path)  # builds the .fai
try:
    run_pipeline_core(helpers.write_vcf(["chr10\t5\t.\tA\tG\t.\t.\t.\n"]), "hg19", "masked", fa_path,
                      on_progress=quiet)
    check("RefSeq-named FASTA stops with a clear message", False, "no error raised")
except ReferenceMismatchError as exc:
    check("RefSeq-named FASTA stops with a clear message", "names them differently" in str(exc), str(exc)[:120])

# --- SnpEff annotation per gene --------------------------------------------------------
ann = ("ANN=G|intron_variant|MODIFIER|FLCN|FLCN|transcript|NM_144997.7|protein_coding|1/13|c.-24-394A>G|||||,"
       "G|upstream_gene_variant|MODIFIER|OTHER|OTHER|transcript|NM_000001.1|protein_coding||c.-500A>G|||||")
per_gene = _gene_annotations(ann, {})
anns = ["G|FLCN|0.48|0.00|0.00|0.00|1|2|3|4", "G|OTHER|0.10|0.00|0.00|0.00|1|2|3|4",
        "G|RP11-45M22.4|0.48|0.00|0.00|0.00|1|2|3|4"]
merged = {r.gene: r.snpeff_transcript for r in merge_variant("chr17", 1, "A", "G", None, anns, snpeff_ann=per_gene)}
check("SnpEff: each gene row gets its own gene's transcript; an unannotated gene stays blank",
      merged == {"FLCN": "NM_144997.7", "OTHER": "NM_000001.1", "RP11-45M22.4": None}, str(merged))
single = merge_variant("chr17", 1, "A", "G", None, ["G|RENAMED|0.1|0|0|0|1|2|3|4"], snpeff_ann=per_gene)
check("SnpEff: a single-gene variant under another symbol keeps SnpEff's overall pick",
      single[0].snpeff_transcript == "NM_144997.7")

if helpers.SNPEFF_DIR is None:
    skip("SnpEff per gene, for real", "no installer/thirdparty/snpeff -- run installer/fetch-thirdparty.ps1")
else:
    # A position where FLCN and the lncRNA RP11-45M22.4 overlap, from SpliceAI's own gene table.
    spans = {}
    with open(helpers.SPLICEAI_PKG / "spliceai" / "annotations" / "grch38.txt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] in ("FLCN", "RP11-45M22.4"):
                spans[f[0]] = (f[1], int(f[3]), int(f[4]))
    chrom = spans["FLCN"][0] if spans["FLCN"][0].startswith("chr") else "chr" + spans["FLCN"][0]
    lo, hi = max(s[1] for s in spans.values()), min(s[2] for s in spans.values())
    variant = helpers.snv_lines(hg38, chrom, range((lo + hi) // 2, hi))[:1]
    rows = run_pipeline_core(helpers.write_vcf(variant, ["##reference=hg38"]), "hg38", "masked", str(helpers.HG38),
                             use_snpeff=True, snpeff_dir=str(helpers.SNPEFF_DIR),
                             mane_dir=str(helpers.MANE_DIR) if helpers.MANE_DIR else None, on_progress=quiet)
    tx = {r.gene: r.snpeff_transcript for r in rows}
    check("SnpEff per gene, for real: FLCN row has its NM_ transcript, RP11-45M22.4 row doesn't borrow it",
          (tx.get("FLCN") or "").startswith("NM_") and "RP11-45M22.4" in tx and not tx["RP11-45M22.4"], str(tx))

# --- the ClinVar demo VCF against its output from version 1.1.1 -------------------
# tests/data/clinvar-splice-demo-hg38.{masked,raw}.tsv were written by the CLI of
# 1.1.1 (commit a38a6dd) before any change for the reference score database:
# masked with SnpEff and MANE, raw likewise. Raw is rerun without SnpEff (it adds
# nothing to the scores), so only its score columns are compared.
FROZEN_COLUMNS = [c for c in FIELDNAMES if c != "zygosity"]   # added in 1.2.0
SCORE_COLUMNS = [c for c in FROZEN_COLUMNS if not c.startswith("snpeff_")]
for mode in ("masked", "raw"):
    with_snpeff = mode == "masked" and helpers.SNPEFF_DIR is not None
    rows = run_pipeline_core(str(helpers.REPO / "examples" / "clinvar-splice-demo-hg38.vcf"), "hg38", mode,
                             str(helpers.HG38), use_snpeff=with_snpeff,
                             snpeff_dir=str(helpers.SNPEFF_DIR) if with_snpeff else None,
                             mane_dir=str(helpers.MANE_DIR) if with_snpeff and helpers.MANE_DIR else None,
                             on_progress=quiet)
    out = helpers.write_vcf([]).replace(".vcf", ".tsv")
    write_rows(rows, out)
    columns = FROZEN_COLUMNS if with_snpeff else SCORE_COLUMNS

    def table(path, columns=columns):
        with open(path, newline="") as fh:
            return [[r[c] for c in columns] for r in csv.DictReader(fh, delimiter="\t")]

    got, frozen = table(out), table(helpers.REPO / "tests" / "data" / f"clinvar-splice-demo-hg38.{mode}.tsv")
    check(f"ClinVar demo, {mode}: same output as 1.1.1" + (" (with SnpEff/MANE)" if with_snpeff else " (scores)"),
          got == frozen, f"{sum(a != b for a, b in zip(got, frozen))} of {len(frozen)} rows differ" if got != frozen else
          ", ".join(f"{r[4]} {r[13]}" for r in got))
    if mode == "masked" and not with_snpeff:
        skip("ClinVar demo, masked: SnpEff columns", "no installer/thirdparty/snpeff")

# --- zygosity from GT, and scoring only chosen genes (1.2.0) -----------------------
from spliceai_pipeline.transcripts import GeneRegions  # noqa: E402
from spliceai_pipeline.vcfio import zygosities  # noqa: E402
from spliceai_pipeline.writer import read_rows  # noqa: E402

cases = {("GT:AD", "0/1:5,5", 1): ["het"], ("GT", "1|1", 1): ["hom"], ("GT", "1", 1): ["hemi"],
         ("GT", "0/0", 1): [None], ("GT", "./.", 1): [None], ("AD", "5,5", 1): [None], ("", "", 1): [None],
         ("GT", "1/2", 2): ["het", "het"], ("GT", "2/2", 2): [None, "hom"], ("DP:GT", "30:0|1", 1): ["het"]}
check("zygosity: het / hom / hemi / none read from GT, per ALT",
      all(zygosities(f, v, n) == want for (f, v, n), want in cases.items()))

regions = GeneRegions("hg38", ["brca2", "CFTR", "NOTAGENE"])
check("chosen genes: SpliceAI's spelling, unknown names listed",
      regions.genes == {"BRCA2", "CFTR"} and regions.unknown == ["NOTAGENE"], f"{regions.genes} {regions.unknown}")
check("chosen genes: inside BRCA2 (chr13 and 13), not ABCA4's position, the margin widens",
      regions.contains("chr13", 32370955) and regions.contains("13", 32370955)
      and not regions.contains("chr1", 94011395) and not regions.contains("chr13", 32315507 - 50)
      and regions.contains("chr13", 32315507 - 50, margin=100))

demo = [l for l in open(helpers.REPO / "examples" / "clinvar-splice-demo-hg38.vcf") if not l.startswith("#")]
gts = ["0/1", "1/1", "0|1", "1/1", "0/1", "1", "./.", "0/1"]
with_gt = [l.rstrip("\n") + "\tGT\t" + g + "\n" for l, g in zip(demo, gts)]
header = ["##reference=hg38", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"]
path = helpers.write_vcf([]).replace(".vcf", "_s.vcf")
with open(path, "w") as fh:  # helpers.VCF_COLUMNS has no FORMAT/sample columns
    fh.write("".join(["##fileformat=VCFv4.2\n"] + [h + "\n" for h in header]
                     + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"] + with_gt))
messages = []
rows = run_pipeline_core(path, "hg38", "masked", str(helpers.HG38), genes=["APC", "BRCA2", "CHEK2"],
                         on_progress=lambda m, *a: messages.append(m))
with open(helpers.REPO / "tests" / "data" / "clinvar-splice-demo-hg38.masked.tsv", newline="") as fh:
    frozen = {r["gene"]: r for r in csv.DictReader(fh, delimiter="\t")}
check("chosen genes: only their rows (APC's neighbour CTC-554D6.1 left out)",
      [r.gene for r in rows] == ["APC", "BRCA2", "CHEK2"], str([r.gene for r in rows]))
check("chosen genes: the same scores as scoring everything",
      all(str(r.max_score) == frozen[r.gene]["max_score"] and str(r.DS_AG) == frozen[r.gene]["DS_AG"]
          and str(r.DS_DL) == frozen[r.gene]["DS_DL"] for r in rows))
check("chosen genes: only 3 of 8 variants normalized and scored",
      any("3 variant allele(s) inside the chosen genes" in m for m in messages))
check("zygosity on the rows: APC 1/1 hom, BRCA2 ./. blank, CHEK2 0/1 het",
      [r.zygosity for r in rows] == ["hom", None, "het"], str([r.zygosity for r in rows]))
out = path.replace(".vcf", ".tsv")
write_rows(rows, out)
check("zygosity: written and read back", [r.zygosity for r in read_rows(out)] == ["hom", None, "het"])
try:
    run_pipeline_core(path, "hg38", "masked", str(helpers.HG38), genes=["NOTAGENE"], on_progress=quiet)
    check("chosen genes: none known stops the run", False)
except ValueError as exc:
    check("chosen genes: none known stops the run", "None of the chosen genes" in str(exc))

# --- the installed program's "SpliceAI missing" message: no pip command -----------
import spliceai_gui.worker as worker  # noqa: E402
sys.frozen = True
try:
    msg = worker._spliceai_missing_message()
finally:
    del sys.frozen
check("installed program's 'SpliceAI missing' message points to the in-app button, no pip",
      "Download and install" in msg and "pip" not in msg)

# --- optional: your own VCF against its known-good output (never committed) -------
# The records scored are those at the positions the known-good file lists, so a
# small known-good file picks the variants that matter out of a whole exome.
# --- three single regression variants: ClinVar demo variants (CHEK2, MYBPC3 in hg19; BRCA2 in hg38) ------
# Build, position, REF, ALT and the expected unrounded masked delta scores (DS_AG, DS_AL, DS_DG,
# DS_DL), checked as numbers within the spec 5.6 tolerance with the local FASTA: no VCF needed.
REGRESSION_VARIANTS = [
    ("hg19", "chr22", 29121230, "C", "T", (0.000000, 0.000000, 0.940884, 0.997426)),
    ("hg19", "chr11", 47364709, "C", "T", (0.982925, 0.403659, 0.003322, 0.000000)),
    ("hg38", "chr13", 32370955, "G", "A", (0.368275, 0.920486, 0.000307, 0.000021)),
]
from spliceai_pipeline.refdb_validate import TOLERANCE, delta_scores  # noqa: E402
from spliceai_pipeline.score import LiveScorer  # noqa: E402

regression_annotators, bad = {}, []
for build, chrom, pos, ref, alt, want in REGRESSION_VARIANTS:
    if build not in regression_annotators:
        regression_annotators[build] = LiveScorer(str(helpers.HG19 if build == "hg19" else helpers.HG38), build,
                                                  "masked")._annotator
    _, live = delta_scores(regression_annotators[build], chrom, pos, ref, alt, masked=True)
    if not any(all(abs(a - b) <= TOLERANCE for a, b in zip(g["ds"], want)) for g in live.values()):
        bad.append(f"{build} {chrom}:{pos} {ref}>{alt}: {[[round(x, 6) for x in g['ds']] for g in live.values()]}")
check("three single regression variants: their unrounded masked delta scores, as numbers (spec 5.6 tolerance)",
      not bad, "; ".join(bad))

own_vcf, own_expected = os.environ.get("SPLICEAI_TEST_VCF"), os.environ.get("SPLICEAI_TEST_EXPECTED")
if not (own_vcf and own_expected):
    skip("own VCF vs known-good output", "set SPLICEAI_TEST_VCF and SPLICEAI_TEST_EXPECTED to use yours")
else:
    build = os.environ.get("SPLICEAI_TEST_BUILD", "hg38")
    with (gzip.open if is_gzip(own_vcf) else open)(own_vcf, "rt", encoding="utf-8-sig", errors="replace") as fh:
        lines = fh.readlines()
    with open(own_expected, newline="") as fh:
        wanted = {(r["chrom"].removeprefix("chr"), r["pos"]) for r in csv.DictReader(fh, delimiter="\t")}
    wanted |= {("chr" + c, p) for c, p in wanted}
    head = [ln.rstrip("\n") for ln in lines if ln.startswith("##")]
    subset = [ln for ln in lines if not ln.startswith("#") and tuple(ln.split("\t", 2)[:2]) in wanted]
    rows = run_pipeline_core(helpers.write_vcf(subset, head), build, "masked",
                             str(helpers.HG19 if build == "hg19" else helpers.HG38),
                             use_snpeff=helpers.SNPEFF_DIR is not None,
                             snpeff_dir=str(helpers.SNPEFF_DIR) if helpers.SNPEFF_DIR else None,
                             mane_dir=str(helpers.MANE_DIR) if helpers.MANE_DIR else None, on_progress=quiet)
    with open(own_expected, newline="") as fh:
        expected = {(r["chrom"], r["pos"], r["ref"], r["alt"], r["gene"]): r for r in csv.DictReader(fh, delimiter="\t")}
    # Compared as numbers (spec 5.6 rule), not as printed strings: a score next
    # to a rounding boundary (0.995332 prints as 1.00) can't flip the check.
    # The scores are recomputed unrounded from inside get_delta_scores. A
    # known-good file with 6-decimal (unrounded) values is held to TOLERANCE;
    # one with printed 2-decimal values gets half a print step more.
    from spliceai_pipeline.refdb_validate import TOLERANCE, compare, delta_scores
    from spliceai_pipeline.score import LiveScorer

    annotator = LiveScorer(str(helpers.HG19 if build == "hg19" else helpers.HG38), build, "masked")._annotator
    ds_cols, dp_cols = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"], ["DP_AG", "DP_AL", "DP_DG", "DP_DL"]
    got = {(r.chrom, str(r.pos), r.ref, r.alt, r.gene): r for r in rows}
    bad, live_by_variant = [], {}
    for key, e in expected.items():
        r = got.get(key)
        if r is None:
            bad.append(f"{key[4]} {key[0]}:{key[1]}: not in the output")
            continue
        if e.get("source") and r.source != e["source"]:
            bad.append(f"{key[4]} {key[0]}:{key[1]}: source {r.source}, expected {e['source']}")
        if e["DS_AG"] == "":
            continue
        variant = key[:4]
        if variant not in live_by_variant:
            live_by_variant[variant] = delta_scores(annotator, r.chrom, r.pos, r.ref, r.alt, masked=True)[1]
        live = live_by_variant[variant][r.gene]
        unrounded = all(len(e[c].partition(".")[2]) > 2 for c in ds_cols)
        problems, _flips = compare(live["ds"], live["dp"], live["series"], [float(e[c]) for c in ds_cols],
                                   [int(e[c]) for c in dp_cols], TOLERANCE if unrounded else TOLERANCE + 0.005)
        if [getattr(r, c) for c in ds_cols] != [round(v, 2) for v in live["ds"]]:
            problems.append("the printed scores aren't the unrounded ones rounded")
        if problems:
            bad.append(f"{key[4]} {key[0]}:{key[1]}: {'; '.join(problems)}")
    unexpected = sorted(set(got) - set(expected))
    check(f"own VCF ({len(subset)} records at the known-good file's positions) matches its known-good output, "
          "as numbers", subset and not bad and not unexpected,
          f"{len(rows)} rows, {len(bad)} differ, {len(unexpected)} not in the known-good file; {bad[:3]}")

helpers.finish()
