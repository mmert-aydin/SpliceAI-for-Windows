"""Pipeline checks: the six fixes of 2026-09-11 and the ANKRD26 known value.

    .venv\\Scripts\\python tests\\test_pipeline.py
"""
import csv
import os
import random
import sys

import helpers
from helpers import check, skip

helpers.setup()

from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline.cli import ReferenceMismatchError, run_pipeline_core  # noqa: E402
from spliceai_pipeline.merge import merge_variant  # noqa: E402
from spliceai_pipeline.normalize import normalize_variant, resolve_contig  # noqa: E402
from spliceai_pipeline.snpeff import _gene_annotations  # noqa: E402
from spliceai_pipeline.vcfio import read_vcf_records  # noqa: E402


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
own_vcf, own_expected = os.environ.get("SPLICEAI_TEST_VCF"), os.environ.get("SPLICEAI_TEST_EXPECTED")
if not (own_vcf and own_expected):
    skip("own VCF vs known-good output", "set SPLICEAI_TEST_VCF and SPLICEAI_TEST_EXPECTED to use yours")
else:
    build = os.environ.get("SPLICEAI_TEST_BUILD", "hg38")
    with open(own_vcf) as fh:
        lines = fh.readlines()
    head = [ln.rstrip("\n") for ln in lines if ln.startswith("##")]
    subset = [ln for ln in lines if not ln.startswith("#") and ln.strip()][:20]
    rows = run_pipeline_core(helpers.write_vcf(subset, head), build, "masked",
                             str(helpers.HG19 if build == "hg19" else helpers.HG38),
                             use_snpeff=helpers.SNPEFF_DIR is not None,
                             snpeff_dir=str(helpers.SNPEFF_DIR) if helpers.SNPEFF_DIR else None,
                             mane_dir=str(helpers.MANE_DIR) if helpers.MANE_DIR else None, on_progress=quiet)
    with open(own_expected, newline="") as fh:
        expected = {(r["chrom"], r["pos"], r["ref"], r["alt"], r["gene"]): r for r in csv.DictReader(fh, delimiter="\t")}
    cols = ["DS_AG", "DP_AG", "DS_AL", "DP_AL", "DS_DG", "DP_DG", "DS_DL", "DP_DL", "max_score", "source"]
    diffs = [r for r in rows
             if (e := expected.get((r.chrom, str(r.pos), r.ref, r.alt, r.gene))) is None
             or ["" if getattr(r, c) is None else str(getattr(r, c)) for c in cols] != [e[c] for c in cols]]
    check("own VCF (first 20 records) matches its known-good output", not diffs,
          f"{len(rows)} rows, {len(diffs)} differ")

helpers.finish()
