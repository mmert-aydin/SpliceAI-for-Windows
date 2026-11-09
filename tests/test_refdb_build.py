"""The reference score database builder (spec 5.2-5.5, Phase 1): transcript
matching and masks, coordinate conversions, what each gene stores, storage
precision, and a real build that is killed mid-way and then resumed.

    .venv\\Scripts\\python tests\\test_refdb_build.py
"""
import glob
import gzip
import json
import os
import subprocess
import sys
import tempfile
import time

import helpers
from helpers import check

helpers.watchdog(900)
helpers.setup()

import numpy as np  # noqa: E402

from spliceai_pipeline import predict, refdb, transcripts  # noqa: E402
from spliceai_pipeline.transcripts import Gene, Transcript  # noqa: E402

REFSEQ = helpers.THIRDPARTY / "refseq"
MANE = helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz"


def tx(acc, symbol, exons, chrom="chr1", strand="+", cds=None):
    return Transcript(acc, symbol, chrom, strand, exons[0][0], exons[-1][1], cds, tuple(exons))


# --- reading the tables: 0-based half-open -> 1-based inclusive ---------------------
tmp = tempfile.mkdtemp(prefix="spliceai_refdb_test_")
table = os.path.join(tmp, "t.txt.gz")
with gzip.open(table, "wt") as fh:
    fh.write("0\tNM_1.1\tchr1\t+\t100\t400\t150\t350\t2\t100,300,\t200,400,\t0\tAAA\tcmpl\tcmpl\t0,2,\n")
    fh.write("0\tNR_2.1\tchr1\t-\t100\t400\t400\t400\t1\t100,\t400,\t0\tBBB\tnone\tnone\t-1,\n")
    fh.write("0\tXM_3.1\tchr1\t+\t100\t400\t150\t350\t1\t100,\t400,\t0\tCCC\tcmpl\tcmpl\t0,\n")
    fh.write("0\tNM_4.1\tchr1_KI270706v1_random\t+\t100\t400\t150\t350\t1\t100,\t400,\t0\tDDD\tcmpl\tcmpl\t0,\n")
read = transcripts.read_table(table)
check("read_table keeps curated NM_/NR_ on chr1-22/X/Y only (no XM_, no alt contigs)",
      [t.accession for t in read] == ["NM_1.1", "NR_2.1"])
t = read[0]
check("read_table: exons, transcript and CDS become 1-based inclusive",
      (t.start, t.end, t.exons, t.cds) == (101, 400, ((101, 200), (301, 400)), (151, 350)), str(t))
check("read_table: CDS length is counted over exons only", t.cds_length == 100, str(t.cds_length))
check("read_table: a non-coding transcript has no CDS", read[1].cds is None and read[1].cds_length == 0)

ann = os.path.join(tmp, "a.txt")
with open(ann, "w") as fh:
    fh.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
    fh.write("GENE1\t1\t-\t99\t400\t99,299,\t200,400,\n")
g = transcripts.read_spliceai_genes(ann)[0]
check("read_spliceai_genes: TX_START and EXON_START are 0-based in SpliceAI's table",
      (g.start, g.end, g.exons, g.strand) == (100, 400, ((100, 200), (300, 400)), "-"), str(g))
check("exonic_overlap counts shared bases", transcripts.exonic_overlap(((1, 10), (21, 30)), ((5, 25),)) == 11)

# --- matching transcripts to SpliceAI's genes ------------------------------------------
genes = [
    # An old symbol: its MANE transcript has a new one, and starts 100 bases earlier.
    Gene("OLDNAME", "1", "+", 1001, 5000, ((1001, 1200), (3001, 3200), (4801, 5000))),
    # A read-through in SpliceAI's table, sharing exons with two MANE genes.
    Gene("RT", "1", "+", 20001, 30000, ((20001, 20100), (29901, 30000))),
    # No MANE: two own NM_, and a host gene's transcript spanning it in its intron.
    Gene("NOMANE", "2", "+", 1001, 3000, ((1001, 1100), (2901, 3000))),
    # Only a host gene's transcript, around it in an intron.
    Gene("NESTED", "3", "+", 5001, 6000, ((5001, 5100), (5901, 6000))),
]
table = [
    tx("NM_000001.1", "NEWNAME", [(901, 1200), (3001, 3200), (4801, 5000)]),
    tx("NM_000002.1", "NEWNAME", [(1001, 1200), (8801, 9000)]),                  # alternative, reaches 9000
    tx("NM_000003.1", "NEWNAME", [(3001, 3200), (12001, 12200)]),                # read-through into NEIGHBOR
    tx("NM_000004.1", "NEIGHBOR", [(12001, 12200), (13001, 13500)]),
    tx("NM_000005.1", "NEWNAME", [(1001, 1200), (6001, 7000)], strand="-"),      # other strand
    tx("NM_000010.1", "AAA", [(20001, 20100), (21001, 21100)]),
    tx("NM_000011.1", "BBB", [(29001, 29100), (29951, 30000)]),
    tx("NM_000020.1", "NOMANE", [(1001, 1100), (2901, 3000)], chrom="chr2", cds=(1001, 3000)),
    tx("NM_000021.1", "NOMANE", [(1001, 1100), (2901, 3000)], chrom="chr2", cds=(1050, 2950)),
    tx("NM_000022.1", "HOST", [(1, 100), (9901, 10000)], chrom="chr2", cds=(1, 10000)),
    tx("NM_000030.1", "HOST3", [(1, 100), (9901, 10000)], chrom="chr3", cds=(1, 10000)),
]
mane = {a: (a + ".1", s) for a, s in
        (("NM_000001", "NEWNAME"), ("NM_000004", "NEIGHBOR"), ("NM_000010", "AAA"), ("NM_000011", "BBB"))}
by_name = {g.name: g for g in transcripts.assign(genes, table, mane)}

g = by_name["OLDNAME"]
check("an old SpliceAI symbol finds its MANE transcript by locus (shared exons)",
      g.canonical.accession == "NM_000001.1" and g.rule == transcripts.RULE_MANE, f"{g.canonical} {g.rule}")
check("the gene's other NM_ on its strand join the mask; the other strand doesn't",
      [t.accession for t in g.transcripts] == ["NM_000001.1", "NM_000002.1"], str([t.accession for t in g.transcripts]))
check("an NM_ sharing exons with another gene's MANE transcript is left out, with the reason",
      g.excluded and g.excluded[0][0] == "NM_000003.1" and "NEIGHBOR" in g.excluded[0][1], str(g.excluded))
check("mask = SpliceAI's interval joined with those transcripts", g.mask == (901, 9000), str(g.mask))
g = by_name["RT"]
check("a locus shared with two MANE genes takes the one with most shared exon bases, and lists the other",
      g.canonical.accession == "NM_000010.1" and g.other_mane == ["BBB"], f"{g.canonical} {g.other_mane}")
g = by_name["NOMANE"]
check("no MANE: the own-symbol NM_ with the longest CDS, not the host gene's longer one",
      g.canonical.accession == "NM_000020.1" and g.rule == transcripts.RULE_LONGEST_CDS, f"{g.canonical} {g.rule}")
g = by_name["NESTED"]
check("a gene inside another's intron doesn't get the host's transcript; mask stays SpliceAI's",
      g.canonical is None and g.rule == transcripts.RULE_NONE and g.mask == (5001, 6000), f"{g.canonical} {g.mask}")

# --- the real tables ---------------------------------------------------------------------
real = {g.name: g for g in transcripts.load_genes("hg38", REFSEQ, MANE)}
check("hg38: BRCA2 and CHEK2 get their MANE Select transcripts",
      real["BRCA2"].canonical.accession == "NM_000059.4" and real["CHEK2"].canonical.accession == "NM_007194.4")
check("hg38: SpliceAI's FAM132A is matched to C1QTNF12 by locus",
      real["FAM132A"].canonical.symbol == "C1QTNF12", str(real["FAM132A"].canonical))
real19 = {g.name: g for g in transcripts.load_genes("hg19", REFSEQ, MANE)}
check("hg19: ANKRD26 gets the MANE Select NM_ as aligned to GRCh37",
      real19["ANKRD26"].canonical.accession == "NM_014915.3", str(real19["ANKRD26"].canonical))

# --- what a gene stores ---------------------------------------------------------------
def plan(start, end, mask, chrom_len=1_000_000):
    gene = Gene("G", "1", "+", start, end, ((start, end),))
    if mask != (start, end):
        gene.transcripts = [tx("NM_9.1", "G", [mask])]
    return refdb.plan_gene(gene, chrom_len)


p = plan(100_000, 200_000, (100_000, 200_000))
check("mask unmoved: stored = mask +/- 50, all comparable with the VCF tab",
      (p["start"], p["end"], p["vcf_comparable"], p["chrom_end_padded"]) == (99_950, 200_050, [99_950, 200_050], []), str(p))
p = plan(100_000, 200_000, (90_000, 200_000))
check("mask moved at the start: comparable from 5,000 bases inside SpliceAI's start",
      (p["start"], p["vcf_comparable"]) == (89_950, [105_000, 200_050]), str(p))
p = plan(100_000, 200_000, (100_000, 260_000))
check("mask moved at the end: comparable to 5,000 bases inside SpliceAI's end",
      (p["end"], p["vcf_comparable"]) == (260_050, [99_950, 195_000]), str(p))
p = plan(2_000, 30_000, (2_000, 30_000), chrom_len=36_000)
check("near the chromosome start: positions whose context runs off it are marked, and not comparable",
      (p["start"], p["end"], p["chrom_end_padded"], p["vcf_comparable"]) == (1_950, 30_050, [[1_950, 5_000]], [5_001, 30_050]), str(p))
p = plan(20, 3_000, (20, 3_000), chrom_len=1_000_000)
check("...and the stored region never starts before position 1", (p["start"], p["chrom_end_padded"]) == (1, [[1, 3_050]]), str(p))
p = plan(10_000, 34_000, (10_000, 34_000), chrom_len=36_000)
check("...and at the far end (the last 5,000 bases of the chromosome)",
      (p["end"], p["chrom_end_padded"], p["vcf_comparable"]) == (34_050, [[31_001, 34_050]], [9_950, 31_000]), str(p))

# --- storage precision -----------------------------------------------------------------
probs = np.random.default_rng(1).random((100_000, 2)).astype(np.float32) ** 4
probs[:10] = [[0.0, 1.0], [0.99976, 0.9998], [0.0004, 0.0006], [0.05, 0.5], [1e-6, 0.123456]] * 2
for name, bound in (("float32", 0.0), ("float16", 2 ** -12), ("q001", 0.0005 + 1e-6)):
    back = refdb.decode(refdb.encode(probs, name), name)
    err = float(np.abs(back - probs).max())
    check(f"{name}: stored and read back within {bound:g}", err <= bound, f"max error {err:.3g}")

# --- a real build: killed mid-way, then resumed --------------------------------------
chr21 = sorted((g for g in real.values() if g.chrom == "21" and 4_000 <= g.mask[1] - g.mask[0] <= 25_000),
               key=lambda g: g.start)
names = [g.name for g in chr21[:4]]
extra = chr21[4].name
out = os.path.join(tmp, "db")
common = ["--build", "hg38", "--fasta", str(helpers.HG38), "--out", out,
          "--refseq-dir", str(REFSEQ), "--mane-dir", str(helpers.MANE_DIR)]
env = dict(os.environ, PYTHONPATH=str(helpers.SPLICEAI_PKG), TF_CPP_MIN_LOG_LEVEL="3")
env["SPLICEAI_CLI_LOG"] = os.path.join(tmp, "cli.log")

# The first run goes through the exe's entry point (run_app.py --build-refdb).
proc = subprocess.Popen([sys.executable, "-W", "ignore", str(helpers.REPO / "run_app.py"), "--build-refdb",
                         *common, "--genes", ",".join(names)], env=env, cwd=helpers.REPO)
deadline = time.time() + 300
while time.time() < deadline and proc.poll() is None and not glob.glob(os.path.join(out, "_work", "chr21", "*.npz")):
    time.sleep(0.02)
proc.kill()
proc.wait()
done_before = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(out, "_work", "chr21", "*.npz")))
check("killed after the first finished gene: some genes done, chromosome not packed yet",
      1 <= len(done_before) < len(names) and not os.path.exists(os.path.join(out, "chr21.h5")), str(done_before))

run = subprocess.run([sys.executable, "-W", "ignore", "-m", "spliceai_pipeline.refdb_build", *common,
                      "--genes", ",".join(names)], env=env, cwd=helpers.REPO, capture_output=True, text=True)
log = run.stdout
check("resumed run finishes", run.returncode == 0, (run.stderr or log)[-300:])
check("...continuing after the genes already finished, not redoing them",
      f"{len(done_before)} already finished" in log
      and not any(f" {n}: " in log for n in done_before), [ln for ln in log.splitlines() if "chr21" in ln][:3])
db = refdb.Database(out)
manifest = json.load(open(os.path.join(out, "manifest.json")))
check("manifest lists the four genes, nothing left unpacked",
      sorted(manifest["genes"]) == sorted(names) and not manifest["unpacked_genes"] and manifest["scope"].startswith("partial"))
check("manifest: SHA-256 and size of chr21.h5 are right",
      manifest["files"] == [{"file": "chr21.h5", "size": os.path.getsize(os.path.join(out, "chr21.h5")),
                             "sha256": refdb.sha256_of(os.path.join(out, "chr21.h5"))}])
check("manifest records the five model hashes and the annotation files",
      manifest["model_sha256"] == predict.MODEL_SHA256 and len(manifest["annotation"]) == 3)
check("manifest states the mask rule and how the VCF-comparable stretch is drawn",
      manifest["mask_rule"] == refdb.MASK_RULE and manifest["vcf_comparable_rule"] == refdb.VCF_COMPARABLE_RULE)
check("the build kept the PC awake (it says so in its log)", "Keeping the PC awake" in log)

from spliceai_pipeline.refdb_build import ES_CONTINUOUS, ES_SYSTEM_REQUIRED, keep_awake  # noqa: E402
keep_awake(True)
held = keep_awake(False)
released = keep_awake(False)
check("keep-awake asks Windows for ES_CONTINUOUS | ES_SYSTEM_REQUIRED, and lets go again",
      held == ES_CONTINUOUS | ES_SYSTEM_REQUIRED and released == ES_CONTINUOUS, f"{held:#x}, then {released:#x}")

hg38 = predict.Predictor(str(helpers.HG38), "hg38")
ref = hg38.annotator.ref_fasta["chr21"]
worst = 0.0
for name in names:
    stored, plan_ = db.read_gene("21", name)
    expected = refdb.plan_gene(real[name], len(ref))
    fresh = hg38.predict(lambda a, b: ref[a - 1:b].seq, len(ref), plan_["strand"], *plan_["mask"],
                         plan_["start"], plan_["end"])
    worst = max(worst, float(np.abs(stored - fresh).max()))
    check(f"{name}: stored plan is the gene's plan", plan_ == expected)
check("stored values (resumed and fresh genes alike) equal a fresh prediction within float16",
      worst <= 2 ** -12 + 1e-6, f"max |diff| {worst:.3g}")

before = db.read_gene("21", names[0])[0]
run = subprocess.run([sys.executable, "-W", "ignore", "-m", "spliceai_pipeline.refdb_build", *common,
                      "--genes", extra], env=env, cwd=helpers.REPO, capture_output=True, text=True)
manifest = json.load(open(os.path.join(out, "manifest.json")))
check("a later gene on a packed chromosome is added to it; the others are kept as they were",
      run.returncode == 0 and sorted(manifest["genes"]) == sorted(names + [extra])
      and np.array_equal(before, db.read_gene("21", names[0])[0]), (run.stderr or run.stdout)[-300:])

run = subprocess.run([sys.executable, "-W", "ignore", "-m", "spliceai_pipeline.refdb_build", *common,
                      "--genes", extra, "--precision", "q001"], env=env, cwd=helpers.REPO, capture_output=True, text=True)
check("continuing with other settings (precision) is refused, and says why",
      run.returncode != 0 and "other settings (precision)" in run.stderr, run.stderr[-300:])

# --- the validation's comparison (spec 5.6) --------------------------------------------
from spliceai_pipeline import refdb_validate as rv  # noqa: E402

flat = [np.zeros(101)] * 4
problems, flips = rv.compare([1e-9, 0, 0, 0], [0, 0, 0, 0], flat, [-2e-11, 0, 0, 0], [0, 0, 0, 0])
check("validation: 0.00 against -0.00 is neither a problem nor a rounding flip", not problems and not flips)
problems, flips = rv.compare([0.995332, 0, 0, 0], [2, 0, 0, 0], flat, [0.99495, 0, 0, 0], [2, 0, 0, 0])
check("validation: 0.995332 against 0.994950 (0.00038 apart) is out of tolerance, and listed as a 1.00 / 0.99 flip",
      problems and flips, f"{problems} {flips}")
tie = [np.r_[np.zeros(99), 0.5, 0.5 - 1e-6]] + [np.zeros(101)] * 3
problems, _ = rv.compare([0.5, 0, 0, 0], [49, 0, 0, 0], tie, [0.5, 0, 0, 0], [50, 0, 0, 0])
check("validation: a DP difference between two near-equal positions isn't a problem", not problems, str(problems))
clear = [np.r_[np.zeros(100), 0.5]] + [np.zeros(101)] * 3
problems, _ = rv.compare([0.5, 0, 0, 0], [50, 0, 0, 0], clear, [0.5, 0, 0, 0], [49, 0, 0, 0])
check("validation: a DP difference at a clear maximum with DS >= 0.05 is", bool(problems), str(problems))

manifest = json.load(open(os.path.join(out, "manifest.json")))
done = 0
for name, g in manifest["genes"].items():
    if not g["vcf_comparable"] or g["vcf_comparable"][1] - g["vcf_comparable"][0] < 400:
        continue
    lo, hi = g["vcf_comparable"]
    for pos, kind in ((lo + 120, "SNV"), (lo + 200, "1-base insertion"), (hi - 120, "2-base deletion")):
        base = ref[pos - 1:pos + 1].seq.upper()
        ref_allele, alt = {"SNV": (base[0], "A" if base[0] != "A" else "C"),
                           "1-base insertion": (base[0], base[0] + "T"), "2-base deletion": (base[:2], base[0])}[kind]
        row = rv.check_variant(hg38.annotator, db, manifest["genes"], name, "chr21", pos, ref_allele, alt)
        check(f"validation, {name} {kind}: database reference gives get_delta_scores' delta scores "
              "(raw and masked; recomputation checked against SpliceAI's own output)",
              row is not None and not row["raw"]["problems"] and not row["masked"]["problems"],
              f"max |dDS| {max(row['raw']['max_abs_diff'], row['masked']['max_abs_diff']):.2g}" if row else "not checked")
        done += 1
    unmarked = dict(manifest["genes"], **{name: dict(g, vcf_comparable=None)})
    base = ref[lo + 119:lo + 120].seq.upper()
    args = (hg38.annotator, db, unmarked, name, "chr21", lo + 120, base, "A" if base != "A" else "C")
    row = rv.check_variant(*args, strict=False)
    check("validation: a variant outside the VCF-comparable stretch is skipped, and compared anyway with strict=False",
          rv.check_variant(*args) is None and row is not None and not row["masked"]["problems"])
    break
if not done:
    helpers.skip("validation on a built gene", "no built gene has a VCF-comparable stretch of 400 bases")

# --- model files are checked ---------------------------------------------------------
fake = os.path.join(tmp, "fake-spliceai")
os.makedirs(os.path.join(fake, "models"))
for name in predict.MODEL_SHA256:
    with open(os.path.join(fake, "models", name), "wb") as fh:
        fh.write(b"not a model")
real_dir, predict.spliceai_package_dir = predict.spliceai_package_dir, lambda: fake
try:
    predict.verify_models()
    check("model files that aren't the pinned ones are refused", False, "no error")
except predict.ModelMismatchError as exc:
    check("model files that aren't the pinned ones are refused", "spliceai1.h5" in str(exc), str(exc)[:100])
finally:
    predict.spliceai_package_dir = real_dir

helpers.finish()
