"""The reference score database in the app (spec 5.7, 5.8, Phase 3): a build
steered the way the window steers it (pause, resume, stop, a gene asked for
first, progress), then the reader with its refusals, the live fallback and
its cache, the verification against a validation set, and import.

    .venv\\Scripts\\python tests\\test_refdb_reader.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

import helpers
from helpers import check

helpers.watchdog(1200)
helpers.setup()

import h5py  # noqa: E402
import numpy as np  # noqa: E402

from spliceai_pipeline import predict, refdb, transcripts  # noqa: E402

REFSEQ = helpers.THIRDPARTY / "refseq"
MANE = helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz"
MODELS = predict.verify_models()

genes = transcripts.load_genes("hg38", REFSEQ, MANE)
small = lambda chrom: sorted((g for g in genes if g.chrom == chrom and 4_000 <= g.mask[1] - g.mask[0] <= 25_000),
                             key=lambda g: g.start)
chr21, chr22 = small("21"), small("22")
names = [g.name for g in chr21[:6]]
asked = chr22[0].name           # "Build this gene now", on another chromosome
live_only = chr21[8].name       # never built: computed live by the fallback

tmp = tempfile.mkdtemp(prefix="spliceai_refdb_reader_")
out = os.path.join(tmp, "hg38")
genes_file = os.path.join(tmp, "genes.txt")
with open(genes_file, "w", encoding="utf-8") as fh:
    fh.write("\n".join(names) + "\nNOT_A_GENE_XYZ\n")
cmd = [sys.executable, "-W", "ignore", "-m", "spliceai_pipeline.refdb_build", "--build", "hg38",
       "--fasta", str(helpers.HG38), "--out", out, "--refseq-dir", str(REFSEQ), "--mane-dir", str(helpers.MANE_DIR),
       "--genes-file", genes_file, "--ignore-unknown", "--chunk", "4000", "--no-keep-awake"]
env = dict(os.environ, PYTHONPATH=str(helpers.SPLICEAI_PKG), TF_CPP_MIN_LOG_LEVEL="3")


def progress():
    try:
        return refdb.read_json(os.path.join(out, "progress.json"))
    except (OSError, ValueError):
        return {}


def work():
    return sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(out, "_work", "*", "*.npz")))


def wait(what, until, seconds=240):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if until():
            return True
        time.sleep(0.1)
    check(f"waiting for {what}", False, f"not within {seconds} s")
    return False


def control(state):
    with open(os.path.join(out, "control.txt"), "w", encoding="utf-8") as fh:
        fh.write(state)


# --- a build steered through its files ------------------------------------------------
proc = subprocess.Popen(cmd, env=env, cwd=helpers.REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
wait("the first finished gene", lambda: work())
with open(os.path.join(out, "priority.txt"), "w", encoding="utf-8") as fh:
    fh.write(asked + "\n")
control("pause")
if wait("the build to pause", lambda: progress().get("state") == "paused", 60):
    held = work()
    time.sleep(3)
    p = progress()
    check("paused: no gene finishes while paused, and progress.json says so",
          work() == held and p["state"] == "paused", f"{held} -> {work()}")
    check("progress.json counts bases done of the whole scope (4 kb-25 kb genes, plus 100 flanking bases)",
          0 < p["bases_done"] < p["bases_total"] and p["bases_total"] == sum(
              g.mask[1] - g.mask[0] + 101 for g in chr21[:6]), f"{p['bases_done']:,} of {p['bases_total']:,}")
    reader = refdb.Reader(out)
    loose = held[0]
    plan = reader.plan(loose)
    check("a finished gene can be read before its chromosome is packed",
          plan is not None and reader.values(loose).shape == (plan["end"] - plan["start"] + 1, 2))
    loose_values = reader.values(loose)
control("run")
wait("the gene asked for", lambda: os.path.exists(os.path.join(out, "_work", f"chr{chr22[0].chrom}", f"{asked}.npz")))
control("stop")
try:
    proc.wait(120)
except subprocess.TimeoutExpired:
    proc.kill()
log = proc.stdout.read()
check("stopped on request: exits normally, progress.json says stopped, finished genes kept",
      proc.returncode == 0 and progress().get("state") == "stopped" and asked in work(), f"exit {proc.returncode}, {work()}")
check("the gene asked for was built first, though it's on another chromosome and outside the list",
      f"chr22 on request {asked}:" in log, [ln for ln in log.splitlines() if "request" in ln][:2])
check("an unknown gene in the list is left out with a note, not a stop",
      "left out: NOT_A_GENE_XYZ" in log)

os.remove(os.path.join(out, "control.txt"))
run = subprocess.run(cmd, env=env, cwd=helpers.REPO, capture_output=True, text=True)
manifest = refdb.read_json(os.path.join(out, "manifest.json"))
check("resumed: finishes, and packs the chromosome of the gene asked for as well",
      run.returncode == 0 and sorted(manifest["genes"]) == sorted(names + [asked])
      and progress().get("state") == "finished" and not work(), (run.stderr or run.stdout)[-300:])
reader = refdb.Reader(out)
check("a gene reads the same before and after packing (rounded alike)",
      np.array_equal(loose_values, reader.values(loose)))

# --- reading -----------------------------------------------------------------------
check("check(): the database is accepted for hg38 with the installed models",
      refdb.check(out, "hg38", MODELS)["build"] == "hg38")
check("checked.json remembers the files found intact", sorted(refdb.read_json(os.path.join(out, "checked.json")))
      == ["chr21.h5", "chr22.h5"])
db = refdb.Database(out)
name = names[0]
plan = reader.plan(name)
full = reader.values(name)
check("Reader.values gives what the builder stored", np.array_equal(full, db.read_gene("21", name)[0]))
a, b = plan["start"] + 100, plan["start"] + 199
check("Reader.values(start, end): that stretch, both ends included", np.array_equal(reader.values(name, a, b), full[100:200]))
try:
    reader.values(name, plan["start"] - 1, plan["start"] + 10)
    check("a stretch outside the gene's stored region is refused", False)
except ValueError as exc:
    check("a stretch outside the gene's stored region is refused", "is stored at" in str(exc), str(exc))
check("genes_at finds the gene by position (bare or chr-prefixed), and filters by strand",
      name in reader.genes_at("chr21", a, b) and name in reader.genes_at("21", a, a, plan["strand"])
      and name not in reader.genes_at("21", a, b, "+" if plan["strand"] == "-" else "-"))
check("a database is None where there is none", refdb.check(os.path.join(tmp, "nothing"), "hg38", MODELS) is None)


def refused(why, folder, build="hg38", models=MODELS, expect=""):
    try:
        refdb.check(folder, build, models)
        check(f"refused: {why}", False, "accepted")
    except refdb.DatabaseRefused as exc:
        check(f"refused: {why}, and says so", expect in str(exc), str(exc)[:160])


def copy_db(tag, change=None):
    folder = os.path.join(tmp, tag)
    shutil.copytree(out, folder)
    if change:
        change(folder)
    return folder


def edit_manifest(folder, **fields):
    path = os.path.join(folder, "manifest.json")
    refdb.write_json(path, dict(refdb.read_json(path), **fields))


refused("another genome build", out, build="hg19", expect="is for hg38, not hg19")
refused("other model files", out, models=dict(MODELS, **{"spliceai3.h5": "0" * 64}),
        expect="other SpliceAI model files (spliceai3.h5")
refused("another format", copy_db("format", lambda f: edit_manifest(f, format=99)), expect="format 99")
refused("a failed verification", copy_db("failed", lambda f: edit_manifest(
    f, local_verification={"passed": False, "date": "2026-10-27", "problems": ["X's values differ"]})),
    expect="failed its verification on 2026-10-27 (X's values differ)")


def flip_byte(folder):
    path = os.path.join(folder, "chr21.h5")
    with open(path, "r+b") as fh:
        fh.seek(os.path.getsize(path) // 2)
        byte = fh.read(1)
        fh.seek(-1, 1)
        fh.write(bytes([byte[0] ^ 0xFF]))


refused("a damaged file (same size, other bytes), though checked.json had it as intact",
        copy_db("damaged", flip_byte), expect="damaged: chr21.h5")
refused("a missing file", copy_db("missing", lambda f: os.remove(os.path.join(f, "chr22.h5"))),
        expect="incomplete: chr22.h5 is missing")

# --- the live fallback and its cache --------------------------------------------------
cache = os.path.join(tmp, "cache")
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=reader, cache_dir=cache,
                               refseq_dir=str(REFSEQ), mane_dir=str(helpers.MANE_DIR))
values, plan_, source = scores.gene(name)
check("fallback: a gene in the database comes from it", source == "database" and np.array_equal(values, full))
t0 = time.time()
values, plan_, source = scores.gene(live_only)
first = time.time() - t0
ref = scores._predictor.annotator.ref_fasta["chr21"]
fresh = scores._predictor.predict(lambda a, b: ref[a - 1:b].seq, len(ref), plan_["strand"], *plan_["mask"],
                                  plan_["start"], plan_["end"])
check("fallback: a gene the database lacks is computed live, the database's way (same plan, float16)",
      source == "computed live" and plan_ == refdb.plan_gene(next(g for g in genes if g.name == live_only), len(ref))
      and np.array_equal(values, refdb.decode(refdb.encode(fresh, "float16"), "float16")))
again = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=reader, cache_dir=cache,
                              refseq_dir=str(REFSEQ), mane_dir=str(helpers.MANE_DIR))
t0 = time.time()
values2, _, source2 = again.gene(live_only)
check("...and cached: the next time it is read back, without loading the models",
      source2 == "computed live" and np.array_equal(values, values2) and again._predictor is None,
      f"{first:.1f} s computed, {time.time() - t0:.2f} s from the cache")
path = os.path.join(cache, f"{live_only}.npz")
with np.load(path) as z:
    stale = dict(z)
stale["settings"] = json.dumps(dict(json.loads(str(stale["settings"])), precision="float32"))
np.savez_compressed(path, **stale)
again._settings = None
_, _, source3 = again.gene(live_only)
check("a cached gene made with other settings is computed again", source3 == "computed live" and again._predictor is not None)
check("the cache's size can be shown", refdb.folder_size(cache) > 0, f"{refdb.folder_size(cache):,} bytes")
try:
    scores.gene("NOT_A_GENE_XYZ")
    check("fallback: an unknown gene is a KeyError", False)
except KeyError as exc:
    check("fallback: an unknown gene is a KeyError, saying so", "not in SpliceAI's hg38 gene table" in str(exc))

# --- verification against a validation set (spec 5.8) -----------------------------------
def validation_set(path, chosen, change=None):
    with h5py.File(path, "w") as vs:
        vs.attrs.update(build="hg38", precision="float16", model_sha256=json.dumps(MODELS))
        for n in chosen:
            v, plan, _ = scores.gene(n)     # live_only: from the fallback's cache
            if change:
                v = change(v)
            group = vs.create_group(n)
            for i, track in enumerate(("acceptor", "donor")):
                group.create_dataset(track, data=refdb.encode(v[:, i], "float16"))
            group.attrs["plan"] = json.dumps(plan)
    return path


record = refdb.verify_against_set(out, validation_set(os.path.join(tmp, "vs-same.h5"), names[:3] + [live_only]))
check("verification: genes in both match (a gene the database lacks is passed over)",
      record["passed"] is True and record["genes"] == 3 and record["max_abs_diff"] == 0.0, str(record)[:200])
nudge = lambda v: v + np.float32(2 ** -10)
record = refdb.verify_against_set(out, validation_set(os.path.join(tmp, "vs-off.h5"), names[:2], change=nudge))
check("verification: values further apart than two float16 rounding steps fail, naming each gene",
      record["passed"] is False and all(any(n in p for p in record["problems"]) for n in names[:2]),
      str(record["problems"])[:200])
nudge = lambda v: v + np.float32(2 ** -12)
record = refdb.verify_against_set(out, validation_set(os.path.join(tmp, "vs-near.h5"), names[:2], change=nudge))
check("verification: values one float16 rounding apart (another PC's float32) pass", record["passed"] is True,
      f"max |diff| {record['max_abs_diff']:.6f}")
record = refdb.verify_against_set(out, validation_set(os.path.join(tmp, "vs-none.h5"), [live_only]))
check("verification: no gene in both is 'not verified', not a pass", record["passed"] is None and record["genes"] == 0)

refdb.record_verification(out, refdb.verify_against_set(out, os.path.join(tmp, "vs-same.h5")))
subprocess.run(cmd, env=env, cwd=helpers.REPO, capture_output=True, text=True)
check("the verification stays in the manifest when a build changes nothing",
      refdb.read_json(os.path.join(out, "manifest.json")).get("local_verification", {}).get("passed") is True)

# The databases this PC built in Phase 2, if they are there (read only).
real_set = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb")
vs_path = os.path.join(real_set, "validation", "hg38", "validation-set-hg38.h5")
if os.path.exists(vs_path) and os.path.exists(os.path.join(real_set, "hg38", "manifest.json")):
    record = refdb.verify_against_set(os.path.join(real_set, "hg38"), vs_path)
    check("the real hg38 database passes against its validation set",
          record["passed"] is True and record["genes"] >= 50, f"{record['genes']} genes, max |diff| {record['max_abs_diff']:.3g}")
else:
    helpers.skip("the real hg38 database against its validation set", "not built on this PC")

# --- import --------------------------------------------------------------------------------
root = os.path.join(tmp, "imported")
dest = refdb.import_database(out, MODELS, root=root)
check("import: copied to <root>/hg38 and accepted there",
      dest == os.path.join(root, "hg38") and refdb.check(dest, "hg38", MODELS)["build"] == "hg38"
      and not os.path.exists(dest + ".importing") and not os.path.exists(os.path.join(dest, "progress.json")))
seen = []
dest = refdb.import_database(out, MODELS, root=root, on_progress=lambda d, t: seen.append((d, t)))
check("import over an existing database replaces it, reporting progress",
      refdb.check(dest, "hg38", MODELS) and seen and seen[-1][0] == seen[-1][1] and not os.path.exists(dest + ".replaced"))
for why, folder, expect in (
        ("no manifest", os.path.join(tmp, "nothing-here"), "no manifest.json"),
        ("a damaged copy (nothing is replaced)", os.path.join(tmp, "damaged"), "damaged: chr21.h5"),
        ("the database already in use", dest, "already the hg38 database")):
    os.makedirs(folder, exist_ok=True)
    before = refdb.read_json(os.path.join(dest, "manifest.json"))
    try:
        refdb.import_database(folder, MODELS, root=root)
        check(f"import refused: {why}", False, "imported")
    except refdb.DatabaseRefused as exc:
        check(f"import refused: {why}", expect in str(exc) and refdb.read_json(os.path.join(dest, "manifest.json")) == before,
              str(exc)[:120])

# --- benchmark ------------------------------------------------------------------------------
bench = os.path.join(tmp, "bench.json")
run = subprocess.run([*cmd[:cmd.index("--out")], "--benchmark", bench, "--refseq-dir", str(REFSEQ), "--mane-dir",
                      str(helpers.MANE_DIR)], env=env, cwd=helpers.REPO, capture_output=True, text=True)
result = refdb.read_json(bench) if os.path.exists(bench) else {}
check("benchmark: times about 600,000 bases of chromosome 21 and writes the speed, building nothing",
      run.returncode == 0 and result.get("bases", 0) >= 600_000 and result.get("bases_per_second", 0) > 0,
      f"{result.get('bases_per_second', 0):,} bases/s in {result.get('seconds')} s" if result else run.stderr[-300:])

helpers.finish()
