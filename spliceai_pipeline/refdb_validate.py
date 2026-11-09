"""Validation of a reference score database (spec 5.6), and its validation set
(spec 5.8).

    python -m spliceai_pipeline.refdb_validate --db refdb\\hg38 --fasta hg38.fa --out validation\\hg38 ^
        --refseq-dir installer\\thirdparty\\refseq --mane-dir installer\\thirdparty\\mane

1. Coverage: every gene of the annotation stored or skipped with a reason;
   each gene stored exactly as plan_gene() says (mask +/- 50 bases, covering
   all its transcripts), its arrays that long; every file's size and SHA-256
   as in the manifest; nothing left unpacked.
2. Stored values against a fresh prediction (predict.py, another chunk size)
   at every position, for >= 200 genes stratified by size and strand plus
   genes overlapping another; the largest absolute difference is reported.
3. Delta scores with the stored reference against SpliceAI's own
   get_delta_scores (raw, and masked from the same predictions), for SNVs,
   1-base insertions and 1-4 base deletions inside each gene's VCF-comparable
   stretch, near gene ends and at the public regression variants: unrounded DS
   within TOLERANCE, DP identical where DS >= 0.05 and the top position beats
   the runner-up by more than 2 x TOLERANCE, every 2-decimal flip listed.
Then about 50 of the checked genes (stratified, at most 200 kb each) are
written as the validation set: their fresh predictions, float16.

delta_scores() is also what the patient regression checks use
(tests/test_pipeline.py): the unrounded numbers get_delta_scores computes,
checked against its own printed output every time.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

# Half of float16's largest step below 1 (the stored precision), plus 1e-5.
TOLERANCE = 2 ** -12 + 1e-5
COV = 101   # get_delta_scores' window with D = 50
HALF = COV // 2
CHECK_CHUNK = 37_777    # not the builder's chunk size, on purpose


class _Recorder:
    """Stands in for one model inside get_delta_scores and keeps its outputs."""

    def __init__(self, model):
        self.model, self.outputs = model, []

    def predict(self, x, *args, **kwargs):
        y = self.model.predict(x, *args, **kwargs)
        self.outputs.append(y)
        return y


def _adjust_alt(y_alt, ref_len, alt_len):
    """get_delta_scores' indel handling of the alternative prediction."""
    if ref_len > 1 and alt_len == 1:
        y_alt = np.concatenate([y_alt[:, :HALF + alt_len], np.zeros((1, ref_len - alt_len, 3)),
                                y_alt[:, HALF + alt_len:]], axis=1)
    elif ref_len == 1 and alt_len > 1:
        y_alt = np.concatenate([y_alt[:, :HALF], np.max(y_alt[:, HALF:HALF + alt_len], axis=1)[:, None, :],
                                y_alt[:, HALF + alt_len:]], axis=1)
    return y_alt


def scores_from(y_ref, y_alt, dist_exon, masked):
    """(DS [AG, AL, DG, DL], DP [same], the four delta series) as
    get_delta_scores computes them, unrounded. y_ref, y_alt: (1, 101, 3)."""
    y = np.concatenate([y_ref, y_alt])
    series = [y[1, :, 1] - y[0, :, 1], y[0, :, 1] - y[1, :, 1], y[1, :, 2] - y[0, :, 2], y[0, :, 2] - y[1, :, 2]]
    idx = [int(s.argmax()) for s in series]
    ds = []
    for j, (s, i) in enumerate(zip(series, idx)):
        # Gains are masked at an annotated exon boundary, losses away from one.
        at_boundary = i - HALF == dist_exon
        m = bool(masked) and (at_boundary if j in (0, 2) else not at_boundary)
        ds.append(float(s[i] * (1 - m)))
    return ds, [i - HALF for i in idx], series


def delta_scores(annotator, chrom, pos, ref, alt, masked):
    """Runs SpliceAI's get_delta_scores for one variant with the models'
    outputs recorded. Returns (its annotation strings, {gene: dict}); each dict
    holds the unrounded "ds" / "dp" / "series", the reference prediction
    "y_ref" (1, 101, 3, genomic order), the adjusted alternative "y_alt" and
    "dist_exon", so scores_from() can recompute them with another reference.
    Raises AssertionError unless the recomputation, rounded, is exactly
    get_delta_scores' own output -- so it can't drift from SpliceAI's."""
    from .score import _Record, get_delta_scores

    models, annotator.models = annotator.models, [_Recorder(m) for m in annotator.models]
    try:
        out = get_delta_scores(_Record(chrom, pos, ref, alt), annotator, 50, 1 if masked else 0)
        recorders = annotator.models
    finally:
        annotator.models = models
    result = {}
    if not out or (len(ref) > 1 and len(alt) > 1):
        return out, result
    genes, strands, idxs = annotator.get_name_and_strand(chrom, pos)
    for k, (gene, strand, i) in enumerate(zip(genes, strands, idxs)):
        y_ref = np.mean([r.outputs[2 * k] for r in recorders], axis=0)
        y_alt = np.mean([r.outputs[2 * k + 1] for r in recorders], axis=0)
        if strand == "-":
            y_ref, y_alt = y_ref[:, ::-1], y_alt[:, ::-1]
        y_alt = _adjust_alt(y_alt, len(ref), len(alt))
        dist_exon = annotator.get_pos_data(i, pos)[2]
        ds, dp, series = scores_from(y_ref, y_alt, dist_exon, masked)
        printed = out[k].split("|")
        mine = [f"{v:.2f}" for v in ds] + [str(p) for p in dp]
        assert printed[1] == gene and printed[2:] == mine, f"recomputed {mine} != get_delta_scores {printed}"
        result[str(gene)] = {"ds": ds, "dp": dp, "series": series, "y_ref": y_ref, "y_alt": y_alt,
                             "dist_exon": dist_exon}
    return out, result


def compare(ds_a, dp_a, series_a, ds_b, dp_b, tol=TOLERANCE):
    """Spec 5.6, a being the reference: ([problems], [2-decimal flips]).
    Numbers, not strings: -0.00 and 0.00 are the same score."""
    problems, flips = [], []
    for j, name in enumerate(("AG", "AL", "DG", "DL")):
        if abs(ds_a[j] - ds_b[j]) > tol:
            problems.append(f"DS_{name} {ds_a[j]:.6f} vs {ds_b[j]:.6f}")
        if round(ds_a[j], 2) + 0.0 != round(ds_b[j], 2) + 0.0:
            flips.append(f"DS_{name} {ds_a[j]:.6f} -> {ds_a[j]:.2f}, {ds_b[j]:.6f} -> {ds_b[j]:.2f}")
        if ds_a[j] >= 0.05 and dp_a[j] != dp_b[j]:
            top = np.sort(series_a[j]) if series_a is not None else None
            if top is None or top[-1] - top[-2] > 2 * tol:
                problems.append(f"DP_{name} {dp_a[j]} vs {dp_b[j]}")
    return problems, flips


PUBLIC_REGRESSION = {
    "hg19": [("ANKRD26", "chr10", 27326999, "T", "C"), ("CHEK2", "chr22", 29121230, "C", "T"),
             ("MYBPC3", "chr11", 47364709, "C", "T"), ("BRCA2", "chr13", 32945092, "G", "A")],
    "hg38": [("ABCA4", "chr1", 94011395, "A", "G"), ("APC", "chr5", 112827247, "G", "A"),
             ("CFTR", "chr7", 117639961, "C", "T"), ("ANKRD26", "chr10", 27038070, "T", "C"),
             ("MYBPC3", "chr11", 47343158, "C", "T"), ("ATM", "chr11", 108309110, "A", "G"),
             ("BRCA2", "chr13", 32370955, "G", "A"), ("CHEK2", "chr22", 28725242, "C", "T")],
}


def _overlapping(genes):
    """Names of genes whose SpliceAI interval overlaps another gene's."""
    spans = sorted((g["chrom"], g["spliceai"][0], g["spliceai"][1], n) for n, g in genes.items())
    found, reach = set(), {}
    for chrom, s, e, n in spans:
        end, other = reach.get(chrom, (0, None))
        if s <= end:
            found |= {n, other}
        if e > end:
            reach[chrom] = (e, n)
    return found


def _sample_genes(genes, rng, per_stratum=25, overlapping=25):
    """[(name, stratum)]: per_stratum genes from each stored-size quartile and
    strand, plus `overlapping` genes that overlap another gene."""
    by_len = sorted(genes, key=lambda n: genes[n]["end"] - genes[n]["start"])
    picked = []
    for q in range(4):
        quartile = by_len[q * len(by_len) // 4:(q + 1) * len(by_len) // 4]
        for strand in "+-":
            pool = [n for n in quartile if genes[n]["strand"] == strand]
            picked += [(str(n), f"size quartile {q + 1}, {strand} strand")
                       for n in rng.choice(pool, min(per_stratum, len(pool)), replace=False)]
    extra = sorted(_overlapping(genes) - {n for n, _ in picked})
    picked += [(str(n), "overlaps another gene") for n in rng.choice(extra, min(overlapping, len(extra)), replace=False)]
    return picked


def _variants(genes, rng, fasta, n_random, build):
    """[(gene, chrom, pos, ref, alt, why)]; the window of each lies inside the
    gene's VCF-comparable stretch."""
    out = []

    def add(name, pos, kind, why):
        g = genes[name]
        if not g["vcf_comparable"]:
            return
        lo, hi = g["vcf_comparable"]
        ref_len = kind if kind >= 2 else 1          # kind: 0 SNV, 1 insertion, 2-5 deletion of kind-1 bases
        if not (lo + HALF <= pos and pos + ref_len - 1 + HALF <= hi and g["spliceai"][0] <= pos <= g["spliceai"][1]):
            return
        contig = "chr" + g["chrom"]
        refb = fasta[contig][pos - 1:pos - 1 + ref_len].seq.upper()
        if set(refb) - set("ACGT"):
            return
        alt = (next(b for b in "ACGT" if b != refb) if kind == 0 else
               refb + "ACGT"[int(rng.integers(4))] if kind == 1 else refb[0])
        out.append((name, contig, pos, refb, alt, why))

    names = sorted(n for n, g in genes.items() if g["vcf_comparable"])
    for _ in range(n_random * 3):
        if len(out) >= n_random:
            break
        name = names[int(rng.integers(len(names)))]
        lo, hi = genes[name]["vcf_comparable"]
        kind = int(rng.choice([0, 0, 0, 0, 1, 2, 3, 4, 5]))
        add(name, int(rng.integers(lo, hi + 1)), kind, ["SNV", "1-base insertion"][kind] if kind < 2 else f"{kind - 1}-base deletion")
    unmoved = sorted(n for n in names if genes[n]["mask"] == genes[n]["spliceai"])
    for name in rng.choice(unmoved, min(n_random // 4, len(unmoved)), replace=False):
        s, e = genes[name]["spliceai"]
        add(str(name), s + int(rng.integers(0, 60)), int(rng.integers(0, 3)), "within 60 bases of a gene start")
        add(str(name), e - int(rng.integers(0, 60)), 0, "within 60 bases of a gene end")
    for name, contig, pos, refb, alt in PUBLIC_REGRESSION[build]:
        out.append((name, contig, pos, refb, alt, "public regression variant"))
    return out


def check_variant(annotator, db, genes, name, chrom, pos, ref, alt, strict=True):
    """Delta-score check of one variant for one gene: dict, or None when the
    variant's window isn't inside the gene's VCF-comparable stretch (with
    strict=False: inside the stored region; for information only)."""
    g = genes.get(name)
    stretch = g and (g["vcf_comparable"] if strict else [g["start"], g["end"]])
    if not stretch:
        return None
    lo, hi = stretch
    if not (lo <= pos - HALF and pos + len(ref) - 1 + HALF <= hi):
        return None
    out, res = delta_scores(annotator, chrom, pos, ref, alt, masked=False)
    if name not in res:
        return None
    live = res[name]
    stored, plan = db.read_gene(g["chrom"], name)
    window = stored[pos - HALF - plan["start"]:pos + HALF + 1 - plan["start"]]
    y_ref_db = np.zeros((1, COV, 3), np.float32)
    y_ref_db[0, :, 1:] = window
    row = {"gene": name, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "spliceai": out}
    for mode in ("raw", "masked"):
        ds_live, dp_live, series = scores_from(live["y_ref"], live["y_alt"], live["dist_exon"], mode == "masked")
        ds_db, dp_db, _ = scores_from(y_ref_db, live["y_alt"], live["dist_exon"], mode == "masked")
        problems, flips = compare(ds_live, dp_live, series, ds_db, dp_db)
        row[mode] = {"live": [round(v, 6) for v in ds_live], "db": [round(v, 6) for v in ds_db],
                     "dp_live": dp_live, "dp_db": dp_db, "max_abs_diff": max(abs(a - b) for a, b in zip(ds_live, ds_db)),
                     "problems": problems, "flips": flips}
    return row


def main(argv=None):
    p = argparse.ArgumentParser(description="Validate a reference score database and write its validation set.")
    p.add_argument("--db", required=True, help="Database folder")
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True, help="Folder for the report and the validation set")
    p.add_argument("--refseq-dir", default=None)
    p.add_argument("--mane-dir", default=None)
    p.add_argument("--variants", type=int, default=300, help="Random variants for the delta-score check")
    p.add_argument("--seed", type=int, default=56)
    p.add_argument("--extra-variants", help="Local TSV gene, chrom, pos, ref, alt (e.g. patient cases): checked "
                                            "too, but reported only in --extra-report, never in --out")
    p.add_argument("--extra-report", help="Where to write the --extra-variants results")
    args = p.parse_args(argv)

    import h5py

    from . import predict, refdb, transcripts
    from .mane import find_cached_summary
    from .normalize import resolve_contig

    os.makedirs(args.out, exist_ok=True)
    db = refdb.Database(args.db)
    db.settings = refdb.read_json(db.path("build.json"))
    manifest = refdb.read_json(db.path("manifest.json"))
    build, genes = manifest["build"], manifest["genes"]
    rng = np.random.default_rng(args.seed)
    report = {"build": build, "database": os.path.abspath(args.db), "started": time.strftime("%Y-%m-%d %H:%M:%S"),
              "tolerance": TOLERANCE, "seed": args.seed}
    t_all = time.time()

    def log(line):
        print(line, flush=True)

    # --- 1. coverage ------------------------------------------------------------
    predictor = predict.Predictor(args.fasta, build)
    fasta = predictor.annotator.ref_fasta
    contigs = set(fasta.keys())
    annotation = {g.name: g for g in transcripts.load_genes(build, args.refseq_dir, find_cached_summary(args.mane_dir))}
    problems = []
    skipped = set(manifest["skipped_genes"])
    missing = sorted(set(annotation) - set(genes) - skipped)
    unknown = sorted(set(genes) - set(annotation))
    if missing:
        problems.append(f"{len(missing)} annotation genes neither stored nor skipped, e.g. {missing[:5]}")
    if unknown:
        problems.append(f"{len(unknown)} stored genes not in the annotation, e.g. {unknown[:5]}")
    if manifest["unpacked_genes"]:
        problems.append(f"{len(manifest['unpacked_genes'])} genes not packed")
    for f in manifest["files"]:
        path = db.path(f["file"])
        if os.path.getsize(path) != f["size"] or refdb.sha256_of(path) != f["sha256"]:
            problems.append(f"{f['file']}: size or SHA-256 differs from the manifest")
    tx_by_acc = {}
    for g in annotation.values():
        for t in g.transcripts:
            tx_by_acc[(g.name, t.accession)] = t
    plan_mismatch, short, uncovered = [], [], []
    for f in manifest["files"]:
        with h5py.File(db.path(f["file"]), "r") as h5:
            for name in h5:
                plan = json.loads(h5[name].attrs["plan"])
                n = len(fasta[resolve_contig(plan["chrom"], contigs)])
                if name in annotation and plan != refdb.plan_gene(annotation[name], n):
                    plan_mismatch.append(name)
                if any(h5[name][t].shape != (plan["end"] - plan["start"] + 1,) for t in ("acceptor", "donor")):
                    short.append(name)
                for acc in plan["transcripts"]:
                    t = tx_by_acc.get((name, acc))
                    if t is None or not (plan["start"] <= t.start and t.end <= plan["end"]):
                        uncovered.append(f"{name} {acc}")
    for label, items in (("stored plan differs from plan_gene()", plan_mismatch),
                         ("array length differs from the stored region", short),
                         ("transcript not inside the stored region", uncovered)):
        if items:
            problems.append(f"{len(items)} {label}, e.g. {items[:5]}")
    report["coverage"] = {"annotation_genes": len(annotation), "stored": len(genes), "skipped": len(skipped),
                          "files": len(manifest["files"]), "problems": problems}
    log(f"coverage: {len(genes):,} stored + {len(skipped)} skipped of {len(annotation):,} genes, "
        f"{len(manifest['files'])} files; problems: {problems or 'none'}")

    # --- 2. stored values against a fresh prediction ------------------------------------
    sample = _sample_genes(genes, rng)
    fresh_values, rows, t0 = {}, [], time.time()
    for k, (name, stratum) in enumerate(sample, 1):
        g = genes[name]
        stored, plan = db.read_gene(g["chrom"], name)
        record = fasta[resolve_contig(g["chrom"], contigs)]
        fresh = predictor.predict(lambda a, b, r=record: r[a - 1:b].seq, len(record), plan["strand"], *plan["mask"],
                                  plan["start"], plan["end"], chunk=CHECK_CHUNK)
        diff = np.abs(stored - fresh)
        edge = np.concatenate([diff[:5050], diff[-5050:]])
        rows.append({"gene": name, "stratum": stratum, "strand": plan["strand"], "bases": len(fresh),
                     "max_abs_diff": float(diff.max()), "max_abs_diff_first_last_5000": float(edge.max())})
        if len(fresh) <= 200_000:
            fresh_values[name] = (fresh, plan, stratum)
        if k % 25 == 0:
            log(f"  fresh prediction {k}/{len(sample)} genes ({time.time() - t0:.0f} s)")
    worst = max(rows, key=lambda r: r["max_abs_diff"])
    by_stratum = {}
    for r in rows:
        by_stratum[r["stratum"]] = max(by_stratum.get(r["stratum"], 0.0), r["max_abs_diff"])
    report["stored_vs_fresh"] = {"genes": len(rows), "bases": sum(r["bases"] for r in rows),
                                 "max_abs_diff": worst["max_abs_diff"], "worst_gene": worst["gene"],
                                 "max_abs_diff_first_last_5000": max(r["max_abs_diff_first_last_5000"] for r in rows),
                                 "max_by_stratum": by_stratum, "passed": worst["max_abs_diff"] <= TOLERANCE,
                                 "seconds": round(time.time() - t0)}
    log(f"stored vs fresh: {len(rows)} genes, {report['stored_vs_fresh']['bases']:,} bases, max |diff| "
        f"{worst['max_abs_diff']:.3g} ({worst['gene']})")

    # --- 3. delta scores -------------------------------------------------------------
    t0 = time.time()
    annotator = predictor.annotator
    checked, outside = [], []
    for name, chrom, pos, ref, alt, why in _variants(genes, rng, fasta, args.variants, build):
        row = check_variant(annotator, db, genes, name, chrom, pos, ref, alt)
        if row is None:
            where = "gene not in the database" if name not in genes else "outside the VCF-comparable stretch"
            outside.append(f"{name} {chrom}:{pos} {ref}>{alt} ({why}; {where})")
            continue
        row["why"] = why
        checked.append(row)
    flips = [f"{r['gene']} {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} {m}: {f}"
             for r in checked for m in ("raw", "masked") for f in r[m]["flips"]]
    failed = [f"{r['gene']} {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} {m}: {'; '.join(r[m]['problems'])}"
              for r in checked for m in ("raw", "masked") if r[m]["problems"]]
    report["delta_scores"] = {
        "variants": len(checked), "by_kind": {w: sum(r["why"] == w for r in checked) for w in sorted({r["why"] for r in checked})},
        "max_abs_diff": max((r[m]["max_abs_diff"] for r in checked for m in ("raw", "masked")), default=0.0),
        "dp_checked": sum(1 for r in checked for m in ("raw", "masked") for j in range(4) if r[m]["live"][j] >= 0.05),
        "problems": failed, "rounding_flips": flips, "not_checked": outside,
        "passed": not failed, "seconds": round(time.time() - t0)}
    log(f"delta scores: {len(checked)} variants, max |dDS| {report['delta_scores']['max_abs_diff']:.3g}, "
        f"{len(failed)} problems, {len(flips)} rounding flips, {len(outside)} not checked (listed in the report)")

    if args.extra_variants:
        extra = []
        with open(args.extra_variants, encoding="utf-8") as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    name, chrom, pos, ref, alt = line.split()[:5]
                    row = check_variant(annotator, db, genes, name, chrom, int(pos), ref, alt)
                    if row is None and name in genes:
                        row = check_variant(annotator, db, genes, name, chrom, int(pos), ref, alt, strict=False)
                        if row:
                            row["outside_vcf_comparable"] = True
                    where = "gene not in the database" if name not in genes else "outside the stored region"
                    extra.append(row or {"gene": name, "chrom": chrom, "pos": pos, "not_checked": where})
        if args.extra_report:
            with open(args.extra_report, "w", encoding="utf-8") as fh:
                json.dump(extra, fh, indent=1, default=str)
        checked = [r for r in extra if "raw" in r]
        log(f"extra variants: {len(checked)} of {len(extra)} checked "
            f"({sum(1 for r in checked if r.get('outside_vcf_comparable'))} outside the VCF-comparable stretch, for "
            f"information), problems: {sum(1 for r in checked if r['raw']['problems'] or r['masked']['problems'])}")

    # --- validation set ----------------------------------------------------------------
    chosen, per = [], {}
    for name in sorted(fresh_values, key=lambda n: (fresh_values[n][2], n)):
        stratum = fresh_values[name][2]
        if stratum.startswith("size") and per.get(stratum, 0) < 7:
            per[stratum] = per.get(stratum, 0) + 1
            chosen.append(name)
    path = os.path.join(args.out, f"validation-set-{build}.h5")
    with h5py.File(path + ".tmp", "w") as out:
        out.attrs["build"] = build
        out.attrs["precision"] = "float16"
        out.attrs["source_manifest_sha256"] = refdb.sha256_of(db.path("manifest.json"))
        out.attrs["model_sha256"] = json.dumps(manifest["model_sha256"])
        out.attrs["note"] = ("Fresh predictions (predict.py) of genes checked against the database; compare a local "
                             "build with them within 2 x 2^-12 + 1e-5 (both float16).")
        for name in chosen:
            fresh, plan, _ = fresh_values[name]
            grp = out.create_group(name)
            for i, track in enumerate(("acceptor", "donor")):
                d = refdb.encode(fresh[:, i], "float16")
                grp.create_dataset(track, data=d, chunks=(min(len(d), 1 << 16),), compression="gzip",
                                   compression_opts=4, shuffle=True)
            grp.attrs["plan"] = json.dumps(plan)
    refdb.replace(path + ".tmp", path)
    report["validation_set"] = {"file": os.path.basename(path), "genes": len(chosen), "size": os.path.getsize(path),
                                "sha256": refdb.sha256_of(path),
                                "bases": sum(len(fresh_values[n][0]) for n in chosen)}
    log(f"validation set: {len(chosen)} genes, {os.path.getsize(path) / 1e6:.1f} MB -> {path}")

    passed = not problems and report["stored_vs_fresh"]["passed"] and report["delta_scores"]["passed"]
    report.update({"passed": passed, "finished": time.strftime("%Y-%m-%d %H:%M:%S"), "seconds": round(time.time() - t_all)})
    with open(os.path.join(args.out, f"validation-{build}.json"), "w", encoding="utf-8") as fh:
        json.dump({"report": report, "stored_vs_fresh_genes": rows, "delta_score_variants": checked}, fh, indent=1,
                  default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    summary = {"passed": passed, "date": report["finished"], "tolerance": TOLERANCE,
               "stored_vs_fresh": {k: report["stored_vs_fresh"][k] for k in ("genes", "bases", "max_abs_diff")},
               "delta_scores": {k: report["delta_scores"][k] for k in ("variants", "max_abs_diff")}
               | {"problems": len(failed), "rounding_flips": len(flips)},
               "coverage_problems": len(problems), "validation_set": report["validation_set"]}
    db.write_manifest(list(annotation), extra={"validation": summary})
    log(f"{'PASSED' if passed else 'FAILED'} in {report['seconds']} s; summary written to the manifest")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
