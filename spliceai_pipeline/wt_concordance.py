"""Genome-wide wild-type concordance (spec 14 C1): for every gene with a
canonical transcript, the WT prediction from the reference score database
(no model runs; the tab's default parameters), compared with the canonical
transcript and with every curated NM_ of the gene.

Per gene: the canonical exons predicted identically, whether the whole
transcript is reproduced exactly, and how many predicted exons match an exon
of any curated NM_ of the gene. A predicted exon matches an annotated one
when the splice sites it has agree: its first base (an acceptor) unless it is
the first exon, its last base (a donor) unless it is the last -- the
transcript start and end are anchors, not predictions.

Every canonical exon not predicted identically gets one category, the first
that applies:
- "first or last exon edge": it is the transcript's first or last exon;
- "canonical site below the site threshold": SpliceAI gives its acceptor or
  donor less than the site threshold, so no path can use it;
- "alternative exon": the prediction there is annotated in another curated
  NM_ of the gene (its exon matches one, or that NM_ skips this exon too);
- "other".

The per-gene table (TABLE_NAME, shipped in spliceai_pipeline/data) is what the
tab shows for a gene. Run from source:

    .venv\\Scripts\\python -m spliceai_pipeline.wt_concordance --build hg38 \\
        --table spliceai_pipeline\\data\\wt_concordance_hg38.tsv.gz --report report.md
"""
import argparse
import csv
import gzip
import io
import os
import sys
import time

from . import cdna, cdna_analysis, refdb, splice_paths, transcripts
from .mane import find_cached_summary

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
TABLE_NAME = "wt_concordance_{build}.tsv.gz"
FIELDS = ["gene", "accession", "strand", "length", "exons", "identical", "first_identical", "last_identical",
          "whole_exact", "predicted_exons", "predicted_annotated", "first_annotated", "last_annotated", "categories",
          "differing"]
NUMBERS = [f for f in FIELDS if f not in ("gene", "accession", "strand", "categories", "differing")]
CATEGORIES = ("first or last exon edge", "canonical site below the site threshold", "alternative exon", "other")
LENGTH_BINS = [(0, 10e3, "< 10 kb"), (10e3, 30e3, "10-30 kb"), (30e3, 100e3, "30-100 kb"), (100e3, 300e3, "100-300 kb"),
               (300e3, 1e6, "300 kb - 1 Mb"), (1e6, 1e12, ">= 1 Mb")]
EXON_BINS = [(2, 4, "2-3"), (4, 7, "4-6"), (7, 11, "7-10"), (11, 21, "11-20"), (21, 41, "21-40"), (41, 10**9, "> 40")]


def in_order(t):
    return list(t.exons) if t.strand == "+" else list(reversed(t.exons))


def _matches(e, k, n, x, strand):
    """Predicted exon e (number k of n, transcript order) has the splice sites of annotated exon x."""
    first = lambda f: f[0] if strand == "+" else f[1]
    last = lambda f: f[1] if strand == "+" else f[0]
    return (k == 0 or first(e) == first(x)) and (k == n - 1 or last(e) == last(x))


def gene_concordance(gene, values, start, settings=None):
    """The dict of one FIELDS row (categories: {category: count}) for gene,
    whose reference values (n, 2) start at genomic position `start`."""
    return _concordance(gene, values, start, settings)[0]


def _concordance(gene, values, start, settings=None):
    """(gene_concordance()'s dict, the concordance rows, whether each predicted exon is annotated)."""
    settings = settings or cdna_analysis.Settings()
    c = gene.canonical
    anchor = cdna_analysis.anchors(gene)
    path = splice_paths.top_paths(values[:, 0], values[:, 1], start, c.strand, anchor["transcript start"],
                                  anchor["transcript end"], k=1, threshold=settings.threshold,
                                  min_exon=settings.min_exon, min_intron=settings.min_intron)[0]
    ref, pred = in_order(c), list(path.exons)
    rows, _ = cdna.concordance(ref, pred, c.strand)
    others = [in_order(t) for t in gene.transcripts if t.accession.startswith("NM_") and t.accession != c.accession]
    annotated = [x for t in [ref] + others for x in t]
    n = len(pred)
    hit = [any(_matches(e, k, n, x, c.strand) for x in annotated) for k, e in enumerate(pred)]
    same = {r["ref"] for r in rows if r["ref"] is not None and r["status"] == "identical"}
    at = lambda x, kind: float(values[x - start, kind]) if 0 <= x - start < len(values) else 0.0
    categories = {}
    for row in rows:
        if row["ref"] is None or row["status"] == "identical":
            continue
        i = row["ref"] - 1
        e = ref[i]
        acc, don = (e[0], e[1]) if c.strand == "+" else (e[1], e[0])
        if i in (0, len(ref) - 1):
            kind = CATEGORIES[0]
        elif at(acc, 0) < settings.threshold or at(don, 1) < settings.threshold:
            kind = CATEGORIES[1]
        elif _alternative(row, pred, others, e, c.strand):
            kind = CATEGORIES[2]
        else:
            kind = CATEGORIES[3]
        categories[kind] = categories.get(kind, 0) + 1
    row = {"gene": gene.name, "accession": c.accession, "strand": c.strand, "length": c.end - c.start + 1,
           "exons": len(ref), "identical": len(same), "first_identical": int(1 in same),
           "last_identical": int(len(ref) in same), "whole_exact": int(len(same) == len(ref) and n == len(ref)),
           "predicted_exons": n, "predicted_annotated": sum(hit), "first_annotated": int(hit[0]),
           "last_annotated": int(hit[-1]), "categories": categories,
           # the canonical exons not predicted identically, e.g. "7,12" (the tab names them)
           "differing": [k for k in range(1, len(ref) + 1) if k not in same]}
    return row, rows, hit


def _alternative(row, pred, others, e, strand):
    if row["other"] is not None:
        k = row["other"] - 1
        return any(_matches(pred[k], k, len(pred), x, strand) for t in others for x in t)
    # Skipped: another NM_ spans this exon without containing it.
    return any(min(map(min, t)) < e[0] and e[1] < max(map(max, t))
               and not any(x[0] <= e[1] and e[0] <= x[1] for x in t) for t in others)


def write_table(rows, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(FIELDS)
        for r in rows:
            w.writerow([";".join(f"{k}={v}" for k, v in sorted(r[f].items())) if f == "categories"
                        else ",".join(map(str, r[f])) if f == "differing" else r[f] for f in FIELDS])


def read_table(path):
    """{gene: row} of a table written by write_table (numbers as int)."""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for f in NUMBERS:
                r[f] = int(r[f])
            r["categories"] = {k: int(v) for k, v in (p.split("=") for p in r["categories"].split(";") if p)}
            r["differing"] = [int(k) for k in (r.get("differing") or "").split(",") if k]
            out[r["gene"]] = r
    return out


_loaded = {}


def for_gene(build, name, data_dir=None):
    """The shipped row of gene `name` on `build`, or None (no table, or not in it)."""
    path = os.path.join(data_dir or DATA_DIR, TABLE_NAME.format(build=build))
    if path not in _loaded:
        _loaded[path] = read_table(path) if os.path.isfile(path) else {}
    return _loaded[path].get(name)


def run(build, reader, genes, on_gene=None):
    """Rows for every gene with a canonical transcript that the database has; (rows, skipped [(gene, why)])."""
    rows, skipped = [], []
    for g in genes:
        if g.canonical is None:
            continue
        plan = reader.plan(g.name)
        if plan is None:
            skipped.append((g.name, "not in the database"))
            continue
        try:
            rows.append(gene_concordance(g, reader.values(g.name), plan["start"]))
        except ValueError as exc:
            skipped.append((g.name, str(exc)))
        if on_gene:
            on_gene(len(rows) + len(skipped))
    return rows, skipped


SWEEP = (25, 15, 10, 6, 3)      # minimum internal exon lengths compared (spec 17 item 4)


def sweep(reader, genes, values=SWEEP, on_gene=None):
    """Spec 17 item 4: the WT prediction of every multi-exon gene with each minimum internal exon length, from
    the database. {min_exon: totals} of exons (a), identical (a), internal, internal_identical, exact (c),
    predicted, annotated (b), extra (predicted exons with no canonical counterpart), false_extra (those that match
    no curated NM_ either), short_found[m] (canonical internal exons shorter than m predicted identically);
    and {m: canonical exons shorter than m}, internal and all."""
    totals = {v: dict.fromkeys(("genes", "exons", "identical", "internal", "internal_identical", "exact", "predicted",
                                "annotated", "extra", "false_extra"), 0) | {"short_found": dict.fromkeys(values, 0)}
              for v in values}
    shorter = {v: {"internal": 0, "all": 0} for v in values}
    done = 0
    for g in genes:
        plan = None if g.canonical is None else reader.plan(g.name)
        if plan is None or len(g.canonical.exons) < 2:
            continue
        ref = in_order(g.canonical)
        lengths = [b - a + 1 for a, b in ref]
        for v in values:
            shorter[v]["all"] += sum(n < v for n in lengths)
            shorter[v]["internal"] += sum(n < v for n in lengths[1:-1])
        data = reader.values(g.name)
        for v in values:
            try:
                row, rows, hit = _concordance(g, data, plan["start"], cdna_analysis.Settings(min_exon=v))
            except ValueError:
                continue
            t = totals[v]
            same = {r["ref"] for r in rows if r["ref"] is not None and r["status"] == "identical"}
            extra = [r["other"] for r in rows if r["ref"] is None and r["other"] is not None]
            t["genes"] += 1
            t["exons"] += row["exons"]
            t["identical"] += row["identical"]
            t["internal"] += row["exons"] - 2
            t["internal_identical"] += len(same - {1, row["exons"]})
            t["exact"] += row["whole_exact"]
            t["predicted"] += row["predicted_exons"]
            t["annotated"] += row["predicted_annotated"]
            t["extra"] += len(extra)
            t["false_extra"] += sum(not hit[k - 1] for k in extra)
            for m in values:
                t["short_found"][m] += sum(lengths[k - 1] < m for k in same if 1 < k < row["exons"])
        done += 1
        if on_gene:
            on_gene(done)
    return totals, shorter


def sweep_report(build, totals, shorter, seconds):
    """The comparison of spec 17 item 4 for one build, as Markdown."""
    values = sorted(totals, reverse=True)
    base = totals[max(values)]
    out = [f"### {build}", "", f"{base['genes']:,} multi-exon genes with a canonical transcript, from the database, "
           f"{seconds / 60:.1f} minutes.", "",
           "| Minimum internal exon | Exon identity (a) | Internal exon identity | Whole transcript exact (c) | "
           "Predicted exons annotated (b) | Extra exons | False extra exons | Canonical exons under 25 nt found |",
           "|---|---|---|---|---|---|---|---|"]
    for v in values:
        t = totals[v]
        out.append(f"| {v} nt | {t['identical']:,} of {t['exons']:,} ({t['identical'] / t['exons']:.2%}) | "
                   f"{t['internal_identical'] / t['internal']:.2%} | {t['exact']:,} of {t['genes']:,} "
                   f"({t['exact'] / t['genes']:.2%}) | {t['annotated'] / t['predicted']:.2%} | {t['extra']:,} | "
                   f"{t['false_extra']:,} | {t['short_found'].get(25, 0):,} of {shorter[25]['internal']:,} |")
    out += ["", "| Canonical exons shorter than | Internal exons | All exons (with first and last) |", "|---|---|---|"]
    out += [f"| {v} nt | {shorter[v]['internal']:,} | {shorter[v]['all']:,} |" for v in values]
    return "\n".join(out) + "\n"


def _bin(value, bins):
    return next(label for lo, hi, label in bins if lo <= value < hi)


def report(build, rows, skipped, seconds):
    """The Markdown report of spec 14 C1 for one build."""
    multi = [r for r in rows if r["exons"] >= 2]
    out = io.StringIO()
    say = lambda line="": out.write(line + "\n")
    pct = lambda a, b: f"{a / b:.1%}" if b else "--"
    exons = sum(r["exons"] for r in multi)
    same = sum(r["identical"] for r in multi)
    pred = sum(r["predicted_exons"] for r in multi)
    annotated = sum(r["predicted_annotated"] for r in multi)
    exact = sum(r["whole_exact"] for r in multi)
    say(f"### {build}")
    say()
    say(f"{len(rows):,} genes with a canonical transcript ({len(rows) - len(multi):,} of them single-exon, left out "
        f"below: nothing to splice); {len(skipped):,} skipped; {seconds / 60:.1f} minutes.")
    say()
    say("| Measure | Count | Rate |")
    say("|---|---|---|")
    say(f"| (a) canonical exons predicted identically | {same:,} of {exons:,} | {pct(same, exons)} |")
    say(f"| (b) predicted exons matching an exon of any curated NM_ of the gene | {annotated:,} of {pred:,} | "
        f"{pct(annotated, pred)} |")
    say(f"| (c) genes whose whole transcript is reproduced exactly | {exact:,} of {len(multi):,} | "
        f"{pct(exact, len(multi))} |")
    say()
    say("**By exon position:**")
    say()
    say("| Position | (a) canonical exons identical | (b) predicted exons annotated |")
    say("|---|---|---|")
    total = lambda key: sum(r[key] for r in multi)
    inner_a = total("identical") - total("first_identical") - total("last_identical")
    inner_b = total("predicted_annotated") - total("first_annotated") - total("last_annotated")
    inner_n, inner_p = exons - 2 * len(multi), pred - 2 * len(multi)
    for label, a_, n_, b_, p_ in (("first", total("first_identical"), len(multi), total("first_annotated"), len(multi)),
                                  ("internal", inner_a, inner_n, inner_b, inner_p),
                                  ("last", total("last_identical"), len(multi), total("last_annotated"), len(multi))):
        say(f"| {label} | {a_:,} of {n_:,} ({pct(a_, n_)}) | {b_:,} of {p_:,} ({pct(b_, p_)}) |")
    say()
    for title, key, bins in (("gene length (canonical transcript)", "length", LENGTH_BINS),
                             ("exon count", "exons", EXON_BINS)):
        say(f"**By {title}:**")
        say()
        say("| Bin | Genes | (a) exons identical | (b) predicted exons annotated | (c) whole transcript exact |")
        say("|---|---|---|---|---|")
        for lo, hi, label in bins:
            sel = [r for r in multi if _bin(r[key], bins) == label]
            if not sel:
                continue
            e, s = sum(r["exons"] for r in sel), sum(r["identical"] for r in sel)
            p_, a_ = sum(r["predicted_exons"] for r in sel), sum(r["predicted_annotated"] for r in sel)
            x = sum(r["whole_exact"] for r in sel)
            say(f"| {label} | {len(sel):,} | {pct(s, e)} | {pct(a_, p_)} | {pct(x, len(sel))} |")
        say()
    say("**Canonical exons not predicted identically, by category:**")
    say()
    say("| Category | Exons | Share |")
    say("|---|---|---|")
    total = sum(sum(r["categories"].values()) for r in multi)
    for k in CATEGORIES:
        n = sum(r["categories"].get(k, 0) for r in multi)
        say(f"| {k} | {n:,} | {pct(n, total)} |")
    say()
    if skipped:
        say("Skipped: " + "; ".join(f"{g} ({why})" for g, why in skipped[:20]) + (" ..." if len(skipped) > 20 else ""))
        say()
    return out.getvalue()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", choices=("hg19", "hg38"), required=True)
    p.add_argument("--db", default=None, help="Database folder (default: the app's own)")
    p.add_argument("--refseq-dir", default=None)
    p.add_argument("--mane-dir", default=None)
    p.add_argument("--table", help="Per-gene table to write (.tsv.gz)")
    p.add_argument("--sweep", action="store_true", help="Compare minimum internal exon lengths instead (spec 17 item 4)")
    p.add_argument("--report", help="Markdown report to write")
    a = p.parse_args(argv)
    t0 = time.time()
    reader = refdb.Reader(a.db or os.path.join(refdb.default_root(), a.build))
    genes = transcripts.load_genes(a.build, a.refseq_dir, find_cached_summary(a.mane_dir))
    total = sum(g.canonical is not None for g in genes)

    def progress(n):
        if n % 1000 == 0:
            print(f"{n:,} of {total:,} genes, {time.time() - t0:,.0f} s", flush=True)

    if a.sweep:
        text = sweep_report(a.build, *sweep(reader, genes, on_gene=progress), time.time() - t0)
    else:
        rows, skipped = run(a.build, reader, genes, progress)
        write_table(rows, a.table)
        text = report(a.build, rows, skipped, time.time() - t0)
    print(text)
    if a.report:
        with open(a.report, "w", encoding="utf-8") as fh:
            fh.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
