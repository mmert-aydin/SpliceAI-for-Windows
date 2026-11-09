"""Spec 18: the workbook and the numbers of the Phase 14 benchmark, from benchmark\\results.

    python tests\\benchmark_report.py      (any Python with openpyxl; it needs no model)

Writes SpliceAI_cDNA_benchmark.xlsx and report\\ (report.md, report.html,
summary.json) in benchmark\\, next to the repository. Definitions are in benchmark.py; the two here:
- a DS-only call is spliceogenic at DS >= 0.2; DS is the largest of the four VCF-tab scores;
- confidence intervals are exact (Clopper-Pearson), 95%.
"""
import base64
import collections
import csv
import html
import json
import math
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = Path(os.environ.get("SPLICEAI_BENCHMARK") or REPO.parent / "benchmark")
DS_LINE, DS_HIGH = 0.2, 0.5
CONTEXT = 5000
TISSUE = re.compile(r"muscle|retina|photoreceptor|brain|neuron|cardiac|heart|testis|liver|kidney", re.I)
CLASSES = ("insertion", "deletion", "duplication", "inversion", "complex", "pseudoexon SNV", "other SNV")


def tail_low(k, n, p):
    """P(X >= k) for Binomial(n, p)."""
    return math.fsum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))


def clopper_pearson(k, n, alpha=0.05):
    """Exact two-sided interval of k successes in n."""
    if n == 0:
        return None, None

    def solve(f, want):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(mid) < want:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    low = 0.0 if k == 0 else solve(lambda p: tail_low(k, n, p), alpha / 2)
    high = 1.0 if k == n else solve(lambda p: 1 - (1 - tail_low(k + 1, n, p)), 1 - alpha / 2)
    return low, high


def load():
    cases = {}
    with open(OUT / "cases.jsonl", encoding="utf-8") as fh:
        for line in fh:
            c = json.loads(line)
            cases[c["id"]] = c
    ds = {}
    path = OUT / "ds_output.tsv"
    if path.exists():
        with open(path, encoding="utf-8", newline="") as fh:
            for d in csv.DictReader(fh, delimiter="	"):
                ds.setdefault(d["chrom"], []).append(d)
    rows = []
    for case_id, c in cases.items():
        p = OUT / "results" / f"{case_id}.json"
        r = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {"error": "not run"}
        d = ds_row(c, ds) if c.get("vcf") else {}
        score, why = None, ""
        if d:
            try:
                score = float(d.get("max_score", ""))
            except ValueError:
                why = d.get("source", "")
                why = "SpliceAI gave no score" if why in ("", "live", "none") else why
        else:
            why = "no VCF record (several events, or too long)" if not c.get("vcf") else "no row: too long for the pipeline's output"
        rows.append(score_case(c, r, score, why))
    return rows


def ds_row(c, ds):
    """The VCF tab's row of a case: the pipeline left-aligns an indel, so the same shape nearest the position, and of the
    case's gene when a variant has rows for several."""
    chrom, pos, ref, alt = c["vcf"]
    near = [d for d in ds.get(chrom, ()) if abs(int(d["pos"]) - pos) <= 300 and len(d["ref"]) == len(ref) and len(d["alt"]) == len(alt)]
    if not near:
        return {}
    best = min(abs(int(d["pos"]) - pos) for d in near)
    near = [d for d in near if abs(int(d["pos"]) - pos) == best]
    return next((d for d in near if d.get("gene") == c["gene"]), near[0])


def score_case(c, r, ds, why):
    out = {"case": c, "result": r, "ds": ds, "ds_why": why, "error": r.get("error")}
    if r.get("error"):
        return out
    truth = not c.get("negative")
    comp = r["compare"]
    out.update(truth=truth, call=r["spliceogenic"], strong=r["spliceogenic_strong"], yes=r["spliceogenic"] == truth,
               ds_yes=None if ds is None else (ds >= DS_LINE) == truth, ds_high_yes=None if ds is None else (ds >= DS_HIGH) == truth,
               event1=comp["event_rank1"], event5=comp["event_top"], exact1=comp["exact_rank1"], exact5=comp["exact_top"])
    return out


def group(rows):
    """[(name, rows)] the classes, the MET cluster apart, native and minigene apart."""
    plain = [r for r in rows if not r["error"]]
    out = []
    for cls in CLASSES:
        some = [r for r in plain if r["case"]["class"] == cls and not r["case"].get("cluster") and not r["case"].get("negative")]
        out.append((cls, some))
    out.append(("deletion: MET exon 14 cluster", [r for r in plain if r["case"].get("cluster")]))
    core = [r for r in plain if r["case"]["class"] in ("insertion", "deletion", "duplication", "inversion", "complex")
            and not r["case"].get("cluster") and not r["case"].get("negative")]
    out.append(("all structural and mid-size, without the MET cluster", core))
    out.append(("  the same, native RNA", [r for r in core if r["case"]["native"]]))
    out.append(("  the same, minigene", [r for r in core if not r["case"]["native"]]))
    out.append(("  the same, the variant does not itself remove the exon", [r for r in core if not r["case"].get("removes_exon")]))
    for name, lo, hi in (("under 10 bases", 0, 9), ("10 to 49 bases", 10, 49), ("50 to 99 bases", 50, 99), ("100 to 999 bases", 100, 999),
                         ("1,000 bases or more", 1000, 10 ** 9)):
        out.append((f"  the same, {name}", [r for r in core if lo <= size_of(r["case"]) <= hi]))
    snv = [r for r in plain if r["case"]["class"] in ("pseudoexon SNV", "other SNV") and not r["case"].get("negative")]
    out.append(("SNVs, native RNA", [r for r in snv if r["case"]["native"]]))
    out.append(("SNVs, minigene", [r for r in snv if not r["case"]["native"]]))
    return out


def size_of(c):
    sizes = []
    for g in c["g"].split("\n"):
        m = re.search(r"g\.(\d+)(?:_(\d+))?([a-z]+)?([ACGT]*)", g)
        if not m:
            continue
        span = int(m.group(2) or m.group(1)) - int(m.group(1)) + 1
        sizes.append(len(m.group(4)) if m.group(3) == "ins" else max(span, len(m.group(4))))
    return max(sizes, default=1)


def level(rows, key):
    known = [r for r in rows if r.get(key) is not None]
    k = sum(bool(r[key]) for r in known)
    return k, len(known)


METRICS = (("Spliceogenic, engine", "yes"), ("Spliceogenic, DS >= 0.2", "ds_yes"), ("Spliceogenic, DS >= 0.5", "ds_high_yes"),
           ("Event, rank 1", "event1"), ("Event, top 5", "event5"), ("Exact junctions, rank 1", "exact1"),
           ("Exact junctions, top 5", "exact5"))


def paired(rows):
    """On the cases with a DS: (engine right, cases, engine right and DS wrong, DS right and engine wrong)."""
    both = [r for r in rows if r.get("ds_yes") is not None]
    return (sum(r["yes"] for r in both), len(both), sum(r["yes"] and not r["ds_yes"] for r in both),
            sum(r["ds_yes"] and not r["yes"] for r in both))


def text(k, n):
    if not n:
        return "-"
    lo, hi = clopper_pearson(k, n)
    return f"{k}/{n} = {100 * k / n:.0f}% ({100 * lo:.0f}-{100 * hi:.0f})"


def why_missed(r):
    c, x = r["case"], r["result"]
    d = x.get("distance")
    if c.get("negative"):
        return "other: the engine calls a change where the experiment found none"
    gained = any(k in o for o in c["outcomes"] for k in ("exons", "donor", "acceptor", "pe_end", "pe_start"))
    if d is not None and d > CONTEXT and gained:      # a skipped exon's far junction says nothing of the model's reach
        return f"beyond the model's 5,000-base context: the variant is {d:,} bases from the documented site"
    if not c["native"]:
        return "minigene context: the assay's outcome may not be the gene's"
    if not x.get("local_same", True):
        return "WT discordant in that region"
    if c.get("partial") and x["top"] and x["top"][0]["max_delta"] is None:
        return "other: a partial effect, the normal transcript keeps rank 1"
    if TISSUE.search(c.get("tissue") or "") and not any(s[3] - s[2] >= 0.1 for s in x.get("signal") or []):
        return f"tissue-specific outcome possible: seen in {c['tissue']}, no gain at the site"
    sites = x.get("signal") or []
    if sites and max(s[3] - s[2] for s in sites) >= 0.1:
        return "other: SpliceAI raises the documented site, but the best paths don't use it"
    if sites:
        return "other: SpliceAI gives the documented site no gain"
    return "other: the source gives no site to read the signal at"


def write_workbook(rows, groups):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill

    wb = Workbook()
    head_font, fill = Font(bold=True, color="FFFFFF"), PatternFill("solid", fgColor="1F3A5F")

    def sheet(ws, header, data, widths=None):
        ws.append(header)
        for cell in ws[1]:
            cell.font, cell.fill = head_font, fill
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        for row in data:
            ws.append(["" if v is None else v for v in row])
        ws.freeze_panes = "A2"
        for i, w in enumerate(widths or (), 1):
            ws.column_dimensions[ws.cell(1, i).column_letter].width = w

    ws = wb.active
    ws.title = "Summary"
    ws.append(["SpliceAI Sequence to cDNA engine against RNA-confirmed variants. Agent-curated, not yet verified: "
               "every record awaits the user's check. Research use only."])
    ws.append([f"{sum(not r['error'] for r in rows)} cases run; 95% exact (Clopper-Pearson) intervals. "
               "The yes/no level counts spliceogenic variants only, so it is a sensitivity."])
    ws.append([])
    ws.append(["Group", "Cases"] + [m for m, _ in METRICS] + ["Engine on the cases DS can score", "Engine yes, DS no", "DS yes, engine no"])
    for cell in ws[4]:
        cell.font, cell.fill = head_font, fill
        cell.alignment = Alignment(wrap_text=True, vertical="top")
    for name, some in groups:
        k, n, only, other = paired(some)
        ws.append([name, len(some)] + [text(*level(some, key)) for _, key in METRICS] + [text(k, n), only, other])
    negatives = [r for r in rows if not r["error"] and r["case"].get("negative")]
    ws.append([])
    ws.append(["Variants an experiment found to have no splicing effect", len(negatives),
               f"engine right on {sum(r['yes'] for r in negatives)} of {len(negatives)}"])
    ws.column_dimensions["A"].width = 52
    for col in "CDEFGHIJ":
        ws.column_dimensions[col].width = 22

    header = ["ID", "Status", "Gene", "Variant as written", "hg38 (g.)", "Assay", "Tissue", "Native RNA", "Partial effect", "Exon removed by the variant",
              "Reported cDNA effect (paraphrase)", "Where in the source", "First author", "Journal", "Year", "PMID", "PMCID", "DOI", "Title", "Curation note",
              "Rank-1 transcript", "Rank-1 change vs WT", "Top-5 changes", "Relative score", "Gap to rank 2", "Max delta", "Weak evidence",
              "Gene WT concordance", "Local concordance", "Frame", "NMD", "DS (VCF tab)", "DS note",
              "Spliceogenic: truth", "Spliceogenic: engine", "Correct", "DS >= 0.2 correct", "Event rank 1", "Event top 5", "Exact rank 1", "Exact top 5"]
    done = [r for r in rows if not r["error"]]
    for cls in CLASSES:
        ws = wb.create_sheet(cls[:31])
        data = []
        for r in done:
            c, x = r["case"], r["result"]
            if c["class"] != cls:
                continue
            t = x["top"]
            data.append([c["id"], c["status"], c["gene"], c["written"], c["g"], c["assay"], c["tissue"], c["native"], c.get("partial"), bool(c.get("removes_exon")),
                         c["effect"], c["where"], c.get("first_author"), c.get("journal"), c.get("year"), c.get("pmid"), c.get("pmcid"), c.get("doi"), c.get("title"),
                         "; ".join(y for y in (c.get("note"), c.get("numbering"), c.get("cluster")) if y),
                         x["transcript"], t[0]["change_vs_wt"][:300], " || ".join(y["change_vs_wt"][:120] for y in t), t[0]["relative"], x["gap"], t[0]["max_delta"], x["weak"],
                         f"{x['gene_identical']}/{x['gene_exons']}" + (", whole transcript" if x["gene_whole_exact"] else ""), " ".join(x["local"])[:300],
                         t[0]["frame"], t[0]["nmd"], r["ds"], r["ds_why"], r["truth"], r["call"], r["yes"], r["ds_yes"],
                         r["event1"], r["event5"], "n/a" if r["exact1"] is None else r["exact1"], "n/a" if r["exact5"] is None else r["exact5"]])
        sheet(ws, header, data, [38, 18, 10, 34, 34])

    ws = wb.create_sheet("Failures")
    data = []
    for r in done:
        c, x = r["case"], r["result"]
        missed = not r["event5"] if not c.get("negative") else not r["yes"]
        near = not missed and r["exact5"] is False
        if missed or near:
            data.append([c["id"], c["class"], "exact junctions" if near else "event", c["gene"], c["written"], c["effect"][:200], x["top"][0]["change_vs_wt"][:200],
                         "right event, another base" if near else why_missed(r), x.get("distance"),
                         "; ".join(f"{s[0]} {s[1]:,}: {s[2]} -> {s[3]}" for s in x.get("signal") or []), c.get("cluster", "")])
    sheet(ws, ["ID", "Class", "Missed at (top 5)", "Gene", "Variant", "Documented", "Rank 1 predicted", "Likely reason", "Distance to the site (bases)",
               "SpliceAI at the documented sites (reference -> sample)", "Cluster"], data, [38, 14, 14, 10, 30, 50, 50, 50, 12, 44, 12])

    ws = wb.create_sheet("Excluded")
    with open(OUT / "excluded.tsv", encoding="utf-8", newline="") as fh:
        ex = list(csv.DictReader(fh, delimiter="\t"))
    sheet(ws, ["Gene", "Variant as written", "Class", "Why excluded", "Reported effect", "Source"],
          [[e["gene"], e["written"], e["cls"], e["reason"], e["effect"], e["src"]] for e in ex], [14, 44, 12, 60, 40, 14])

    ws = wb.create_sheet("Sources")
    log = json.loads((OUT / "sources" / "meta" / "search_log.json").read_text(encoding="utf-8"))
    sheet(ws, ["Source", "Search name", "Date", "Hits", "Terms"],
          [[e["source"], e["name"], e["date"], e["hits"], e["query"]] for e in log]
          + [["DBASS5, DBASS3 (CSV download)", "whole tables", "2026-10-31", "693 + 381 records", "https://dbass.soton.ac.uk"],
             ["NCBI efetch db=pmc", "full texts", "2026-10-31", "9,297 papers", "the abstract-selected papers, then the focused searches"],
             ["VariantValidator REST", "HGVS to hg38", "2026-10-31", "", "rest.variantvalidator.org"],
             ["UCSC liftOver chains", "hg19, hg18 to hg38", "2026-10-31", "", "hgdownload.soe.ucsc.edu"]], [30, 16, 12, 16, 120])
    wb.save(OUT / "SpliceAI_cDNA_benchmark.xlsx")


# --- the report ---------------------------------------------------------------------------------------------------

BANNER = ("Research use only: predictions from SpliceAI, not measurements. All records are agent-curated and not yet "
          "verified; nothing here is a clinical result.")


def render(blocks):
    """(Markdown, a self-contained HTML page) of [("h1"/"h2"/"p", text), ("list", [..]), ("table", head, rows), ("img", path, caption)]."""
    md, page = [], []
    for b in blocks:
        kind = b[0]
        if kind in ("h1", "h2"):
            md.append("#" * int(kind[1]) + " " + b[1])
            page.append(f"<{kind}>{html.escape(b[1])}</{kind}>")
        elif kind == "p":
            md.append(b[1])
            page.append(f"<p>{html.escape(b[1])}</p>")
        elif kind == "list":
            md.append("\n".join(f"- {x}" for x in b[1]))
            page.append("<ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in b[1]) + "</ul>")
        elif kind == "table":
            md.append("\n".join(["| " + " | ".join(b[1]) + " |", "|" + "---|" * len(b[1])]
                                + ["| " + " | ".join(str(x) for x in r) + " |" for r in b[2]]))
            page.append("<table><tr>" + "".join(f"<th>{html.escape(x)}</th>" for x in b[1]) + "</tr>"
                        + "".join("<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in r) + "</tr>" for r in b[2]) + "</table>")
        elif kind == "img":
            md.append(f"![{b[2]}](../images/{b[1].name})\n\n*{b[2]}*")
            data = base64.b64encode(b[1].read_bytes()).decode("ascii")
            page.append(f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(b[2])}"><figcaption>{html.escape(b[2])}</figcaption></figure>')
    style = ("body{font-family:Segoe UI,Arial,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#1c2733;line-height:1.5}"
             "table{border-collapse:collapse;margin:1em 0;font-size:13px}th,td{border:1px solid #c9d2dc;padding:4px 8px;text-align:left;vertical-align:top}"
             "th{background:#1f3a5f;color:#fff}img{max-width:100%;border:1px solid #c9d2dc}figure{margin:1.5em 0}figcaption{font-size:13px;color:#46566a}"
             "h1{color:#1f3a5f}h2{color:#1f3a5f;margin-top:1.6em}")
    return "\n\n".join(md) + "\n", ("<!DOCTYPE html><html><head><meta charset='utf-8'><title>SpliceAI Sequence to cDNA benchmark</title>"
                                     f"<style>{style}</style></head><body>" + "".join(page) + "</body></html>")


def write_report(rows, groups):
    by = dict(groups)
    core, all_rows = by["all structural and mid-size, without the MET cluster"], [r for r in rows if not r["error"]]
    pct = lambda kn: f"{100 * kn[0] / kn[1]:.0f}%" if kn[1] else "-"
    row = lambda name: [name.strip(), len(by[name])] + [text(*level(by[name], key)) for _, key in METRICS]
    k1, n1 = level(core, "event1")
    k5, n5 = level(core, "event5")
    x1, xn = level(core, "exact1")
    x5, xn5 = level(core, "exact5")
    yes, yes_n = level(core, "yes")
    both, both_n, only, other = paired(core)
    ds, ds_n = level(core, "ds_yes")
    negatives = [r for r in all_rows if r["case"].get("negative")]
    with open(OUT / "excluded.tsv", encoding="utf-8", newline="") as fh:
        excluded = list(csv.DictReader(fh, delimiter="\t"))
    no_seq = sum("unknown inserted sequence" in e["reason"] for e in excluded)
    failures = collections.Counter()
    for r in all_rows:
        if not r["case"].get("negative") and not r["event5"]:
            why = why_missed(r)
            failures[why if why.startswith("other") else why.split(":")[0]] += 1
    figures = [line.strip() for line in (OUT / "figures.txt").read_text().splitlines() if line.strip()]
    missed = [r for r in core if not r["event5"]]
    no_gain = sum("no gain" in why_missed(r) for r in missed)
    plain = (f"We tested whether the Sequence to cDNA engine, which runs SpliceAI over a sample's whole gene and rebuilds its most likely "
             f"transcripts, can foresee how real variants change RNA. We gathered {len(all_rows)} variants with experimental RNA evidence, "
             f"{len(core)} of them structural or mid-size and not part of one repeated site. For those {len(core)}, the engine named the right "
             f"kind of change at the right exon or intron as its first answer in {pct((k1, n1))} of cases and within its top five in {pct((k5, n5))}. "
             f"It placed the splice junctions to the exact base as its first answer in {pct((x1, xn))} of the cases where the paper gives them, "
             f"and within its top five in {pct((x5, xn5))}. It is an aid for planning RNA tests, not a replacement for them: it missed {len(missed)} of "
             f"the {len(core)} even in its top five, and in {no_gain} of those SpliceAI showed no gain at the site the patient's RNA uses.")
    blocks = [("h1", "SpliceAI Sequence to cDNA engine against RNA-confirmed variants"), ("p", BANNER),
              ("h2", "Plain-language summary"), ("p", plain),
              ("h2", "What was tested"),
              ("list", [f"{len(all_rows)} cases run with the default settings, no tuning: " + ", ".join(
                  f"{n} {c}" for c, n in collections.Counter(r["case"]["class"] for r in all_rows).most_common()) + ".",
                        f"{sum(bool(r['case'].get('cluster')) for r in all_rows)} deletions are one site, the MET exon 14 acceptor, somatic, from one paper's table; they are shown apart so they don't carry the deletion class.",
                        f"{sum(bool(r['case'].get('dbass')) for r in all_rows)} cases come from DBASS, placed on hg38 by matching their bases, not by trusting a coordinate.",
                        f"{len(negatives)} cases are variants an experiment found to have no splicing effect.",
                        f"{len(excluded)} variants were excluded and listed, {no_seq} of them because the inserted sequence is not known (most Alu, SVA and LINE-1 reports)."]),
              ("h2", "Results"),
              ("p", "Each cell: right / cases = percent (exact 95% interval). Spliceogenic counts documented spliceogenic variants only, so it is a sensitivity. "
                    "Event: the same kind of change at the same exon or intron, and not already in the WT prediction. Exact: every new junction of a documented outcome, base for base; "
                    "cases whose source gives no junction are left out of that column."),
              ("table", ["Group", "Cases"] + [m for m, _ in METRICS], [row(name) for name, _ in groups]),
              ("h2", "Against the DS-only call"),
              ("p", f"On the {both_n} structural and mid-size cases the VCF tab can score, the engine calls {both} spliceogenic and DS >= 0.2 calls {ds}: the engine alone is right on {only}, "
                    f"DS alone on {other}. Another {len(core) - both_n} of the {len(core)} get no DS from the VCF tab, mostly because the reference allele is too long for SpliceAI; "
                    "the engine scores them. For SNVs DS >= 0.2 is more sensitive than the engine's rank 1, which asks for more than a change in one score: a whole transcript that outranks the normal one. "
                    "DS names no transcript, so it has no event or junction level to compare."),
              ("h2", "Why cases were missed"),
              ("p", "Likely reasons, by rule, for the cases with no matching event in the top five; each case is in the Failures sheet."),
              ("table", ["Likely reason", "Cases"], [[k, v] for k, v in failures.most_common()]),
              ("h2", "Examples"), ("p", "From the viewer's own exports: the exon schematic and the stacked view at the event.")]
    cases = {r["case"]["id"]: r for r in all_rows}
    for case_id in figures:
        r = cases[case_id]
        c, t = r["case"], r["result"]["top"][0]
        good = "predicted" if r["event1"] else "missed"
        blocks.append(("p", f"{c['gene']} {c['written'][:70]} ({good}). Documented: {c['effect'][:220]} Rank 1: {t['change_vs_wt'][:140]}."))
        for name in sorted(p.name for p in (OUT / "images").glob(f"{case_id}_*.png")):
            blocks.append(("img", OUT / "images" / name, f"{c['gene']}: {'exon schematic' if 'schematic' in name else 'stacked view'}"))
    blocks += [("h2", "Limits of this benchmark"),
               ("list", ["Every record is agent-curated and not yet verified. Sources are cited by section, table or figure with a paraphrase; the papers' sentences are not reproduced.",
                         "Published cases favour clear effects, so these figures can overstate results on a new, unselected variant.",
                         f"Specificity isn't measured: only {len(negatives)} no-effect variants are here. A proper estimate needs many more.",
                         f"The insertion and inversion classes are small ({len(by['insertion'])} and {len(by['inversion'])} cases), because most mobile-element reports don't give the inserted bases; their intervals are very wide.",
                         "Exon-removing deletions are easy: the exon's bases are gone. The row without them is the fairer figure.",
                         "SpliceAI has no tissue information. Sample RNA comes from blood, fibroblasts, muscle, retina-like cells and minigenes; a miss may be a true difference between tissues.",
                         "DBASS records don't give assay and tissue case by case, and most of its pseudoexons are known by one end only, so the exact level asks for that base alone.",
                         "Some junctions were derived from reported sizes, and marked so in the curation note. Supplementary files were not read, and papers with no open full text were read at abstract level."])]
    md, page = render(blocks)
    (OUT / "report" / "report.md").write_text(md, encoding="utf-8")
    (OUT / "report" / "report.html").write_text(page, encoding="utf-8")


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    rows = load()
    groups = group(rows)
    write_workbook(rows, groups)
    os.makedirs(OUT / "report", exist_ok=True)
    write_report(rows, groups)
    summary = {name: {"cases": len(some), **{key: level(some, key) for _, key in METRICS}} for name, some in groups}
    summary["errors"] = [[r["case"]["id"], r["error"]] for r in rows if r["error"]]
    summary["paired"] = {name: paired(some) for name, some in groups}
    negatives = [r for r in rows if not r["error"] and r["case"].get("negative")]
    summary["negatives"] = [sum(r["yes"] for r in negatives), len(negatives)]
    (OUT / "report" / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    for name, some in groups:
        print(f"{name[:56]:56} {len(some):4} " + " | ".join(text(*level(some, key)) for _, key in METRICS))
    print("errors", len(summary["errors"]), summary["errors"][:5])


if __name__ == "__main__":
    main()
