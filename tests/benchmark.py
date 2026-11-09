"""Spec 18 (Phase 14): the Sequence to cDNA engine against variants whose splicing
effect was shown on RNA.

    .venv\\Scripts\\python tests\\benchmark.py check      every case read by the tab's own check, no model
    .venv\\Scripts\\python tests\\benchmark.py run        the engine, default settings, one JSON per case; resumes
    .venv\\Scripts\\python tests\\benchmark.py ds         the VCF tab's DS for every case it can score
    .venv\\Scripts\\python tests\\benchmark.py images ID...  the viewer's exports of the cases named

The case list and everything made from it stay in benchmark\\, next to the
repository and never in it (SPLICEAI_BENCHMARK): the records quote third-party
sources. cases.jsonl has one case a line:

    id, class, gene, build ("hg38"), events (the "From coordinates" box's text),
    vcf ([chrom, pos, ref, alt] or null), outcomes (below), partial, native, ...

An outcome is one documented aberrant transcript, as changes to the canonical
one. Positions are "c.123+45" on the canonical transcript or "g.123456" on hg38;
an exon is its canonical number, or any position inside it:

    skip [exon, ...]         the exons are missing
    skip_range [a, b]        every exon between these two positions is missing
    exons [[first, last]]    an extra exon, its first and last base in transcript order
    donor [position]         the last base of the exon that holds it, or before it
    acceptor [position]      the first base of the exon that holds it, or after it
    retain [exon]            the intron after that exon is kept
    dup [exon, ...]          a tandem copy of those exons
    inserted N               an exon of N bases made of inserted sequence (0: any size)
    site [[kind, position]]  the same, where the paper gives no junction: pseudoexon, donor or acceptor
                             at that place; the exact level is then not applicable (None)
    pe_end, pe_start p       a pseudoexon known by one end: its last or first base. The exact level
                             asks for that base only

The three levels: spliceogenic (rank 1 differs from the WT prediction at a splice
site), the event (the same kind of change at the same exon or intron, and new: a
change the WT prediction already has isn't the variant's) and the exact junctions
(every new junction of the outcome, base for base). Any outcome of a case counts.
See compare().
"""
import csv
import json
import os
import sys
import time
from pathlib import Path

import helpers
from helpers import check

OUT = Path(os.environ.get("SPLICEAI_BENCHMARK") or helpers.REPO.parent / "benchmark")
CASES, RESULTS, IMAGES = OUT / "cases.jsonl", OUT / "results", OUT / "images"
TOP = 5
STRONG = 0.5            # the tab's weak-evidence line
BUILD = "hg38"


# --- the comparison: no model, no Qt ---------------------------------------------------------------------

def first_base(e, plus):
    return e[0] if plus else e[1]


def last_base(e, plus):
    return e[1] if plus else e[0]


def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def changes(exons, canonical, plus):
    """{(kind, where)}: how a chain of exons (transcript order, reference
    positions, None for an exon of inserted bases) differs from the canonical
    chain. kind: skip, donor, acceptor (exon number), retain, pseudoexon,
    inserted (intron number, the one before canonical exon n + 1)."""
    out = set()
    real = [e for e in exons if e]
    for k, c in enumerate(canonical, 1):
        if not any(overlaps(e, c) for e in real):
            out.add(("skip", k))
    for i, e in enumerate(exons):
        if e is None:
            before = next((x for x in reversed(exons[:i]) if x), None)
            n = max((k for k, c in enumerate(canonical, 1) if before and (
                first_base(c, plus) * (1 if plus else -1) <= last_base(before, plus) * (1 if plus else -1))), default=0)
            out.add(("inserted", n))
            continue
        hit = [k for k, c in enumerate(canonical, 1) if overlaps(e, c)]
        if not hit:
            n = sum(last_base(c, plus) * (1 if plus else -1) < first_base(e, plus) * (1 if plus else -1) for c in canonical)
            if 0 < n < len(canonical):
                out.add(("pseudoexon", n))
            continue
        out |= {("retain", k) for k in hit[:-1]}
        if first_base(e, plus) != first_base(canonical[hit[0] - 1], plus) and hit[0] > 1 and e is not real[0]:
            out.add(("acceptor", hit[0]))
        if last_base(e, plus) != last_base(canonical[hit[-1] - 1], plus) and hit[-1] < len(canonical) and e is not real[-1]:
            out.add(("donor", hit[-1]))
    return out


def junctions(exons, plus):
    """{(donor, acceptor)} of a chain, reference positions; none across an exon of inserted bases."""
    return {(last_base(a, plus), first_base(b, plus)) for a, b in zip(exons, exons[1:]) if a and b}


def expected_chain(outcome, canonical, plus, position):
    """(exon chain, {(kind, exon)} of copies and inserted exons) of an outcome.
    position(text) gives a "c." or "g." position's genomic position."""
    sign = 1 if plus else -1

    def number(x):
        if isinstance(x, int):
            return x
        p = position(x)
        return next(k for k, c in enumerate(canonical, 1) if c[0] <= p <= c[1])

    def holding(p, upstream):
        """The exon with p in it, else the nearest on that side of p."""
        for k, c in enumerate(canonical, 1):
            if c[0] <= p <= c[1]:
                return k
        side = [k for k, c in enumerate(canonical, 1) if (last_base(c, plus) * sign < p * sign if upstream
                                                         else first_base(c, plus) * sign > p * sign)]
        return side[-1] if upstream else side[0]

    chain = {k: list(c) for k, c in enumerate(canonical, 1)}
    extra = [tuple(sorted((position(a), position(b)))) for a, b in outcome.get("exons", ())]
    for x in outcome.get("donor", ()):
        p = position(x)
        k = holding(p, True)
        chain[k][1 if plus else 0] = p
    for x in outcome.get("acceptor", ()):
        p = position(x)
        k = holding(p, False)
        chain[k][0 if plus else 1] = p
    for x in outcome.get("retain", ()):
        k = number(x)
        chain[k + 1] = [min(chain[k][0], chain[k + 1][0]), max(chain[k][1], chain[k + 1][1])]
        chain.pop(k)
    for x in outcome.get("skip", ()):
        chain.pop(number(x), None)
    if "skip_range" in outcome:
        lo, hi = sorted(position(x) for x in outcome["skip_range"])
        for k, c in enumerate(canonical, 1):
            if lo <= c[0] and c[1] <= hi:
                chain.pop(k, None)
    exons = sorted([tuple(e) for e in chain.values()] + extra, key=lambda e: e[0] * sign)
    marks = {("dup", number(x)) for x in outcome.get("dup", ())}
    if "inserted" in outcome:
        marks.add(("inserted", outcome["inserted"]))
    return exons, marks


def site(kind, p, canonical, plus):
    """(kind, where) of an event known only by a place: as changes() would name it."""
    sign = 1 if plus else -1
    inside = next((k for k, c in enumerate(canonical, 1) if c[0] <= p <= c[1]), None)
    before = sum(last_base(c, plus) * sign < p * sign for c in canonical)
    if kind == "donor":
        return kind, inside or before
    if kind == "acceptor":
        return kind, inside or before + 1
    return kind, before


def match(outcome, exons, marks, canonical, plus, position, wild=()):
    """(same event, same junctions) of one documented outcome and one predicted
    chain with its marks: {("dup", k): length} and {("inserted", n): lengths}.
    Same junctions is None when the outcome names a site without its junction.
    wild: the changes the WT prediction has already, which don't count."""
    want, want_marks = expected_chain(outcome, canonical, plus, position)
    named = {site(kind, position(p), canonical, plus) for kind, p in outcome.get("site", ())}
    ends = [(last_base, position(outcome["pe_end"]))] if "pe_end" in outcome else []
    ends += [(first_base, position(outcome["pe_start"]))] if "pe_start" in outcome else []
    kinds = changes(want, canonical, plus) | named | {site("pseudoexon", p, canonical, plus) for _f, p in ends}
    got = changes(exons, canonical, plus) - set(wild)
    same = bool(kinds or want_marks) and kinds <= got
    new = junctions(want, plus) - junctions(canonical, plus)
    exact = None if named else same and new <= junctions(exons, plus) and all(
        any(e and f(e, plus) == p for e in exons) for f, p in ends)
    for kind, x in want_marks:
        if kind == "dup":
            size = canonical[x - 1][1] - canonical[x - 1][0] + 1
            same = same and ("dup", x) in marks
            exact = exact is not None and exact and marks.get(("dup", x)) == size
        else:
            sizes = marks.get(("inserted", 0), ())
            same = same and bool(sizes)
            exact = exact is not None and exact and bool(sizes) and (x == 0 or x in sizes)
    return same, exact


def compare(outcomes, predicted, canonical, plus, position, wild=()):
    """Rank 1 and the top ranks against any documented outcome. predicted: [(exons, marks)], rank order."""
    out = {"event_rank1": False, "exact_rank1": False, "event_top": False, "exact_top": False, "outcome": None}
    if all(o.get("site") for o in outcomes):
        out["exact_rank1"] = out["exact_top"] = None
    for i, o in enumerate(outcomes):
        for rank, (exons, marks) in enumerate(predicted, 1):
            same, exact = match(o, exons, marks, canonical, plus, position, wild)
            if same and out["outcome"] is None:
                out["outcome"] = i
            out["event_top"] = out["event_top"] or same
            if exact is not None:
                out["exact_top"] = out["exact_top"] or exact
            if rank == 1:
                out["event_rank1"] = out["event_rank1"] or same
                if exact is not None:
                    out["exact_rank1"] = out["exact_rank1"] or exact
    return out


# --- the cases through the engine ------------------------------------------------------------------------

def read_cases():
    with open(CASES, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


class Engine:
    """What the tab holds, for hg38: the gene table, the chromosome, the scores (loaded when first needed)."""

    def __init__(self):
        helpers.setup()
        from pyfaidx import Fasta

        from spliceai_pipeline import cdna_analysis as ca
        from spliceai_pipeline.mane import find_cached_summary

        self.table = ca.GeneTable(BUILD, str(helpers.THIRDPARTY / "refseq"), find_cached_summary(str(helpers.MANE_DIR)))
        self.fasta = Fasta(str(helpers.HG38))
        self.scores = None

    def gene(self, case):
        hits = self.table.find(case["gene"])
        if len(hits) != 1:
            raise ValueError(f"{case['gene']}: {len(hits)} genes by this name")
        if hits[0].canonical is None:
            raise ValueError(f"{case['gene']}: no curated transcript")
        return hits[0]

    def compat(self, case):
        from spliceai_pipeline import cdna_analysis as ca

        gene = self.gene(case)
        record = self.fasta["chr" + gene.chrom]
        fetch = lambda a, b: record[a - 1:b].seq
        return gene, fetch, ca.from_coordinates(case["events"], gene, BUILD, fetch, len(record))

    def run(self, case):
        from spliceai_pipeline import cdna_analysis as ca, refdb

        gene, fetch, compat = self.compat(case)
        if self.scores is None:
            self.scores = refdb.ReferenceScores(BUILD, str(helpers.REF_DIR / f"{BUILD}.fa"),
                                                reader=refdb.Reader(os.path.join(refdb.default_root(), BUILD)),
                                                refseq_dir=str(helpers.THIRDPARTY / "refseq"),
                                                mane_dir=str(helpers.MANE_DIR))
        return ca.run_analysis(BUILD, gene, compat, self.scores, fetch, self.table)


def position_reader(result):
    from spliceai_pipeline import cdna

    reference = next(t.cdna for t in result.transcripts if t.kind == "Reference")
    return lambda text: int(text[2:]) if text.startswith("g.") else cdna.c_to_genomic(reference, text)


def predicted_chain(result, t):
    """(exons, marks) of sample transcript t for compare()."""
    from spliceai_pipeline import cdna_analysis as ca

    marks = {}
    inserted = []
    for n, (copy, _at, length) in ca.in_insertion(result, t).items():
        if copy:
            marks[("dup", copy)] = length
        else:
            inserted.append(length)
    if inserted:
        marks[("inserted", 0)] = inserted
    return list(t.ref_exons), marks


def truth_signal(case, result, canonical, plus, position):
    """(distance from the variant to the nearest documented site, [what SpliceAI says at each]): the sites the
    outcomes name, each as [kind, position, reference value, sample value]; for the Failures sheet."""
    wt, w0, pv, p0, cmap = result.wt_values, result.plan["start"], result.patient.values, result.patient.start, result.cmap
    sites = set()
    for o in case["outcomes"]:
        chain, _ = expected_chain(o, canonical, plus, position)
        for d, a in junctions(chain, plus) - junctions(canonical, plus):
            sites |= {("donor", d), ("acceptor", a)}
        if "pe_end" in o:
            sites.add(("donor", position(o["pe_end"])))
        if "pe_start" in o:
            sites.add(("acceptor", position(o["pe_start"])))
    out = []
    for kind, x in sorted(sites):
        k = 1 if kind == "donor" else 0
        h = cmap.to_hap(x)
        ref = float(wt[x - w0, k]) if 0 <= x - w0 < len(wt) else 0.0
        sample = float(pv[h - p0, k]) if h is not None and 0 <= h - p0 < len(pv) else 0.0
        out.append([kind, x, round(ref, 3), round(sample, 3)])
    edges = [p for e in result.edits for p in (e.start, e.end)]
    return min((abs(x - p) for _k, x in sites for p in edges), default=None), out


def summarise(case, result):
    from spliceai_pipeline import cdna_analysis as ca, wt_concordance

    gene, c = result.gene, result.gene.canonical
    plus = c.strand == "+"
    reference = [t for t in result.transcripts if t.kind == "Reference"][0]
    wt = next(t for t in result.transcripts if t.kind == "WT predicted")
    sample = [t for t in result.transcripts if t.kind == "Patient"]
    canonical = [tuple(e) for e in reference.exons]
    position = position_reader(result)
    table = [r for r in ca.table_rows(result) if r["Type"] == "Sample"]
    top = [predicted_chain(result, t) for t in sample[:TOP]]
    wild = changes(wt.ref_exons, canonical, plus)
    outcome = compare(case["outcomes"], top, canonical, plus, position, wild)
    local = ca.local_concordance(result)
    known = wt_concordance.for_gene(BUILD, gene.name)
    lead = ca.lead(result)
    rows = []
    for t, r, (exons, marks) in zip(sample[:TOP], table, top):
        rows.append({"rank": t.rank, "relative": t.relative, "change_vs_wt": t.change_vs_wt, "change_vs_ref": t.change_vs_ref,
                     "max_delta": t.max_delta, "evidence": r["Evidence vs WT"], "frame": r["Frame"], "nmd": r["PTC/NMD"],
                     "changes": sorted(changes(exons, canonical, plus)), "marks": [[*k, v] for k, v in marks.items()],
                     "exons": exons})
    distance, signal = truth_signal(case, result, canonical, plus, position)
    return {"id": case["id"], "class": case["class"], "gene": gene.name, "transcript": c.accession, "strand": c.strand,
            "distance": distance, "signal": signal,
            "exons_canonical": len(canonical), "top": rows, "gap": lead[0] if lead else None,
            "wt_exons": wt.ref_exons, "wt_changes": sorted(wild),
            "spliceogenic": sample[0].max_delta is not None,
            "spliceogenic_strong": sample[0].max_delta is not None and sample[0].max_delta >= STRONG,
            "weak": "weak evidence" in table[0]["Evidence vs WT"], "wt_vs_ref": wt.change_vs_ref,
            "gene_identical": known and known["identical"], "gene_exons": known and known["exons"],
            "gene_whole_exact": known and bool(known["whole_exact"]),
            "local": [d["text"] for d in local],
            "local_same": all(not d["differing"] and not d["introns"] for d in local),
            "where": [d["where"] for d in local], "compare": outcome,
            "notices": result.warnings, "seconds": round(result.seconds, 1)}


def check_all():
    engine = Engine()
    bad = []
    for case in read_cases():
        try:
            gene, _fetch, compat = engine.compat(case)
            if not compat.alignment.events:
                raise ValueError("no event")
        except Exception as e:
            bad.append(case["id"])
            print(f"CANNOT READ {case['id']}: {e}", flush=True)
    check("every case is read by the From coordinates check, on the reference FASTA", not bad, f"{len(bad)}: {bad[:20]}")


def run():
    engine = Engine()
    os.makedirs(RESULTS, exist_ok=True)
    only = set(sys.argv[2:])
    cases = [c for c in read_cases() if (not only or c["id"] in only) and not (RESULTS / f"{c['id']}.json").exists()]
    print(f"{len(cases)} to run", flush=True)
    for i, case in enumerate(cases, 1):
        t0 = time.time()
        try:
            out = summarise(case, engine.run(case))
        except Exception as e:
            out = {"id": case["id"], "error": f"{type(e).__name__}: {e}"}
        path = RESULTS / f"{case['id']}.json"
        with open(f"{path}.tmp", "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, default=str)
        os.replace(f"{path}.tmp", path)
        print(f"{i}/{len(cases)} {case['id']}: {time.time() - t0:.0f} s; "
              f"{out.get('error') or out['top'][0]['change_vs_wt'][:80]}", flush=True)


def ds():
    """The VCF tab's run on every case with a VCF record: masked, scored live, its skip reasons kept."""
    from spliceai_pipeline import cli

    helpers.setup()
    rows = [c for c in read_cases() if c.get("vcf")]
    path = OUT / "ds_input.vcf"
    records = [f"{v[0]}\t{v[1]}\t{c['id']}\t{v[2]}\t{v[3]}\t.\t.\t.\n" for c in rows for v in [c["vcf"]]]
    order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    records.sort(key=lambda r: (order[r.split("\t")[0]], int(r.split("\t")[1])))
    path.write_text(helpers.vcf_text(records))
    cli.main([str(path), "--build", BUILD, "--mode", "masked", "--fasta", str(helpers.HG38), "--skip-precomputed",
              "-o", str(OUT / "ds_output.tsv")])


def images(ids):
    os.environ.setdefault("QT_QPA_PLATFORM", "windows" if sys.platform == "win32" else "offscreen")
    from PySide6.QtWidgets import QApplication

    import scenarios
    from spliceai_gui import theme

    engine = Engine()
    app = QApplication(sys.argv)
    theme.apply_theme(app)
    scenarios.IMAGES = IMAGES
    os.makedirs(IMAGES, exist_ok=True)
    by = {c["id"]: c for c in read_cases()}
    for case_id in ids:
        print(case_id, scenarios.export(app, engine.run(by[case_id]), case_id), flush=True)


if __name__ == "__main__":
    command = sys.argv[1] if len(sys.argv) > 1 else ""
    if command == "check":
        check_all()
    elif command == "run":
        run()
    elif command == "ds":
        ds()
    elif command == "images":
        images(sys.argv[2:])
    else:
        print(__doc__)
