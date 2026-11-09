"""Spec 17 item 6: a test corpus of synthetic structural changes, and what the
Sequence to cDNA engine makes of it.

    .venv\\Scripts\\python tests\\scenarios.py generate     the FASTAs and scenarios.tsv, from a fixed seed
    .venv\\Scripts\\python tests\\scenarios.py check        every file through the check: found = truth? (no model)
    .venv\\Scripts\\python tests\\scenarios.py run          the model, the control and the images; resumes
    .venv\\Scripts\\python tests\\scenarios.py report       docs/phase_reports/phase_13_scenarios.md and report.html

The corpus goes to test_scenarios\\, next to the repository and never in it.
Every sample is the reference genome with changes put in here: synthetic, not
patient data. The files are read as the tab reads a paste
(cdna_analysis.check_compatibility), run as it runs them (run_analysis) and
drawn by the viewer's own export. Truth is written as generated and in its
normal form (left-aligned), which is computed here, not by the aligner, so
the comparison means something.
"""
import csv
import json
import os
import random
import sys
import time
from pathlib import Path

import helpers
from helpers import check

OUT = Path(os.environ.get("SPLICEAI_SCENARIOS") or helpers.REPO.parent / "test_scenarios")
FASTA, RESULTS, REPORT, IMAGES = OUT / "fasta", OUT / "results", OUT / "report", OUT / "report" / "images"
MANIFEST = OUT / "scenarios.tsv"
SEED = 13
FLANK = 2000            # bases of the genome kept on each side of the canonical transcript
BUILDS = ("hg38", "hg19")
# Strands, sizes and exon counts differ on purpose; MYBPC3, RWDD4 and CDH23 have exons of 3 nt.
GENES = ("VHL", "TP53", "DMPK", "MYBPC3", "RWDD4", "MLH1", "BRCA2", "PTEN", "SCN1A", "NF1", "CDH23", "DMD")
SMALL_EXONS = ("MYBPC3", "RWDD4", "CDH23")
DEEP = 1000             # an insertion this far from any exon of any transcript
SPECIFICITY_N = 20      # the tab's default

COLUMNS = ["id", "file", "gene", "build", "strand", "transcript", "scenario", "type", "as_generated", "normal_form"]


def complement(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def introns(exons):
    """[(k, lo, hi)]: intron k, after exon k in transcript order (exons as wt_concordance.in_order gives them)."""
    return [(k, min(a[1], b[1]) + 1, max(a[0], b[0]) - 1) for k, (a, b) in enumerate(zip(exons, exons[1:]), 1)]


class Locus:
    """The gene on one build: its exons, and picks that land on the same place on both builds' copies."""

    def __init__(self, gene, record, in_order):
        c = gene.canonical
        self.gene, self.c, self.exons = gene, c, in_order(c)
        self.plus = c.strand == "+"
        self.lo, self.hi = max(1, c.start - FLANK), min(len(record), c.end + FLANK)
        self.chrom = "chr" + gene.chrom
        self.seq = record[self.lo - 1:self.hi].seq.upper()
        self.introns = introns(self.exons)
        self.all_exons = sorted({tuple(e) for t in gene.transcripts for e in t.exons} | {tuple(e) for e in c.exons})

    def ref(self, a, b):
        return self.seq[a - self.lo:b - self.lo + 1]

    def from_exon(self, x):
        return min(0 if a <= x <= b else min(abs(x - a), abs(x - b)) for a, b in self.all_exons)

    def coding_exons(self):
        """Internal exons (numbered in transcript order) fully inside the CDS, away from the ATG and the stop."""
        lo, hi = self.c.cds or (self.c.start, self.c.end)
        return [k for k, (a, b) in enumerate(self.exons, 1) if 1 < k < len(self.exons) and lo < a and b < hi]

    def around(self, k):
        """Lengths of the introns before and after exon k."""
        return (next(h - l + 1 for i, l, h in self.introns if i == k - 1),
                next(h - l + 1 for i, l, h in self.introns if i == k))

    def apply(self, events):
        """The sample: the events put into the sequence, in transcript orientation, as gene exports give it."""
        seq, at = [], self.lo
        for e in sorted(events, key=lambda e: e["start"]):
            seq.append(self.ref(at, e["start"] - 1))
            if e["kind"] == "inversion":
                seq.append(complement(self.ref(e["start"], e["end"])))
            elif e["kind"] != "deletion":
                seq.append(e["inserted"])
            at = e["end"] + 1
        seq.append(self.ref(at, self.hi))
        joined = "".join(seq)
        return joined if self.plus else complement(joined)

    def normal(self, e):
        """Truth in its normal form: indels left-aligned on the + strand; an inversion's ends must differ from it."""
        kind, a, b, ins = e["kind"], e["start"], e["end"], e.get("inserted", "")
        if kind == "deletion":
            while self.ref(a - 1, a - 1) == self.ref(b, b):
                a, b = a - 1, b - 1
        elif kind in ("insertion", "duplication"):
            while self.ref(a - 1, a - 1) == ins[-1]:
                a, ins = a - 1, ins[-1] + ins[:-1]
            b = a - 1
        return {"kind": kind, "start": a, "end": b, "inserted": ins}


def text_of(e):
    size = e["end"] - e["start"] + 1 if e["end"] >= e["start"] else len(e["inserted"])
    if e["end"] < e["start"]:
        return f"{e['kind']} {size} bp before {e['start']}: {e['inserted']}" if size <= 60 else \
            f"{e['kind']} {size} bp before {e['start']}: {e['inserted'][:20]}...{e['inserted'][-20:]}"
    return f"{e['kind']} {size} bp {e['start']}-{e['end']}"


def bases(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def deep_insertion(locus, rng, keep_off=()):
    """(event, intron): 30-50 random bases at least DEEP from any exon of any transcript."""
    room = [i for i in locus.introns if i[2] - i[1] + 1 >= 2 * DEEP + 100]
    n = rng.randint(30, 50)
    for _ in range(5000):
        k, lo, hi = rng.choice(room)
        at = lo + DEEP + int(rng.random() * (hi - lo - 2 * DEEP))
        if locus.from_exon(at) >= DEEP and locus.from_exon(at - 1) >= DEEP \
                and all(abs(at - x) > 500 for x in keep_off):
            return {"kind": "insertion", "start": at, "end": at - 1, "inserted": bases(rng, n)}, k
    raise ValueError("no deep place for an insertion")


def junction_deletion(locus, rng, side, k=None):
    """Across exon k's acceptor or donor: some of the intron and some of the exon, 60-150 bp."""
    exons = [x for x in locus.coding_exons() if locus.exons[x - 1][1] - locus.exons[x - 1][0] + 1 >= 50
             and min(locus.around(x)) >= 100]
    k = k or rng.choice(exons)
    a, b = locus.exons[k - 1]
    n = rng.randint(60, 150)
    in_exon = min(rng.randint(15, n - 15), b - a + 1 - 15)
    in_intron = n - in_exon
    first_side = (side == "acceptor") == locus.plus     # the junction is at the exon's low coordinate
    start, end = (a - in_intron, a + in_exon - 1) if first_side else (b - in_exon + 1, b + in_intron)
    return {"kind": "deletion", "start": start, "end": end}, k


def exon_span(locus, rng, k, least, most):
    """Exon k with some of the intron on each side (at most a third of it): its low and high coordinates."""
    a, b = locus.exons[k - 1]
    before, after = locus.around(k)
    if not locus.plus:
        before, after = after, before
    return a - min(rng.randint(least, most), before // 3), b + min(rng.randint(least, most), after // 3)


def whole_exon(locus, rng, k):
    start, end = exon_span(locus, rng, k, 30, 150)
    return {"kind": "deletion", "start": start, "end": end}


def inversion(locus, a, b):
    """start..end with ends that aren't complementary, or the inversion's outermost bases wouldn't change."""
    while locus.ref(a, a) == complement(locus.ref(b, b)):
        b += 1
    return {"kind": "inversion", "start": a, "end": b}


def intron_inversion(locus, rng, keep_off=()):
    room = [i for i in locus.introns if i[2] - i[1] + 1 >= 700]
    for _ in range(5000):
        k, lo, hi = rng.choice(room)
        n = min(rng.randint(60, 300), (hi - lo) // 3)
        a = lo + 150 + int(rng.random() * (hi - lo - 300 - n))
        if all(x < a - 300 or x > a + n + 300 for x in keep_off):
            return inversion(locus, a, a + n), k
    raise ValueError("no intron for an inversion")


def exon_inversion(locus, rng):
    k = rng.choice([x for x in locus.coding_exons() if min(locus.around(x)) >= 250])
    return inversion(locus, *exon_span(locus, rng, k, 60, 200)), k


def tandem_duplication(locus, rng):
    """Exon k and some of the intron on each side, copied right behind itself."""
    k = rng.choice([x for x in locus.coding_exons() if min(locus.around(x)) >= 120])
    start, end = exon_span(locus, rng, k, 20, 100)
    return {"kind": "duplication", "start": end + 1, "end": end, "inserted": locus.ref(start, end)}, k


def scenarios(locus, name):
    """[(scenario, type, events)]: the same picks on both builds (the rng is seeded without the build)."""
    rng = lambda what: random.Random(f"{SEED}:{name}:{what}")
    size = lambda e: e["end"] - e["start"] + 1 if e["end"] >= e["start"] else len(e["inserted"])
    out = [("NEG_unchanged", "negative control", [])]
    ins, k = deep_insertion(locus, rng("insertion"))
    out.append((f"INS{size(ins)}_deep_intron{k}", "deep intronic insertion", [ins]))
    for side in ("acceptor", "donor"):
        e, k = junction_deletion(locus, rng("junction " + side), side)
        out.append((f"DEL{size(e)}_exon{k}_junction_{side}", f"junction deletion ({side} side)", [e]))
    coding = locus.coding_exons()
    lengths = {k: locus.exons[k - 1][1] - locus.exons[k - 1][0] + 1 for k in coding}
    k = min(coding, key=lambda x: (lengths[x], x)) if name in SMALL_EXONS else rng("whole exon").choice(coding)
    e = whole_exon(locus, rng("whole exon flanks"), k)
    out.append((f"DEL{size(e)}_exon{k}_whole", "whole-exon deletion", [e]))
    e, k = intron_inversion(locus, rng("intron inversion"))
    out.append((f"INV{size(e)}_intron{k}", "inversion in an intron", [e]))
    e, k = exon_inversion(locus, rng("exon inversion"))
    out.append((f"INV{size(e)}_exon{k}", "inversion containing an exon", [e]))
    e, k = tandem_duplication(locus, rng("duplication"))
    out.append((f"DUP{size(e)}_exon{k}_tandem", "tandem exon duplication", [e]))
    r = rng("complex")
    for _ in range(200):            # three changes in different places, none next to another
        k = r.choice(coding)
        deletion = whole_exon(locus, r, k)
        try:
            ins, _i = deep_insertion(locus, r, keep_off=(deletion["start"], deletion["end"]))
            inv, _j = intron_inversion(locus, r, keep_off=(deletion["start"], deletion["end"], ins["start"]))
        except ValueError:
            continue
        out.append(("COMPLEX_del_ins_inv", "complex (three events)", [deletion, ins, inv]))
        return out
    raise ValueError(f"{name}: no room for the complex sample")


def setup():
    helpers.setup()
    from pyfaidx import Fasta

    from spliceai_pipeline import cdna_analysis as ca
    from spliceai_pipeline.mane import find_cached_summary
    from spliceai_pipeline.wt_concordance import in_order

    summary = find_cached_summary(str(helpers.MANE_DIR))
    tables = {b: ca.GeneTable(b, str(helpers.THIRDPARTY / "refseq"), summary) for b in BUILDS}
    fastas = {"hg19": Fasta(str(helpers.HG19)), "hg38": Fasta(str(helpers.HG38))}
    return tables, fastas, in_order


def write_fasta(path, header, seq):
    with open(path, "w", newline="\n") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")


def generate():
    tables, fastas, in_order = setup()
    rows = []
    for name in GENES:
        for build in BUILDS:
            gene = tables[build].find(name)[0]
            locus = Locus(gene, fastas[build]["chr" + gene.chrom], in_order)
            folder = FASTA / name
            os.makedirs(folder, exist_ok=True)
            where = f"{locus.chrom}:{locus.lo}-{locus.hi}"
            orientation = "transcript orientation" + ("" if locus.plus else " (reverse complement of the + strand)")
            write_fasta(folder / f"{name}_{build}_canonical.fa",
                        f"{name} {locus.c.accession} {build} {where} {orientation}, reference", locus.seq if locus.plus
                        else complement(locus.seq))
            for scenario, kind, events in scenarios(locus, name):
                sid = f"{name}_{build}_{scenario}"
                write_fasta(folder / f"{sid}.fa", f"{sid} SYNTHETIC, not patient data; {name} {locus.c.accession} "
                                                  f"{build} {where} {orientation}", locus.apply(events))
                rows.append({"id": sid, "file": f"fasta/{name}/{sid}.fa", "gene": name, "build": build,
                             "strand": gene.strand, "transcript": locus.c.accession, "scenario": scenario, "type": kind,
                             "as_generated": json.dumps(events), "normal_form": json.dumps(
                                 sorted((locus.normal(e) for e in events), key=lambda e: e["start"]))})
            print(name, build, f"{len(locus.exons)} exons, {locus.hi - locus.lo + 1:,} bp", flush=True)
    with open(MANIFEST, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} scenarios in {OUT}")


def read_manifest():
    with open(MANIFEST, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["normal_form"] = json.loads(r["normal_form"])
    return rows


def found_events(compat):
    return [{"kind": e.kind, "start": e.start, "end": e.end, "inserted": "" if e.kind == "inversion" else e.inserted}
            for e in sorted(compat.alignment.events, key=lambda e: e.start)]


def checked(row, tables, fastas):
    """The file through the tab's check: (gene, fetch, compat, what was found, equal to the truth, why not)."""
    from spliceai_pipeline import cdna_analysis as ca, haplotype

    gene = tables[row["build"]].find(row["gene"])[0]
    record = fastas[row["build"]]["chr" + gene.chrom]
    fetch = lambda a, b: record[a - 1:b].seq
    text = open(OUT / row["file"], encoding="ascii").read()
    compat = ca.check_compatibility(text, gene, fetch, len(record))
    found = found_events(compat)
    problems = [f"{level}: {text}" for level, text, _ in ca.check_notices(compat) if level == ca.CRITICAL]
    why = []
    if found != row["normal_form"]:
        why.append("found " + "; ".join(text_of(e) for e in found) + " | truth " + "; ".join(
            text_of(e) for e in row["normal_form"]))
    if not all(e.apply for e in compat.alignment.events):
        why.append("not every difference is applied by default")
    hap, _ = haplotype.apply_events(compat.ref_seq, compat.region[0], compat.alignment.events)
    if (haplotype.reverse_complement(compat.sequence) if compat.alignment.reverse else compat.sequence) not in hap:
        why.append("the events applied don't give the file's sequence")
    return gene, fetch, compat, found, not why, why, problems


def check_all():
    tables, fastas, _ = setup()
    rows = read_manifest()
    t0, misses = time.time(), []
    for row in rows:
        _gene, _fetch, compat, _found, ok, why, problems = checked(row, tables, fastas)
        if not ok or problems or compat.alignment.reverse != (row["strand"] == "-"):
            misses.append(row["id"])
            print(f"MISS {row['id']}: {why} {problems} reverse={compat.alignment.reverse}", flush=True)
    check(f"every one of the {len(rows)} files: the differences found are the truth, exactly", not misses,
          f"{len(misses)} differ: {misses[:12]}")
    print(f"{time.time() - t0:.0f} s")


def export(app, result, sid):
    """The viewer's own exports at 2x: the schematic, and the stacked view centred on each event."""
    from PIL import Image
    from PySide6.QtCore import Qt

    from spliceai_gui.viewer_window import ViewerWindow

    v = ViewerWindow(result)
    v.setAttribute(Qt.WA_DontShowOnScreen)
    v.resize(1500, 900)
    v.show()
    app.processEvents()
    names = [f"{sid}_schematic.png"]
    v.exons.render_png(str(IMAGES / names[0]), 1100)
    v.tabs.setCurrentIndex(1)
    app.processEvents()
    s = v.stacked
    spans = [(a, b) for a, b, _ in sorted(s.data.events)]
    if not spans:                       # the control: the middle exon
        exon = result.gene.canonical.exons[len(result.gene.canonical.exons) // 2]
        spans = [(s.data.column_of_ref(exon[0]), s.data.column_of_ref(exon[1]))]
    for i, (a, b) in enumerate(spans, 1):
        pad = max(150, (b - a) // 2)
        s.show_genomic(a - pad, b + pad)
        names.append(f"{sid}_stacked{i}.png")
        s.render_png(str(IMAGES / names[-1]))
    v.close()
    v.deleteLater()
    app.processEvents()
    for name in names:                  # far fewer than 256 colours are drawn; a third of the size
        Image.open(IMAGES / name).convert("RGB").quantize(256, dither=Image.Dither.NONE).save(IMAGES / name, optimize=True)
    return names


def summary_of(row, result, found, ok, why, problems, controls, images, seconds):
    from spliceai_pipeline import cdna_analysis as ca

    by = {(t.kind, t.rank): t for t in result.transcripts}
    sample, wt = by[("Patient", 1)], by[("WT predicted", 1)]
    table = next(r for r in ca.table_rows(result) if r["Type"] == "Sample" and r["Rank"] == 1)
    local = ca.local_concordance(result)
    lead = ca.lead(result)
    events = ca.event_summary(result)
    return {"id": row["id"], "gene": row["gene"], "build": row["build"], "strand": row["strand"],
            "transcript": row["transcript"], "scenario": row["scenario"], "type": row["type"], "file": row["file"],
            "truth": row["normal_form"], "found": found, "exact": ok, "why": why, "critical": problems,
            "reverse": result.compat.alignment.reverse, "labels": [d["label"] for d in events],
            "origins": [d["origin"] for d in events], "effects": [d["effect"] for d in events],
            "wt_near": [d["text"] for d in local], "wt_near_same": all(not d["differing"] and not d["introns"] for d in local),
            "wt_vs_ref": wt.change_vs_ref, "change_vs_wt": sample.change_vs_wt, "change_vs_ref": sample.change_vs_ref,
            "relative": sample.relative, "gap": lead[0] if lead else None, "max_delta": sample.max_delta,
            "evidence": table["Evidence vs WT"], "weak": "weak evidence" in table["Evidence vs WT"],
            "frame": table["Frame"], "nmd": table["PTC/NMD"], "protein": table["Protein length"],
            "rna": table["RNA check (RT-PCR)"], "exons": table["Exons"], "reference_exons": len(result.gene.canonical.exons),
            "notices": result.warnings, "specificity": [{"text": c.text(), "observed_changed": c.observed_changed,
                                                         "observed_delta": c.observed_delta, "changed": c.changed,
                                                         "n": c.n, "deltas": c.deltas, "note": c.note} for c in controls],
            "images": images, "seconds": seconds}


def run():
    os.environ.setdefault("QT_QPA_PLATFORM", "windows" if sys.platform == "win32" else "offscreen")
    tables, fastas, _ = setup()
    from PySide6.QtWidgets import QApplication

    from spliceai_gui import theme
    from spliceai_pipeline import cdna_analysis as ca, refdb, specificity

    app = QApplication(sys.argv)
    theme.apply_theme(app)
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(IMAGES, exist_ok=True)
    genes = sys.argv[2:] or GENES
    rows = [r for r in read_manifest() if r["gene"] in genes and not (RESULTS / f"{r['id']}.json").exists()]
    scores = {}
    print(f"{len(rows)} to run", flush=True)
    for i, row in enumerate(rows, 1):
        t0 = time.time()
        build = row["build"]
        gene, fetch, compat, found, ok, why, problems = checked(row, tables, fastas)
        if build not in scores:
            reader = refdb.Reader(os.path.join(refdb.default_root(), build))
            scores[build] = refdb.ReferenceScores(build, str(helpers.REF_DIR / f"{build}.fa"), reader=reader,
                                                  refseq_dir=str(helpers.THIRDPARTY / "refseq"),
                                                  mane_dir=str(helpers.MANE_DIR))
        result = ca.run_analysis(build, gene, compat, scores[build], fetch, tables[build])
        ran = time.time() - t0
        controls = []
        if row["type"] == "deep intronic insertion":
            controls = result.specificity = specificity.run(result, scores[build], SPECIFICITY_N)
        images = export(app, result, row["id"])
        out = summary_of(row, result, found, ok, why, problems, controls, images, round(time.time() - t0, 1))
        path = RESULTS / f"{row['id']}.json"
        with open(f"{path}.tmp", "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1, default=str)
        os.replace(f"{path}.tmp", path)
        print(f"{i}/{len(rows)} {row['id']}: run {ran:.0f} s, all {time.time() - t0:.0f} s; "
              f"{out['change_vs_wt'][:70]}", flush=True)


# --- the report ------------------------------------------------------------------------------------------

BANNER = ("Synthetic data, not patient data: every sample here is the public reference genome with changes put in by "
          "a script. Research use only: these are predictions from SpliceAI, not measurements; RNA (cDNA) confirmation "
          "is required for any clinical conclusion.")
LIMITS_FILE = Path(__file__).with_name("scenarios_limits.md")
HEAD = ["Gene", "Build", "Scenario", "True event", "Detected event", "Exact match", "WT vs reference near event",
        "Rank-1 change", "Relative score", "Gap to rank 2", "Max delta", "Weak evidence", "Frame", "NMD", "RNA check"]


def short_event(e):
    size = e["end"] - e["start"] + 1 if e["end"] >= e["start"] else len(e["inserted"])
    word = {"insertion": "INS", "deletion": "DEL", "inversion": "INV", "duplication": "DUP"}.get(e["kind"], e["kind"])
    place = f"before {e['start']:,}" if e["end"] < e["start"] else f"{e['start']:,}-{e['end']:,}"
    return f"{word} {size:,} bp {place}"


def near_text(r):
    if not r["wt_near"]:
        return r["wt_vs_ref"]
    if r["wt_near_same"]:
        return "same"
    prefix = "the WT prediction differs from the reference for "
    return "differs: " + "; ".join(t.split("): ", 1)[-1].replace(prefix, "").rstrip(".") for t in r["wt_near"]
                                   if "differs" in t)


def number(x, digits=3):
    return "" if x is None else f"{x:.{digits}f}"


def cells(r):
    return [r["gene"], r["build"], r["scenario"], "; ".join(short_event(e) for e in r["truth"]) or "none",
            "; ".join(short_event(e) for e in r["found"]) or "none", "yes" if r["exact"] else "NO", near_text(r),
            r["change_vs_wt"], number(r["relative"]), number(r["gap"]), number(r["max_delta"], 2),
            "yes" if r["weak"] else "no", r["frame"], r["nmd"], r["rna"]]


def what_was_done(r):
    gene, strand = r["gene"], "+" if r["strand"] == "+" else "−"
    start = f"This file is the reference sequence of {gene} ({r['transcript']}, {r['build']}, {strand} strand)"
    place = r["labels"][0].split(" · ", 2)[-1] if r["labels"] else ""
    return {"negative control": f"{start}, unchanged: a control, in which the tool must find nothing.",
            "deep intronic insertion": f"{start} with a stretch of random bases put deep inside an intron, far from "
                                       f"any exon: {place}.",
            "junction deletion (acceptor side)": f"{start} with a deletion across the start of an exon, taking the "
                                                 f"splice acceptor with it: {place}.",
            "junction deletion (donor side)": f"{start} with a deletion across the end of an exon, taking the splice "
                                              f"donor with it: {place}.",
            "whole-exon deletion": f"{start} with one whole exon and some of the introns around it deleted: {place}.",
            "inversion in an intron": f"{start} with a stretch of an intron turned around, all exons untouched: {place}.",
            "inversion containing an exon": f"{start} with one exon and some of the introns around it turned around: "
                                            f"{place}.",
            "tandem exon duplication": f"{start} with one exon and some of the introns around it copied right behind "
                                       f"itself: {place}.",
            "complex (three events)": f"{start} with three separate changes in one sample: a deletion of a whole exon, "
                                      "a random insertion deep in an intron and an inversion in an intron."}[r["type"]]


def narrative(r):
    out = [what_was_done(r)]
    n = len(r["truth"])
    if r["exact"]:
        out.append("The tool found nothing, as it should." if not n else
                   "The tool found exactly that: the same type, the same place and the same bases."
                   if n == 1 else f"The tool found all {n} changes exactly: the same types, places and bases.")
    else:
        out.append("The tool did not find exactly what was made: " + "; ".join(r["why"]) + ".")
    if r["critical"]:
        out.append("It flagged a critical problem: " + " ".join(r["critical"]))
    if not r["wt_near"]:
        out.append(f"On this unchanged sequence SpliceAI's top transcript compared with the annotated one: "
                   f"{r['wt_vs_ref'].replace(' | ', '; ')}.")
    elif r["wt_near_same"]:
        out.append("On the unchanged reference, SpliceAI's top transcript agrees with the annotated exons around the "
                   "change, so the prediction there rests on firm ground.")
    else:
        out.append("On the unchanged reference, SpliceAI's top transcript already differs from the annotated one "
                   "around the change, so the prediction there needs caution: " + near_text(r).replace("differs: ", "")
                   + ".")
    if r["change_vs_wt"] == "No change":
        out.append("With the change in place the top transcript stays the same as the prediction for the unchanged "
                   "sequence: the model predicts no effect on splicing here.")
    else:
        out.append(f"With the change in place the top transcript differs from the prediction for the unchanged "
                   f"sequence: {r['change_vs_wt'].replace(' | ', '; ')}.")
    lead = f" and leads the second candidate by {r['gap']:.3f}" if r["gap"] is not None else ""
    out.append(f"It has a relative score of {r['relative']:.3f}{lead}" + (
        f"; the largest change in a splice-site score behind it is {r['max_delta']:.2f}" if r["max_delta"] is not None
        else "") + (", which counts as weak evidence." if r["weak"] else "."))
    if r["nmd"]:
        same = ", the same as for the unchanged sequence" if r["change_vs_wt"] == "No change" else ""
        out.append(f"Reading frame and decay{same}: {r['nmd'].rstrip('.')}.")
    if r["rna"]:
        out.append(f"To check this on RNA: {r['rna']}.")
    for c in r["specificity"]:
        if c["n"]:
            median = sorted(c["deltas"])[len(c["deltas"]) // 2]
            out.append(f"Sequence-specificity control: {c['changed']} of {c['n']} random sequences of the same length "
                       f"and composition changed the top transcript (the real one "
                       f"{'did' if c['observed_changed'] else 'did not'}); the largest score change was "
                       f"{c['observed_delta']:.2f} for the real sequence, median {median:.2f} for the random ones.")
    return " ".join(out)


def by_type(rows):
    kinds = list(dict.fromkeys(r["type"] for r in rows))
    out = []
    for kind in kinds:
        some = [r for r in rows if r["type"] == kind]
        out.append([kind, len(some), sum(r["exact"] for r in some), sum(r["wt_near_same"] if r["wt_near"] else r["wt_vs_ref"] == "No change" for r in some),
                    sum(r["change_vs_wt"] != "No change" for r in some), sum(r["weak"] for r in some)])
    return out


TYPE_HEAD = ["Scenario type", "Files", "Found exactly", "WT matches the reference near the event (control: whole transcript)",
             "Rank-1 differs from the WT prediction", "Weak evidence"]


def md_table(head, rows):
    clean = lambda c: str(c).replace("|", "/").replace("\n", " ")
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
                     + ["| " + " | ".join(clean(c) for c in row) + " |" for row in rows])


def html_table(head, rows):
    import html

    return ("<table><thead><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in head) + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>" for row in rows)
            + "</tbody></table>")


def report():
    import base64
    import html
    import re

    from spliceai_pipeline import splice_paths

    rows = [json.load(open(RESULTS / f"{m['id']}.json", encoding="utf-8")) for m in read_manifest()]
    limits = LIMITS_FILE.read_text(encoding="utf-8")
    stats = by_type(rows)
    exact = sum(r["exact"] for r in rows)
    intro = (f"{len(rows)} synthetic files, {len(GENES)} genes on both builds: {', '.join(GENES)}. Generated by "
             f"`tests/scenarios.py` with the fixed seed {SEED}; the exact truth of each file is in `scenarios.tsv`. "
             f"Each file went through the same code as the Sequence to cDNA tab with the default settings (minimum "
             f"internal exon {splice_paths.MIN_EXON} nt); images are the viewer's own export at 2x. {exact} of {len(rows)} files: the "
             f"differences found equal the truth exactly.")
    md = ["# Phase 13: test scenarios", "", f"**{BANNER}**", "", intro, "",
          "The full report with its images is `test_scenarios\\report\\report.html`, next to the repository; the "
          "images below are in the same folder and show when this file is read on this computer.", "",
          "## Summary by type", "", md_table(TYPE_HEAD, stats), "", limits, "", "## Overview", "",
          md_table(HEAD, [cells(r) for r in rows]), "", "## Scenario by scenario", ""]
    body = []
    for name in GENES:
        some = [r for r in rows if r["gene"] == name]
        body.append(f"<details id='{name}'><summary>{name}: {some[0]['reference_exons']} exons, "
                    f"{'+' if some[0]['strand'] == '+' else '−'} strand, {sum(r['exact'] for r in some)} of {len(some)} "
                    f"found exactly</summary>")
        for r in some:
            md += [f"### {r['id']}", "", narrative(r), "", md_table(HEAD, [cells(r)]), ""]
            md += [f"![{image}](../../../test_scenarios/report/images/{image})" for image in r["images"]] + [""]
            body.append(f"<h3>{html.escape(r['id'])}</h3><p>{html.escape(narrative(r))}</p>"
                        f"<div class='wide'>{html_table(HEAD, [cells(r)])}</div>")
            for image in r["images"]:
                data = base64.b64encode((IMAGES / image).read_bytes()).decode("ascii")
                body.append(f"<img loading='lazy' alt='{html.escape(image)}' src='data:image/png;base64,{data}'>")
        body.append("</details>")
    (helpers.REPO / "docs" / "phase_reports" / "phase_13_scenarios.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    rendered = "".join(f"<h3>{html.escape(x[4:])}</h3>" if x.startswith("### ") else
                       f"<h2>{html.escape(x[3:])}</h2>" if x.startswith("## ") else
                       "<ul>" + "".join("<li>" + html.escape(i[2:]) + "</li>" for i in x.splitlines()) + "</ul>"
                       if x.startswith("- ") else f"<p>{html.escape(x)}</p>"
                       for x in re.split(r"\n\n+", limits.strip()))
    rendered = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", rendered)
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>SpliceAI for Windows: test scenarios</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;max-width:1240px;margin:2em auto;padding:0 1em;color:#222;line-height:1.45}}
.banner{{background:#fff3cd;border:1px solid #d9b44a;padding:.8em 1em;font-weight:600}}
table{{border-collapse:collapse;font-size:12px;margin:.6em 0}}th,td{{border:1px solid #ccc;padding:3px 6px;text-align:left;vertical-align:top}}
th{{background:#eef1f4;position:sticky;top:0}}.wide{{overflow-x:auto}}img{{max-width:100%;display:block;margin:.5em 0;border:1px solid #e3e3e3}}
summary{{font-size:1.15em;font-weight:600;cursor:pointer;padding:.3em 0}}h3{{margin-bottom:.2em}}.over{{overflow:auto;max-height:80vh}}</style></head><body>
<h1>Test scenarios: what the Sequence to cDNA engine makes of known changes</h1><p class="banner">{html.escape(BANNER)}</p>
<p>{html.escape(intro.replace("`", ""))}</p><h2>Summary by type</h2>{html_table(TYPE_HEAD, stats)}{rendered}
<h2>Scenario by scenario</h2><p>Open a gene to see its scenarios, each with a plain description, its row and its images.</p>{"".join(body)}
<h2>Overview</h2><div class="over">{html_table(HEAD, [cells(r) for r in rows])}</div><p class="banner">{html.escape(BANNER)}</p></body></html>"""
    os.makedirs(REPORT, exist_ok=True)
    (REPORT / "report.html").write_text(page, encoding="utf-8")
    print(f"{len(rows)} scenarios, {exact} exact; report.html {os.path.getsize(REPORT / 'report.html') / 1e6:.0f} MB")


if __name__ == "__main__":
    {"generate": generate, "check": check_all, "run": run, "report": report}[sys.argv[1]]()
    if sys.argv[1] == "check":
        helpers.finish()
