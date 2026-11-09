"""Independent review (phase 14): off-by-one, strand and 0-based/1-based problems.
Each section below is one problem found in review; every check here fails until
that problem is fixed, for the reason written in its section. All inputs are
synthetic: no FASTA, no models, no database.

    .venv\\Scripts\\python -W ignore tests\\test_review_phase14.py
"""
import os
import random
from types import SimpleNamespace as NS

import helpers
from helpers import check

helpers.watchdog(300)
helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"

import numpy as np  # noqa: E402

from spliceai_pipeline import align, cdna, cdna_analysis as ca, haplotype, transcripts, vcf_bridge  # noqa: E402
from spliceai_pipeline.align import reverse_complement as rc  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.view_data import ViewData  # noqa: E402


# --- synthetic transcripts and analyses ------------------------------------------------------

def laid_out(tx, cut, strand, cds_first, cds_last):
    """transcript string tx as two exons (split at index cut) around a 30-base intron, on either strand, from
    genomic position 101: (fetch, Transcript, genomic position of transcript index i)."""
    intron = "GT" + "A" * 26 + "AG"
    plus = tx[:cut] + intron + tx[cut:]
    seq = plus if strand == "+" else rc(plus)
    n, start = len(seq), 101

    def geno(i):
        j = i if i < cut else i + len(intron)
        return start + j if strand == "+" else start + n - 1 - j

    exons = sorted([tuple(sorted((geno(0), geno(cut - 1)))), tuple(sorted((geno(cut), geno(len(tx) - 1))))])
    t = transcripts.Transcript("NM_T.1", "T", "chr1", strand, start, start + n - 1,
                               tuple(sorted((geno(cds_first), geno(cds_last)))), tuple(exons))
    return (lambda a, b: seq[a - start:b - start + 1]), t, geno


R0 = 1001


def make_gene(strand, exons_ascending, cds=None):
    tx = transcripts.Transcript("NM_T.1", "T", "chr1", strand, exons_ascending[0][0], exons_ascending[-1][1], cds,
                                tuple(exons_ascending))
    return transcripts.Gene("T", "1", strand, tx.start, tx.end, tuple(exons_ascending), tx, transcripts.RULE_MANE, [tx])


def analysis(ref, gene, events, patient_exons=None):
    """A finished analysis with no model: all scores 0, the transcripts given (default the reference's exons)."""
    fetch = lambda a, b: ref[a - R0:b - R0 + 1]
    region = (R0, R0 + len(ref) - 1)
    a = align.Alignment(len(ref), False, region[0], region[1], 1, len(ref), 10, 10, events=events, method="test")
    compat = ca.Compatibility(ref, region, region, ref, a, {"lines": [], "problems": [], "warnings": []},
                              ca.anchors(gene))
    edits = [Edit.from_event(e) for e in events if e.apply]
    hap, cmap = haplotype.build(ref, R0, edits)
    m = gene.mask
    plan = {"start": m[0] - 50, "end": m[1] + 50, "mask": m, "strand": gene.strand}
    start, end, mask = haplotype.patient_region(cmap, plan)
    pred = haplotype.PatientPrediction(start, end, np.zeros((end - start + 1, 2), np.float32), mask)
    wt = np.zeros((plan["end"] - plan["start"] + 1, 2), np.float32)
    result = ca.Analysis("hg38", gene, ca.Settings(), compat, edits, hap, cmap, plan, wt, "database", pred)
    c = gene.canonical
    exons = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
    lifted = [tuple(sorted((cmap.lift(x, 1), cmap.lift(y, -1)))) for x, y in exons]
    path = lambda e: NS(exons=e, relative=1.0, score=0.0)
    ca.describe(result, fetch, [path(exons)], [path(patient_exons or lifted)])
    return result


rng = random.Random(14)
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
EVENT = lambda kind, start, end, inserted: align.Event(kind, start, "", "", start, end, inserted,
                                                       max(end - start + 1, len(inserted)), apply=True)

# --- 1. The reference's c. numbering comes from translation, not from the annotated CDS (cdna.reference) ----------
# The annotated CDS defines c.1 and the last c. position. With a non-ATG start (MYC, WT1, VEGFA, ...) c_position()
# falls back to n. numbering, still shown as "c."; with an in-frame TGA before the annotated stop (selenoproteins)
# everything past it becomes c.*N.
utr5, utr3 = "GCGCGATCGC", "GGCTTGCAAC"
for strand in "+-":
    cds = "CTG" + "GCA" * 9 + "TAA"
    fetch, t, geno = laid_out(utr5 + cds + utr3, 20, strand, 10, 10 + len(cds) - 1)
    ref = cdna.reference(fetch, t)
    got = [cdna.genomic_to_c(ref, geno(i)) for i in (9, 10, 10 + len(cds) - 1, 10 + len(cds))]
    check(f"1. CTG-start CDS, {strand} strand: c. numbering follows the annotated CDS (c.-1, c.1, last base, c.*1)",
          got == ["-1", "1", str(len(cds)), "*1"], str(got))
    cds = "ATG" + "GCA" * 3 + "TGA" + "GCA" * 3 + "TAA"
    fetch, t, geno = laid_out(utr5 + cds + utr3, 20, strand, 10, 10 + len(cds) - 1)
    ref = cdna.reference(fetch, t)
    got = [cdna.genomic_to_c(ref, geno(i)) for i in (10 + len(cds) - 1, 10 + len(cds))]
    check(f"1. in-frame TGA before the annotated stop, {strand} strand: the annotated stop is the last c. position",
          got == [str(len(cds)), "*1"], str(got))

# --- 2. Heat strips skip the middle pyramid block of a pixel (stacked_view.TrackPainter._scores) --------------------
# At cols_per_px between about 1.5 and 2 times a power of two a pixel spans three pyramid blocks and only the first
# and last are read ("two pyramid entries cover a pixel at most" is false): a peak in the middle block isn't drawn.
from PySide6.QtCore import QRectF  # noqa: E402
from PySide6.QtGui import QColor, QImage, QPainter  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

app = QApplication([])

from spliceai_gui import stacked_view, theme  # noqa: E402
from spliceai_gui.exon_view import ExonView  # noqa: E402


def drawn(cols_per_px, column, n=2000, width=200):
    """A single peak of 1.0 at `column`: is it visible in the heat strip?"""
    d = ViewData()
    d.n, d.threshold = n, 0.05
    d.tracks = {"ref_acc": np.zeros(n, np.float32)}
    d.tracks["ref_acc"][column] = 1.0
    image = QImage(width, 8, QImage.Format_ARGB32_Premultiplied)
    image.fill(QColor("white"))
    p = QPainter(image)
    stacked_view.TrackPainter(d, [], 0.0, cols_per_px, labels=False)._scores(p, QRectF(0, 0, width, 8), "ref_acc",
                                                                             theme.ACCEPTOR, 0, n - 1)
    p.end()
    return any(image.pixelColor(x, 4) != QColor("white") for x in range(width))


for cols_per_px in (3.9, 7.9, 31.9):
    missing = [c for c in range(600) if not drawn(cols_per_px, c)]
    check(f"2. heat strip at {cols_per_px} columns per pixel: every peak is visible", not missing,
          f"{len(missing)} of 600 peaks not drawn, first {missing[:6]}")

# --- 3. The ATG and stop anchors are single bases (cdna_analysis.anchors, align.validate_events) -------------------
# Spec 6.2 and 9: a warning for an event touching the ATG or the stop codon. With the +-1 of validate_events, the G of
# the ATG and the first base of the stop codon are out of reach.
for strand, cds, exons in (("+", (1100, 1399), [(1050, 1500)]), ("-", (1100, 1399), [(1050, 1500)])):
    gene = make_gene(strand, exons, cds)
    a = ca.anchors(gene)
    step = 1 if strand == "+" else -1
    g_of_atg, first_of_stop = a["ATG"] + 2 * step, a["stop codon"] - 2 * step
    for name, pos in (("G of the ATG", g_of_atg), ("first base of the stop codon", first_of_stop)):
        compat = NS(alignment=NS(events=[EVENT("SNV", pos, pos, "A")]), region=(1000, 1600), anchors=a)
        got = ca.chosen_check(compat)["warnings"]
        check(f"3. {strand} strand: an SNV at the {name} is warned about", bool(got), str(got))

# --- 4. c.0 doesn't exist (cdna.c_to_genomic) ---------------------------------------------------------------------
# The viewer's "Go to" takes c.0 for c.-1, c.-0 for c.1 and c.*0 for the stop codon's last base.
fetch, t, _ = laid_out(utr5 + "ATG" + "GCA" * 9 + "TAA" + utr3, 20, "+", 10, 42)
ref = cdna.reference(fetch, t)
for text in ("c.0", "c.-0", "c.*0"):
    try:
        got = f"gave {cdna.c_to_genomic(ref, text)}"
    except ValueError:
        got = ""
    check(f"4. {text} is refused", not got, got)

# --- 5. The frame is taken from the length of everything after the ATG (cdna.summarize, cdna.hgvs_p) -----------------
# A 1-base CDS deletion stays a frameshift whatever happens in the 3'UTR, but -1 there plus -2 in the UTR is a
# multiple of 3, and the transcript reads "in frame" with a delins protein change.
r5 = random.Random(9)
codons = [c for c in (a + b + e for a in "ACGT" for b in "ACGT" for e in "ACGT") if c not in cdna.STOPS]
body = "".join(r5.choice(codons) for _ in range(99))
tail = "".join(r5.choice(codons) for _ in range(40))
seq = "CCCC" + "ATG" + body + "TAA" + tail
one = seq[:100] + seq[101:]
both = one[:330] + one[332:]
as_cdna = lambda s: cdna.translate(cdna.CDNA([(1, len(s))], "+", s, [0]), 4)
r, c = as_cdna(seq), as_cdna(both)
check("5. a CDS frameshift with a compensating 3'UTR change is still a frameshift",
      cdna.summarize(c, r).frame == "frameshift" and "fs" in cdna.hgvs_p(r, c),
      f"{cdna.summarize(c, r).frame}, {cdna.hgvs_p(r, c)}")

# --- 6. A lost stop is an extension (cdna.hgvs_p) ---------------------------------------------------------------
# HGVS: p.(Ter4GlnextTer...). The code gives p.(Ala3_Ter4ins...) and p.(Ter4Asnfs...).
tail = "CAGCGTACTTAGGG"
as_cdna = lambda s: cdna.translate(cdna.CDNA([(1, len(s))], "+", s, [0]), 0)
r = as_cdna("ATGGCTGCATAA" + tail)
for name, s in (("in frame", "ATGGCTGCACAA" + tail), ("by a 1-base deletion", "ATGGCTGCAAA" + tail)):
    got = cdna.hgvs_p(r, as_cdna(s))
    check(f"6. stop lost {name}: described as an extension (Ter4...ext)", got.startswith("p.(Ter4") and "ext" in got, got)

# --- 7. engine_scores: a delins with one side of 1 base (vcf_bridge.engine_scores) ---------------------------------
# get_delta_scores (ref_len > 1, alt_len == 1, or ref_len == 1, alt_len > 1) scores the alt bases at pos; to_hap(pos)
# is None for a replaced base, so the engine gives 0 there, a loss that isn't real.
for r_, alt in (("AG", "C"), ("A", "GC")):
    bases = list(rand(300))
    bases[149:149 + len(r_)] = list(r_)
    region = "".join(bases)
    hap, cmap = haplotype.build(region, 1, [Edit(150, 149 + len(r_), alt, "complex")])
    values = np.zeros((len(hap), 2), np.float32)
    values[149] = (0.7, 0.0)
    wt = np.zeros((len(region), 2), np.float32)
    wt[149] = (0.7, 0.0)
    pred = haplotype.PatientPrediction(1, len(hap), values, (1, len(hap)))
    result = NS(plan={"start": 1, "end": len(region)}, wt_values=wt, patient=pred, cmap=cmap, gene=NS(exons=()))
    ds = vcf_bridge.engine_scores(result, 150, r_, alt, exons=[(1, 2)])[0]
    check(f"7. {r_}>{alt}: no acceptor loss at pos when the alt base scores what the reference did, as "
          "get_delta_scores has it", abs(ds[1]) < 1e-6, f"DS_AL {ds[1]:.3f}")

# --- 8. An exon block is one run of columns even when an inversion splits it (view_data.ViewData._span) --------------
# A sample exon across an inversion's edge lies on two separate runs of columns; the one block drawn covers the
# columns between them too, which aren't in the exon.
ref_seq = rand(3000)
gene = make_gene("+", [(1301, 1400), (1701, 1800), (2101, 2250)])
start, end = 1751, 1900
inverted = rc(ref_seq[start - R0:end - R0 + 1])
result = analysis(ref_seq, gene, [EVENT("inversion", start, end, inverted)])
view = ViewData.from_analysis(result)
row = next(r for r in view.rows if r.group == "patient")
outside = [(b.number, c) for b in row.blocks for c in range(b.first, b.last + 1) if row.membership[c] != b.number]
check("8. no exon block covers columns outside that exon", not outside,
      f"{len(outside)} columns, exons {sorted({n for n, _ in outside})}")

# --- 9. An insertion right after a deletion is marked at its own columns (view_data.ViewData.inserted_after) ---------
# haplotype.columns() puts the inserted columns before the deleted ones; inserted_after(last deleted base) looks after
# that base, finds none and returns the deleted column, so the band and "Zoom to event" are off by the deletion's size.
ref_seq = rand(3000)
gene = make_gene("+", [(1301, 1400), (1701, 1800), (2101, 2250)])
events = [EVENT("deletion", 1500, 1600, ""), EVENT("insertion", 1601, 1600, "ACGTAC")]
view = ViewData.from_analysis(analysis(ref_seq, gene, events))
inserted = np.nonzero(view.kind == haplotype.INSERTED)[0]
got = view.events[1][:2]
check("9. the insertion's event columns are its inserted columns", got == (int(inserted[0]), int(inserted[-1])),
      f"{got}, inserted columns {int(inserted[0])}..{int(inserted[-1])}")

# --- 10. "first" and "last" of a reversed paste are the wrong way round (align.align) --------------------------------
# The notes count along the reverse complement but say "of the sequence": the user's first bases are called the last.
ref_seq = rand(3000)
fetch = lambda c, a, b: ref_seq[a - 1:b]
paste = rand(200) + rc(ref_seq[500:1500])       # the unmatched bases are the first ones of what was pasted
aligned = align.align(paste, ref_seq, 1, fetch, "chr1")
note = next((n for n in aligned.notes if "left out" in n), "")
check("10. a reversed paste: the note names the end of the sequence as pasted",
      aligned.reverse and ("first" in note or "reverse" in note), note)

# --- 11. An exon box takes in the room kept right after it (exon_view.ExonView._row, Axis.x) ---------------------------
# Axis.x(p) at a position with room reads the far side of the room, right for an exon starting there, wrong for an
# exon ending there: the exon is drawn as wide as itself plus the inserted exon, over the inserted exon's own box.
for strand in "+-":
    ref_seq = rand(3000)
    gene = make_gene(strand, [(1301, 1400), (1701, 1800), (2101, 2250)])
    a, b = 1701, 1800
    copy = ref_seq[a - R0:b - R0 + 1]
    events = [EVENT("duplication", b + 1, b, copy)]
    _, cmap = haplotype.build(ref_seq, R0, [Edit.from_event(events[0])])
    inserted = next(s for s in cmap.segments if s.orientation == 0)
    exons = [(1301, 1400), (1701, 1800), (2101, 2250)] if strand == "+" else [(2101, 2250), (1701, 1800), (1301, 1400)]
    lifted = [tuple(sorted((cmap.to_hap(x), cmap.to_hap(y)))) for x, y in exons]
    lifted.insert(lifted.index((cmap.to_hap(a), cmap.to_hap(b))) + (1 if strand == "+" else 0),
                  (inserted.hap_start, inserted.hap_end))
    view = ViewData.from_analysis(analysis(ref_seq, gene, events, lifted))
    widget = ExonView(view)
    widget.resize(1000, 300)
    image = QImage(1000, 300, QImage.Format_ARGB32)
    p = QPainter(image)
    widget.paint(p, 1000)
    p.end()
    box = next(h[0] for h in widget._hit if h[1] == "exon" and h[2][0]["group"] == "reference" and h[2][2] == (a, b))
    dup = next(h[0] for h in widget._hit if h[1] == "note" and "in the inserted bases" in str(h[2]))
    check(f"11. {strand} strand: an exon and its 100-base copy are drawn the same width, side by side",
          abs(box.width() - dup.width()) <= 1 and not box.intersects(dup),
          f"exon {box.x():.1f}+{box.width():.1f}, copy {dup.x():.1f}+{dup.width():.1f}")

helpers.finish()
