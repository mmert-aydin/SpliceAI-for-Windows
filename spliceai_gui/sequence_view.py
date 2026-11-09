"""The viewer's sequence comparison (spec 8.2 C): the cDNAs of the
reference, the wild type's prediction and the chosen patient candidate
aligned (one alignment, built from the two pairwise ones on the reference),
each base on its exon's colour, the first base of every exon underlined (the
junctions), gaps in red; and the proteins aligned, with the first changed
residue, the frameshifted stretch and the new stop marked. Copy buttons put
the sequences on the clipboard as FASTA.
"""
import html

from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QPushButton, QTextBrowser, QVBoxLayout, QWidget

from spliceai_pipeline import cdna, cdna_analysis

from . import theme
from .stacked_view import EXON_COLORS, NOVEL
from .i18n import t

LINE = 90


def merge(ref_a, other_a, ref_b, other_b):
    """Three aligned rows (reference, a, b) from the pairwise alignments
    (ref_a, other_a) and (ref_b, other_b) of the same reference: a gap the
    other alignment lacks becomes a gap in its rows too."""
    out = ["", "", ""]
    i = j = 0
    while i < len(ref_a) or j < len(ref_b):
        ga = i < len(ref_a) and ref_a[i] == "-"
        gb = j < len(ref_b) and ref_b[j] == "-"
        if ga and not gb:
            out[0] += "-"; out[1] += other_a[i]; out[2] += "-"; i += 1          # noqa: E702
        elif gb and not ga:
            out[0] += "-"; out[1] += "-"; out[2] += other_b[j]; j += 1          # noqa: E702
        elif i < len(ref_a) and j < len(ref_b):
            out[0] += ref_a[i] if not ga else "-"
            out[1] += other_a[i]
            out[2] += other_b[j]
            i, j = i + 1, j + 1
        elif i < len(ref_a):
            out[0] += ref_a[i]; out[1] += other_a[i]; out[2] += "-"; i += 1     # noqa: E702
        else:
            out[0] += ref_b[j]; out[1] += "-"; out[2] += other_b[j]; j += 1     # noqa: E702
    return out


def exon_identities(t, reference_exons):
    """For each base of a transcript's cDNA: (identity colour index, first base of an exon?)."""
    ident = []
    for k, e in enumerate(t.ref_exons):
        ident.append(next((n for n, c in enumerate(reference_exons, 1) if e and e[0] <= c[1] and c[0] <= e[1]), 0))
    out = []
    for k, (a, b) in enumerate(t.exons):
        out += [(ident[k], m == 0 and k > 0) for m in range(b - a + 1)]
    return out


class SequenceView(QWidget):
    def __init__(self, result, parent=None):
        super().__init__(parent)
        self.result = result
        layout = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel(t("Sample candidate:")))
        self.candidate = QComboBox()
        patients = [p for p in result.transcripts if p.kind == "Patient"]
        self.candidate.addItems([t("Sample {rank} (relative score {score:.3f}; {evidence})", rank=p.rank,
                                    score=p.relative, evidence=cdna_analysis.evidence(p, result.settings))
                                 for p in patients])
        self.candidate.currentIndexChanged.connect(self.refresh)
        row.addWidget(self.candidate)
        row.addWidget(QLabel(t(cdna_analysis.lead_text(result))))
        row.addStretch(1)
        copy_cdna = QPushButton(t("Copy cDNAs"))
        copy_cdna.clicked.connect(lambda: QGuiApplication.clipboard().setText(self.fasta(protein=False)))
        row.addWidget(copy_cdna)
        copy_protein = QPushButton(t("Copy proteins"))
        copy_protein.clicked.connect(lambda: QGuiApplication.clipboard().setText(self.fasta(protein=True)))
        row.addWidget(copy_protein)
        layout.addLayout(row)
        self.text = QTextBrowser()
        self.text.setFont(theme.mono_font(9))
        layout.addWidget(self.text)
        self.refresh()

    def chosen(self):
        by = {(t.kind, t.rank): t for t in self.result.transcripts}
        return by[("Reference", 0)], by.get(("WT predicted", 1)), by.get(("Patient", self.candidate.currentIndex() + 1))

    def fasta(self, protein):
        ref, wt, pat = self.chosen()
        out = []
        for t, name in ((ref, "Reference"), (wt, "WT_predicted_1"), (pat, f"Sample_{pat.rank}" if pat else "")):
            if t is None:
                continue
            seq = t.cdna.protein if protein else t.cdna.sequence
            out.append(f">{name} {self.result.gene.name}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)))
        return "\n".join(out) + "\n"

    def set_candidate(self, rank):
        self.candidate.setCurrentIndex(max(0, rank - 1))

    def refresh(self, *_):
        ref, wt, pat = self.chosen()
        if wt is None or pat is None:
            self.text.setPlainText(t("No predicted transcripts to compare."))
            return
        r1, w1, _ = cdna.pairwise(ref.cdna.sequence, wt.cdna.sequence)
        r2, p2, _ = cdna.pairwise(ref.cdna.sequence, pat.cdna.sequence)
        rows = merge(r1, w1, r2, p2)
        canonical = ref.exons
        styles = [exon_identities(t, canonical) for t in (ref, wt, pat)]
        names = ["Reference", "WT predicted 1", f"Sample {pat.rank}"]
        parts = [f"<p><b>cDNA</b> -- {html.escape(self.result.gene.name)}; "
                 + t("colours: the reference exon each base's exon matches (grey: novel), underlined: the first base "
                     "of an exon, red: a gap.") + "</p><pre>"]
        counters = [0, 0, 0]
        for start in range(0, len(rows[0]), LINE):
            for k in range(3):
                segment = rows[k][start:start + LINE]
                spans, current, buf = [], None, ""
                for ch in segment:
                    if ch == "-":
                        style = "background:#f4b6b6"
                    else:
                        ident, first = styles[k][counters[k]] if counters[k] < len(styles[k]) else (0, False)
                        color = EXON_COLORS[(ident - 1) % (len(EXON_COLORS) - 1)] if ident else NOVEL
                        style = f"background:{color}33" + (";text-decoration:underline;font-weight:bold" if first else "")
                        counters[k] += 1
                    if style != current and buf:
                        spans.append(f'<span style="{current}">{html.escape(buf)}</span>')
                        buf = ""
                    current = style
                    buf += ch
                if buf:
                    spans.append(f'<span style="{current}">{html.escape(buf)}</span>')
                parts.append(f"{names[k]:<15} {counters[k]:>6}  " + "".join(spans) + "\n")
            parts.append("\n")
        parts.append("</pre>")
        parts.append(self._protein_html(ref, pat))
        self.text.setHtml("".join(parts))

    def _protein_html(self, ref, pat):
        a = ref.cdna.protein + ("*" if ref.cdna.cds_end is not None else "")
        b = pat.cdna.protein + ("*" if pat.cdna.cds_end is not None else "")
        if not a or not b:
            return ("<p><b>" + t("Protein") + "</b>: " + (t("the ATG is lost in this candidate.") if not b else
                                                           t("no reference protein.")) + "</p>")
        top, mid, bottom = cdna.pairwise(a, b, protein=True)
        first = next((i for i, (x, y) in enumerate(zip(top, bottom)) if x != y), None)
        frameshift = cdna.summarize(pat.cdna, ref.cdna).frame == "frameshift"
        out = ["<p><b>" + t("Protein") + "</b> -- "
               + t("reference {ref} aa, sample {rank}: {sample} aa ({hgvs}); bold red: the first changed residue, "
                   "orange: frameshifted, black: the new stop.", ref=len(ref.cdna.protein), rank=pat.rank,
                   sample=len(pat.cdna.protein), hgvs=html.escape(cdna.hgvs_p(ref.cdna, pat.cdna))) + "</p><pre>"]
        for start in range(0, len(top), LINE):
            for name, seq in (("Reference", top), (f"Sample {pat.rank}", bottom)):
                cells = []
                for i, ch in enumerate(seq[start:start + LINE], start):
                    style = ""
                    if name != "Reference" and first is not None:
                        if i == first:
                            style = "color:#a3181f;font-weight:bold"
                        elif frameshift and i > first and ch != "*":
                            style = "background:#ffd9a8"
                        if ch == "*" and i == len(seq.rstrip("-")) - 1:
                            style = "background:#000;color:#fff"
                    cells.append(f'<span style="{style}">{html.escape(ch)}</span>' if style else html.escape(ch))
                out.append(f"{name:<15} " + "".join(cells) + "\n")
            out.append("\n")
        out.append("</pre>")
        return "".join(out)
