"""Spec 16 item 4: "Analyse transcript effect", the confirmation before
variants scored in the VCF tab go to the Sequence-to-cDNA tab.

One variant: the build, the gene (chosen here when the variant lies in
several), its canonical transcript, the variant, its genotype from the VCF,
the estimated time, and that only this variant is applied -- unless "Also
apply other variants in this gene" is ticked and some of them chosen (their
phase is unknown; they go on one haplotype together). Two to ten variants:
a queue, each run on its own and saved as its own results file. The whole
VCF is never processed.
"""
import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QPushButton, QTableWidget, QVBoxLayout,
)

from spliceai_pipeline import vcf_bridge
from spliceai_pipeline.cli import format_elapsed

from . import theme
from .i18n import N_, t

ONLY_THIS = t("Only this variant is applied: the sample is the reference with this one change. The VCF's other "
             "variants are not applied.")
PHASE_NOTE = (t("Their phase is unknown: the ones ticked are applied together with this variant on one haplotype, as "
              "if all were on the same allele."))
NONE_OTHER = N_("The VCF's results have no other scored variant in {gene}.")


def gene_text(g):
    strand = "+ strand" if g.strand == "+" else "− strand"
    c = g.canonical
    return f"{g.name} (chr{g.chrom}, {strand})" + (f", {c.accession}" if c else ", no NM_ transcript")


class TickList(QListWidget):
    """A list of checkable rows: a click anywhere on a row ticks it, not only on its small box."""

    def mouseReleaseEvent(self, event):
        item = self.itemAt(event.position().toPoint())
        before = item.checkState() if item else None
        super().mouseReleaseEvent(event)
        if item is not None and item.checkState() == before:
            item.setCheckState(Qt.Unchecked if before == Qt.Checked else Qt.Checked)


class TranscriptEffectDialog(QDialog):
    def __init__(self, parent, build, rows, genes_of, others_of, genotype_of, estimate_of, mode):
        """rows: the selected variants (merge.ScoreRow, one per variant, at most vcf_bridge.MAX_QUEUE).
        genes_of(row) -> ([transcripts.Gene], [gene names not in the table]); others_of(row, gene) -> the VCF's
        other variants in that gene; genotype_of(row) -> vcf_bridge.Genotype; estimate_of(gene, edits) ->
        (seconds, what it covers); mode: the VCF tab's raw / masked."""
        super().__init__(parent)
        self.build, self.rows, self.mode = build, rows, mode
        self._genes_of, self._others_of, self._estimate_of = genes_of, others_of, estimate_of
        self.genotypes = [genotype_of(r) for r in rows]
        self.found = [genes_of(r) for r in rows]
        self.setWindowTitle(t("Analyse transcript effect"))
        layout = QVBoxLayout(self)
        intro = QLabel(t("The Sequence-to-cDNA tab will build the sample from the reference ({build}) with the "
                         "variant given as coordinates, predict its transcripts and compare them with the reference "
                         "and the WT prediction.", build=build) if len(rows) == 1 else
                       t("The Sequence-to-cDNA tab will build the sample from the reference ({build}) with the {n} "
                         "variants, one at a time, given as coordinates, predict its transcripts and compare them "
                         "with the reference and the WT prediction.", build=build, n=len(rows)))
        intro.setWordWrap(True)
        layout.addWidget(intro)
        self.problem = QLabel("")
        self.problem.setWordWrap(True)
        self.problem.setStyleSheet(theme.STATUS_CRIT)
        buttons = QDialogButtonBox()
        self.compute_button = buttons.addButton(t("Compute"), QDialogButtonBox.AcceptRole)
        buttons.addButton(t("Cancel"), QDialogButtonBox.RejectRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        if len(rows) == 1:
            self._single(layout)
        else:
            self._queue(layout)
        layout.addWidget(self.problem)
        layout.addWidget(buttons)
        self._update()
        self.resize(720, self.sizeHint().height())

    # --- one variant ------------------------------------------------------------------------------
    def _single(self, layout):
        row, (genes, missing) = self.rows[0], self.found[0]
        form = QFormLayout()
        form.addRow(t("Build:"), QLabel(self.build))
        try:
            line = vcf_bridge.event_line(self.build, row.chrom, row.pos, row.ref, row.alt)
        except (ValueError, StopIteration) as exc:
            line = t("(can't be written as an event: {error})", error=exc)
        form.addRow(t("Variant:"), QLabel(f"{row.chrom}:{row.pos:,} {row.ref}>{row.alt}  —  "
                                         + t("entered as {line}", line=line)))
        self.gene_combo = QComboBox()
        self.gene_combo.addItems([gene_text(g) for g in genes])
        self.gene_combo.setVisible(len(genes) > 1)
        self.gene_combo.currentIndexChanged.connect(self._on_gene)
        self.gene_label = QLabel(gene_text(genes[0]) if len(genes) == 1 else "")
        self.gene_label.setVisible(len(genes) <= 1)
        box = QHBoxLayout()
        box.addWidget(self.gene_combo)
        box.addWidget(self.gene_label, stretch=1)
        form.addRow(t("Gene:"), box)
        if len(genes) > 1:
            form.addRow("", self._muted(t("The variant lies in {n} genes: choose the one to analyse.", n=len(genes))))
        if missing:
            form.addRow("", self._muted(t("Also scored in {genes}, which the Sequence-to-cDNA gene table doesn't have.",
                                          genes=", ".join(missing))))
        self.transcript_label = QLabel("")
        form.addRow(t("Canonical transcript:"), self.transcript_label)
        g = self.genotypes[0]
        form.addRow(t("Genotype:"), self._wrapped(t(g.text())))
        self.estimate_label = QLabel("")
        form.addRow(t("Estimated time:"), self.estimate_label)
        layout.addLayout(form)
        self.only_label = self._wrapped(ONLY_THIS)
        layout.addWidget(self.only_label)
        self.also = QCheckBox(t("Also apply other variants in this gene"))
        self.also.setChecked(False)
        self.also.toggled.connect(self._update)
        layout.addWidget(self.also)
        self._others = []
        self.others_list = TickList()
        self.others_list.setMaximumHeight(140)
        self.others_list.itemChanged.connect(self._update)
        layout.addWidget(self.others_list)
        self.phase_label = self._muted(PHASE_NOTE)
        layout.addWidget(self.phase_label)
        self._on_gene(0)

    def _on_gene(self, _index):
        gene = self.gene()
        if gene is None:
            return
        c = gene.canonical
        self.transcript_label.setText(f"{c.accession} ({gene.rule})" if c else t("none: nothing can be predicted"))
        self.others_list.blockSignals(True)
        self.others_list.clear()
        self._others = self._others_of(self.rows[0], gene)
        for r in self._others:
            score = "" if r.max_score is None else f", max score {r.max_score:.2f}"
            item = QListWidgetItem(f"{r.chrom}:{r.pos:,} {r.ref}>{r.alt}{score}")
            item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            item.setCheckState(Qt.Unchecked)
            self.others_list.addItem(item)
        self.others_list.blockSignals(False)
        self._update()

    def gene(self, i=0):
        genes = self.found[i][0]
        if not genes:
            return None
        if len(self.rows) == 1:
            return genes[max(self.gene_combo.currentIndex(), 0)]
        return genes[max(self.queue_combos[i].currentIndex(), 0)]

    def others(self):
        if len(self.rows) > 1 or not self.also.isChecked():
            return []
        return [r for k, r in enumerate(self._others) if self.others_list.item(k).checkState() == Qt.Checked]

    # --- a queue ------------------------------------------------------------------------------------
    def _queue(self, layout):
        layout.addWidget(self._wrapped(t("Each variant is run on its own: only that variant is applied, never the "
                                         "VCF's others. Each result is saved as its own results file in the folder "
                                         "below. At most {n} variants.", n=vcf_bridge.MAX_QUEUE)))
        self.table = QTableWidget(len(self.rows), 3)
        self.table.setHorizontalHeaderLabels([t("Variant"), t("Gene"), t("Genotype")])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.queue_combos = []
        for i, (row, (genes, _missing), g) in enumerate(zip(self.rows, self.found, self.genotypes)):
            self.table.setCellWidget(i, 0, QLabel(f" {row.chrom}:{row.pos:,} {row.ref}>{row.alt} "))
            combo = QComboBox()
            combo.addItems([gene_text(x) for x in genes] or [t("not in the Sequence-to-cDNA gene table")])
            combo.setEnabled(len(genes) > 1)
            combo.currentIndexChanged.connect(self._update)
            self.queue_combos.append(combo)
            self.table.setCellWidget(i, 1, combo)
            self.table.setCellWidget(i, 2, QLabel(f" {g.gt or t(g.meaning)} "))
        self.table.resizeColumnsToContents()
        layout.addWidget(self.table)
        row = QHBoxLayout()
        row.addWidget(QLabel(t("Save the results in:")))
        self.folder_edit = QLineEdit()
        self.folder_edit.textChanged.connect(self._update)
        row.addWidget(self.folder_edit, stretch=1)
        browse = QPushButton(t("Browse..."))
        browse.clicked.connect(self._on_browse)
        row.addWidget(browse)
        layout.addLayout(row)
        self.estimate_label = QLabel("")
        layout.addWidget(self.estimate_label)

    def _on_browse(self):
        folder = QFileDialog.getExistingDirectory(self, t("Save the results in"))
        if folder:
            self.folder_edit.setText(folder)

    def folder(self):
        return self.folder_edit.text().strip() if len(self.rows) > 1 else None

    # --- both ------------------------------------------------------------------------------------------
    def _update(self, *_):
        problems, seconds, what = [], 0.0, set()
        for i, row in enumerate(self.rows):
            gene = self.gene(i)
            if gene is None:
                problems.append(f"{row.chrom}:{row.pos:,}: " + t("none of its genes ({genes}) is in the "
                                                                   "Sequence-to-cDNA gene table for {build}.",
                                                                   genes=row.gene, build=self.build))
                continue
            if gene.canonical is None:
                problems.append(t("{gene} has no NM_ transcript to predict against.", gene=gene.name))
            s, w = self._estimate_of(gene, 1 + len(self.others()))
            seconds, what = seconds + s, what | {w}
        if len(self.rows) == 1:
            on = self.also.isChecked()
            self.others_list.setVisible(on and bool(self._others))
            self.phase_label.setVisible(on)
            gene = self.gene()
            self.phase_label.setText(PHASE_NOTE if self._others else t(NONE_OTHER, gene=gene.name if gene else t("it")))
            n = len(self.others())
            self.only_label.setText(ONLY_THIS if not n else t("This variant and the {n} ticked below are "
                                    "applied together; the VCF's other variants are not.", n=n))
        elif not self.folder():
            problems.append(t("Choose the folder the results files go in."))
        elif not os.path.isdir(self.folder()):
            problems.append(t("{folder} isn't a folder.", folder=self.folder()))
        self.estimate_label.setText((t("about {time} in all", time=format_elapsed(seconds)) if len(self.rows) > 1
                                     else t("about {time}", time=format_elapsed(seconds)))
                                    + (f" ({'; '.join(sorted(what))})" if what else ""))
        self.problem.setText("\n".join(problems))
        self.compute_button.setEnabled(not problems)

    def requests(self):
        """The vcf_bridge.Request of each variant, in order."""
        return [vcf_bridge.Request(self.build, self.gene(i), row, g, self.mode, self.others())
                for i, (row, g) in enumerate(zip(self.rows, self.genotypes))]

    @staticmethod
    def _wrapped(text):
        label = QLabel(text)
        label.setWordWrap(True)
        return label

    @classmethod
    def _muted(cls, text):
        label = cls._wrapped(text)
        label.setStyleSheet(theme.STATUS_MUTED)
        return label
