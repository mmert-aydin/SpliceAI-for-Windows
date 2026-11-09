"""The Sequence-to-cDNA tab (spec 6, 8.1, 8.3, 8.4, 9): a sample's sequence
for one gene, compared with the reference, and the transcripts SpliceAI's
predictions imply for both.

Inputs (build, gene, optional window, the sequence), "Check compatibility"
with its report and the differences table (which events build the sample's
haplotype), the parameters, Run / Cancel with the status area, the results
table, and its exports: CSV, the results bundle (also opened again, or
dropped onto the tab) and bedGraph tracks. All the work happens in
spliceai_pipeline.cdna_analysis, run by cdna_worker in the background.
"""
import html
import os
import time
from types import SimpleNamespace

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFrame, QGroupBox, QHBoxLayout,
    QHeaderView, QLabel, QLineEdit, QMenu, QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QSpinBox,
    QTableWidget, QTableWidgetItem, QTabWidget, QTextBrowser, QToolButton, QVBoxLayout, QWidget,
)

from spliceai_pipeline import (align, bedgraph, bundle, cdna_analysis as ca, haplotype, refdb, specificity, vcf_bridge,
                               transcripts, wt_concordance)
from spliceai_pipeline.cli import format_elapsed
from spliceai_pipeline.mane import find_cached_summary
from spliceai_pipeline.writer import delimiter_for_path, write_rows_ordered

from . import config, theme
from .cdna_help import (CHECK_HELP_TEXT, COORDINATES_HELP_TEXT, COORDS_HELP_TEXT, PARAMETERS_HELP_TEXT,
                        SCORES_HELP_TEXT, CHECKS_HELP_TEXT, HOW_TO_READ_HELP_TEXT, SEQUENCE_NOTE, intro_html)
from .cdna_worker import RunWorker, TaskWorker
from .cpu_meter import ProcessMeter
from .help_button import make_help_button, make_notice_button, notice_icon
from .refdb_wizard import DATABASE_HELP_TEXT, MissingDatabasePanel, live_seconds, mane_dir, refseq_dir
from .stacked_view import EVENT_COLORS
from .i18n import N_, t

SEQUENCE_SUFFIXES = (".fa", ".fasta", ".fna", ".txt")
EVENT_COLUMNS = ["Apply", "Type", "Location", "Reference position", "Size", "Inserted", "Origin", "Note"]
INSERTED_COLUMN = EVENT_COLUMNS.index("Inserted")
TOOLTIP_BASES = 3000    # more would not fit on a screen; the copy action gives them all


def fasta_text(name, sequence, width=60):
    """One FASTA record, `width` bases per line."""
    return f">{name}\n" + "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def sequence_tooltip(sequence, width=60):
    """The whole sequence, in lines; a very long one cut, saying so."""
    shown = sequence[:TOOLTIP_BASES]
    return "\n".join(shown[i:i + width] for i in range(0, len(shown), width)) + (
        "\n... " + t("{n:,} bp in all: right-click to copy them.", n=len(sequence)) if len(sequence) > len(shown)
        else "\n" + t("Right-click to copy."))


class SequenceEdit(QPlainTextEdit):
    """The paste box; a sequence file dropped onto it is loaded."""

    fileDropped = Signal(str)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            super().dragEnterEvent(event)

    dragMoveEvent = dragEnterEvent

    def dropEvent(self, event):
        files = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()]
        if files:
            event.acceptProposedAction()
            self.fileDropped.emit(os.path.normpath(files[0]))
        else:
            super().dropEvent(event)


class IntroDialog(QDialog):
    """Spec 15 item 6: the tab's introduction (cdna_help.intro_html), shown the first time the tab is opened, and
    from the "?" in its header; "Don't show again" is kept in the settings (cdna_show_intro)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t("The Sequence-to-cDNA Tab"))
        self.resize(780, 680)
        v = QVBoxLayout(self)
        self.text = QTextBrowser()
        # No database folder at all: an install without the USB drive -- say first how to get one.
        none = not any(os.path.isfile(os.path.join(refdb.default_root(), b, "manifest.json")) for b in ("hg19", "hg38"))
        self.text.setHtml(intro_html(first_steps=none))
        v.addWidget(self.text)
        row = QHBoxLayout()
        self.dont_show = QCheckBox(t("Don't show again"))
        self.dont_show.setChecked(not config.load()["cdna_show_intro"])
        self.dont_show.toggled.connect(self._save)
        row.addWidget(self.dont_show)
        row.addStretch(1)
        close = QPushButton(t("Close"))
        close.clicked.connect(self.close)
        row.addWidget(close)
        v.addLayout(row)

    def _save(self, on):
        settings = config.load()
        settings["cdna_show_intro"] = not on
        config.save(settings)


class CdnaTab(QWidget):
    def __init__(self, window):
        super().__init__(window)
        self.main = window              # MainWindow: reference FASTAs, MANE folder, the database manager
        self.tables, self.fastas, self.scores = {}, {}, {}
        self.genes, self.gene, self.compat, self.result = [], None, None, None
        self._workers = set()
        self._run_worker = None
        self._run_started = None
        self._skipped_panel = False
        self._vcf_queue, self._pending_vcf, self._run_request, self._auto_run = [], None, None, False
        self._vcf_folder, self._vcf_report = None, []
        self._meter = ProcessMeter(os.getpid())
        self.setAcceptDrops(True)
        settings = config.load()

        layout = QVBoxLayout(self)
        top = QHBoxLayout()                 # no research-use banner: the startup pop-up says it (spec 14 C7)
        self.tab_help_button = make_help_button(self, t("The Sequence-to-cDNA Tab"), "")
        # The way into the tab's introduction: yellow, at the left, so it is seen first.
        self.tab_help_button.setStyleSheet(
            f"QToolButton {{ color: {theme.INK}; background: {theme.YELLOW}; border: 1px solid {theme.YELLOW_LINE}; "
            "border-radius: 9px; font-weight: bold; font-size: 10pt; padding: 1px 8px; }")
        self.tab_help_button.setToolTip(t("The introduction to this tab: what it does, every parameter and result, "
                                        "and how to read the result."))
        self.tab_help_button.clicked.disconnect()
        self.tab_help_button.clicked.connect(lambda: self.show_intro())
        top.addWidget(self.tab_help_button)
        top.addStretch(1)
        layout.addLayout(top)
        self.intro = None
        self._intro_offered = False
        # Spec 5.3: which NM_ stands in for a missing MANE Select transcript, or
        # why nothing can be predicted for the gene.
        self.gene_note = QLabel("")
        self.gene_note.setWordWrap(True)
        self.gene_note.setStyleSheet(theme.pill("warn").replace("QPushButton, QLabel", "QLabel"))
        self.gene_note.setVisible(False)
        layout.addWidget(self.gene_note)

        self.missing_panel = MissingDatabasePanel(window.refdb)
        skip = QPushButton(t("Continue without it (genes computed live, a minute or two each)"))
        skip.clicked.connect(self._on_continue_without_database)
        self.missing_panel.layout().addWidget(skip, 0, Qt.AlignLeft)
        # Top-aligned: the tab is as tall as the VCF tab, and the panel spread over that height pushed its
        # buttons off the bottom of the window.
        layout.addWidget(self.missing_panel, 0, Qt.AlignTop)

        self.body = QWidget()
        body = QVBoxLayout(self.body)
        body.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.body)
        body.addWidget(self._build_inputs(settings))
        body.addWidget(self._build_check())
        body.addWidget(self._build_parameters(settings))
        body.addLayout(self._build_run_row())
        body.addWidget(self._build_results())

        self._timer = QTimer(self)
        self._timer.setInterval(1000)
        self._timer.timeout.connect(self._update_status)
        window.refdb.changed.connect(self._update_database_view)
        self._update_database_view()
        self._update_enabled()

    # --- building the widgets ----------------------------------------------------------------
    def _build_inputs(self, settings):
        group = QGroupBox(t("Sample sequence"))
        g = QVBoxLayout(group)
        row = QHBoxLayout()
        row.addWidget(QLabel(t("Genome build:")))
        self.build_combo = QComboBox()
        self.build_combo.addItems(["hg19", "hg38"])
        self.build_combo.setCurrentIndex(-1)
        self.build_combo.setPlaceholderText(t("choose"))
        self.build_combo.currentTextChanged.connect(self._on_build_changed)
        row.addWidget(self.build_combo)
        row.addStretch(1)
        g.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel(t("Gene:")))
        self.gene_edit = QLineEdit()
        self.gene_edit.setPlaceholderText(t("Gene symbol, e.g. BRCA2"))
        self.gene_edit.setMinimumWidth(self.gene_edit.fontMetrics().horizontalAdvance("Gene symbol, e.g. BRCA2") + 28)
        self.gene_edit.setToolTip(t("SpliceAI's gene symbol or RefSeq's: both are accepted (SpliceAI's annotation uses "
                                  "some older symbols, e.g. FAM132A for C1QTNF12)."))
        self.gene_edit.editingFinished.connect(self._on_gene_entered)
        row.addWidget(self.gene_edit)
        self.gene_choice = QComboBox()
        self.gene_choice.setVisible(False)
        self.gene_choice.currentIndexChanged.connect(self._on_gene_chosen)
        row.addWidget(self.gene_choice)
        row.addSpacing(16)
        row.addWidget(QLabel(t("Coordinates (optional):")))
        self.window_edit = QLineEdit()
        self.window_edit.setPlaceholderText(t("chr7:117,480,025-117,668,665 or chr7:117,559,590"))
        self.window_edit.setToolTip(t("A range restricts the alignment to it; a single position is searched around. "
                                    "1-based, inclusive."))
        self.window_edit.textChanged.connect(self._invalidate_check)
        self.window_edit.editingFinished.connect(self._on_coordinates_entered)
        row.addWidget(self.window_edit, stretch=1)
        row.addWidget(make_help_button(self, t("Coordinates"), COORDINATES_HELP_TEXT))
        g.addLayout(row)
        info = QHBoxLayout()
        self.gene_label = QLabel("")
        self.gene_label.setWordWrap(True)
        info.addWidget(self.gene_label, stretch=1)
        self.build_gene_button = QPushButton(t("Build this gene now"))
        self.build_gene_button.setToolTip(t("Adds this gene to the reference database first (in the background); "
                                          "otherwise Run computes it live."))
        self.build_gene_button.clicked.connect(self._on_build_gene_now)
        self.build_gene_button.setVisible(False)
        info.addWidget(self.build_gene_button)
        info.addWidget(make_help_button(self, t("Reference score database"), DATABASE_HELP_TEXT))
        g.addLayout(info)
        # Spec 14 C6: the sequence pasted, or built from events given as coordinates.
        self.input_tabs = QTabWidget()
        paste = QWidget()
        p = QVBoxLayout(paste)
        self.sequence_edit = SequenceEdit()
        self.sequence_edit.setPlaceholderText(t("Paste the sample's sequence (raw or FASTA), use \"Load FASTA...\", or "
                                              "drag a .fa / .fasta / .txt file here."))
        self.sequence_edit.setMinimumHeight(90)
        self.sequence_edit.textChanged.connect(self._invalidate_check)
        self.sequence_edit.fileDropped.connect(self.load_sequence_file)
        p.addWidget(self.sequence_edit)
        row = QHBoxLayout()
        load = QPushButton(t("Load FASTA..."))
        load.clicked.connect(self._on_load_fasta)
        row.addWidget(load)
        clear = QPushButton(t("Clear"))
        clear.clicked.connect(lambda: self.sequence_edit.clear())
        row.addWidget(clear)
        self.sequence_label = QLabel("")
        row.addWidget(self.sequence_label, stretch=1)
        p.addLayout(row)
        note = QLabel(SEQUENCE_NOTE)
        note.setWordWrap(True)
        note.setStyleSheet(theme.STATUS_MUTED)
        p.addWidget(note)
        self.input_tabs.addTab(paste, t("Paste sequence"))
        self.input_tabs.addTab(self._build_from_coordinates(), t("From coordinates"))
        g.addWidget(self.input_tabs)
        # The next step once a sequence is in: right under it, in the accent colour.
        row = QHBoxLayout()
        self.check_button = QPushButton(t("Check compatibility"))
        self.check_button.setObjectName("runButton")
        self.check_button.clicked.connect(self._on_check)
        row.addWidget(self.check_button)
        row.addWidget(make_help_button(self, t("Check compatibility"), CHECK_HELP_TEXT))
        row.addStretch(1)
        g.addLayout(row)
        return group

    def _build_from_coordinates(self):
        box = QWidget()
        v = QVBoxLayout(box)
        hint = QLabel(t("One event per line. Coordinates are 1-based and inclusive, as in VCF files and genome browsers."))
        hint.setWordWrap(True)
        v.addWidget(hint)
        self.coords_edit = QPlainTextEdit()
        self.coords_edit.setPlaceholderText(t("deletion chr17:17228330-17228462\n"
                                            "tandem duplication chr17:17228330-17228462 (add \"inverted\" for an "
                                            "inverted one)\n"
                                            "inversion chr17:17228330-17228462\n"
                                            "insertion chr17:17228330 ACGTACGT (inserted after that position)\n"
                                            "NC_000017.11:g.17228330_17228462del"))
        self.coords_edit.setMinimumHeight(90)
        v.addWidget(self.coords_edit)
        row = QHBoxLayout()
        self.generate_button = QPushButton(t("Generate sequence"))
        self.generate_button.setToolTip(t("The sample sequence of the gene's region: the reference with these events. "
                                        "It goes into the paste box, and the events into the differences table."))
        self.generate_button.clicked.connect(self._on_generate)
        row.addWidget(self.generate_button)
        self.save_fasta_button = QPushButton(t("Save as FASTA..."))
        self.save_fasta_button.setEnabled(False)
        self.save_fasta_button.clicked.connect(self._on_save_fasta)
        row.addWidget(self.save_fasta_button)
        row.addWidget(make_help_button(self, t("From coordinates"), COORDS_HELP_TEXT))
        self.coords_status = QLabel("")
        self.coords_status.setWordWrap(True)
        row.addWidget(self.coords_status, stretch=1)
        v.addLayout(row)
        return box

    def _build_check(self):
        group = QGroupBox(t("Compatibility"))
        g = QVBoxLayout(group)
        row = QHBoxLayout()
        self.check_status = QLabel("")
        self.check_status.setWordWrap(True)
        row.addWidget(self.check_status, stretch=1)
        icons = QVBoxLayout()
        self.caution_button = make_notice_button(self, "caution", t("Cautions about this input"))
        self.critical_button = make_notice_button(self, "critical", t("Critical problems"))
        icons.addWidget(self.caution_button)
        icons.addWidget(self.critical_button)
        icons.addStretch(1)
        row.addLayout(icons)
        g.addLayout(row)
        # Spec 17 item 3: a compact summary card; the technical lines under a collapsible "Details".
        self.report_label = QLabel("")
        self.report_label.setWordWrap(True)
        self.report_label.setTextFormat(Qt.RichText)
        self.report_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.report_label.setStyleSheet(f"background: {theme.SURFACE}; border: 1px solid {theme.LINE}; padding: 4px;")
        self.report_label.setVisible(False)
        g.addWidget(self.report_label)
        self.details_button = QToolButton()
        self.details_button.setText(t("Details"))
        self.details_button.setCheckable(True)
        self.details_button.setArrowType(Qt.RightArrow)
        self.details_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.details_button.setVisible(False)
        self.details_button.toggled.connect(self._on_details_toggled)
        g.addWidget(self.details_button)
        self.details_label = QLabel("")
        self.details_label.setWordWrap(True)
        self.details_label.setStyleSheet(theme.STATUS_MUTED)
        self.details_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.details_label.setVisible(False)
        g.addWidget(self.details_label)
        self.events_table = QTableWidget(0, len(EVENT_COLUMNS))
        self.events_table.setHorizontalHeaderLabels(EVENT_COLUMNS)
        header = self.events_table.horizontalHeader()
        for col in range(len(EVENT_COLUMNS) - 1):
            header.setSectionResizeMode(col, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(len(EVENT_COLUMNS) - 1, QHeaderView.Stretch)        # Note: the rest of the width
        self.events_table.setWordWrap(True)
        self.events_table.verticalHeader().setVisible(False)
        self.events_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.events_table.setMaximumHeight(240)
        self._event_rows = []           # table row -> index in compat.alignment.events
        self.events_table.itemChanged.connect(self._on_event_toggled)
        self.events_table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.events_table.customContextMenuRequested.connect(self._on_events_menu)
        g.addWidget(self.events_table)
        self.run_anyway = QCheckBox(t("Run anyway (the change touches the ATG or the transcript's start or end; the "
                                    "predicted transcripts may lose that anchor)"))
        self.run_anyway.toggled.connect(self._update_enabled)
        self.run_anyway.setVisible(False)
        g.addWidget(self.run_anyway)
        return group

    def _build_parameters(self, settings):
        group = QGroupBox(t("Parameters"))
        rows = QVBoxLayout(group)
        row = QHBoxLayout()
        rows.addLayout(row)
        row.addWidget(QLabel(t("Site threshold:")))
        self.threshold = QDoubleSpinBox()
        self.threshold.setRange(0.01, 0.99)
        self.threshold.setSingleStep(0.01)
        self.threshold.setValue(settings["cdna_threshold"])
        row.addWidget(self.threshold)
        row.addWidget(make_help_button(self, t("Parameters"), PARAMETERS_HELP_TEXT))
        row.addSpacing(12)
        row.addWidget(QLabel(t("Minimum internal exon:")))
        self.min_exon = QSpinBox()
        self.min_exon.setRange(1, 1000)
        self.min_exon.setValue(settings["cdna_min_exon"])
        self.min_exon.setSuffix(t(" nt"))
        row.addWidget(self.min_exon)
        row.addWidget(QLabel(t("Minimum intron:")))
        self.min_intron = QSpinBox()
        self.min_intron.setRange(1, 10000)
        self.min_intron.setValue(settings["cdna_min_intron"])
        self.min_intron.setSuffix(t(" nt"))
        row.addWidget(self.min_intron)
        row.addStretch(1)
        row = QHBoxLayout()
        rows.addLayout(row)
        row.addWidget(QLabel(t("Apply by default from:")))
        self.apply_min_size = QSpinBox()
        self.apply_min_size.setRange(1, 100000)
        self.apply_min_size.setValue(settings["cdna_apply_min_size"])
        self.apply_min_size.setSuffix(t(" bp"))
        self.apply_min_size.setToolTip(t("Differences this long or longer (and every inversion and duplication) are "
                                       "ticked under Apply when the check finds them."))
        self.apply_min_size.valueChanged.connect(self._reapply_defaults)
        row.addWidget(self.apply_min_size)
        self.advanced = QPushButton(t("Advanced"))
        self.advanced.setCheckable(True)
        row.addWidget(self.advanced)
        row.addStretch(1)
        row = QHBoxLayout()
        rows.addLayout(row)
        self.top_label = QLabel(t("Transcripts per sequence:"))
        self.top = QSpinBox()
        self.top.setRange(1, 10)
        self.top.setValue(settings["cdna_top"])
        for w in (self.top_label, self.top):
            w.setVisible(False)
            row.addWidget(w)
        self.weak_label = QLabel(t("Weak evidence below delta:"))
        self.weak_delta = QDoubleSpinBox()
        self.weak_delta.setRange(0.01, 0.99)
        self.weak_delta.setSingleStep(0.05)
        self.weak_delta.setValue(settings["cdna_weak_delta"])
        self.weak_delta.setToolTip(t("A sample candidate whose differing splice sites all changed by less than this is "
                                   "marked \"weak evidence\" (see Reading the results). It doesn't change the ranking."))
        self.specificity_n_label = QLabel(t("Random stand-ins:"))
        self.specificity_n = QSpinBox()
        self.specificity_n.setRange(5, 200)
        self.specificity_n.setValue(settings["cdna_specificity_n"])
        self.specificity_n.setToolTip(t("How many random stand-ins per change the sequence-specificity control uses."))
        advanced = (self.top_label, self.top, self.weak_label, self.weak_delta, self.specificity_n_label,
                    self.specificity_n)
        for w in advanced[2:]:
            w.setVisible(False)
            row.addWidget(w)
        self.advanced.toggled.connect(lambda on: [w.setVisible(on) for w in advanced])
        row.addStretch(1)
        return group

    def _build_run_row(self):
        box = QVBoxLayout()
        row = QHBoxLayout()
        self.run_button = QPushButton(t("Run"))
        self.run_button.setObjectName("runButton")
        self.run_button.clicked.connect(self._on_run)
        row.addWidget(self.run_button)
        self.cancel_button = QPushButton(t("Cancel"))
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self._on_cancel)
        row.addWidget(self.cancel_button)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1000)
        row.addWidget(self.progress_bar, stretch=1)
        self.progress_label = QLabel("")
        row.addWidget(self.progress_label, stretch=1)
        box.addLayout(row)
        self.status_label = QLabel("")
        self.status_label.setStyleSheet(theme.STATUS_MUTED)
        self.status_label.setWordWrap(True)
        box.addWidget(self.status_label)
        self.estimate_label = QLabel("")
        self.estimate_label.setStyleSheet(theme.STATUS_MUTED)
        box.addWidget(self.estimate_label)
        return box

    def _build_results(self):
        group = QGroupBox(t("Results"))
        g = QVBoxLayout(group)
        row = QHBoxLayout()
        self.show_alternatives = QCheckBox(t("WT alternatives (2 to 5)"))
        self.show_alternatives.setToolTip(t("Also shows the wild type's alternative transcripts, ranks 2 to 5."))
        self.show_alternatives.toggled.connect(self._fill_table)
        row.addWidget(self.show_alternatives)
        row.addStretch(1)
        g.addLayout(row)
        # Reading the results, in three: a row of its own, so it never widens the tab.
        helps = QHBoxLayout()
        for title, text in ((N_("Scores and evidence"), SCORES_HELP_TEXT), (N_("Optional checks"), CHECKS_HELP_TEXT),
                            (N_("How to read the result"), HOW_TO_READ_HELP_TEXT)):
            helps.addWidget(make_help_button(self, t(title), text, label=t(title)))
        helps.addStretch(1)
        g.addLayout(helps)
        tools, row = row, QHBoxLayout()        # the viewer and the control on the first row, the files below
        self.export_csv_button = QPushButton(t("Export CSV..."))
        self.export_csv_button.setToolTip(t("The rows shown, as a CSV or tab-separated file."))
        self.export_csv_button.clicked.connect(self._on_export_csv)
        row.addWidget(self.export_csv_button)
        self.download_button = QPushButton(t("Download results..."))
        self.download_button.setToolTip(t("Everything, in one file this tab (and the viewer) can open again without "
                                        "SpliceAI, the FASTA or the database."))
        self.download_button.clicked.connect(self._on_download)
        row.addWidget(self.download_button)
        self.open_button = QPushButton(t("Open calculated results..."))
        self.open_button.setToolTip(t("Opens a results file saved with \"Download results\" -- or drag it onto this tab."))
        self.open_button.clicked.connect(self._on_open_results)
        row.addWidget(self.open_button)
        self.viewer_button = QPushButton(t("Open viewer"))
        self.viewer_button.setObjectName("viewerButton")
        self.viewer_button.setToolTip(t("The base-level view, the exon schematic and the sequence comparison -- or click "
                                      "a row."))
        self.viewer_button.clicked.connect(lambda: self.open_viewer())
        tools.addWidget(self.viewer_button)
        self.bedgraph_button = QPushButton(t("Export bedGraph..."))
        self.bedgraph_button.setToolTip(t("Reference, sample and delta tracks for IGV or the UCSC browser."))
        self.bedgraph_button.clicked.connect(self._on_export_bedgraph)
        row.addWidget(self.bedgraph_button)
        self.specificity_button = QPushButton(t("Sequence-specificity control"))
        self.specificity_button.setToolTip(t("Optional: repeats the prediction with random stand-ins for each change "
                                           "(see Reading the results). Runs only when clicked."))
        self.specificity_button.clicked.connect(self._on_specificity)
        tools.addWidget(self.specificity_button)
        self.breakpoints_button = QPushButton(t("Breakpoint sensitivity"))
        self.breakpoints_button.setToolTip(t("Optional: moves each breakpoint by 10 and 50 bases either way and says "
                                           "whether the top transcript changes. Runs only when clicked."))
        self.breakpoints_button.clicked.connect(self._on_breakpoints)
        tools.addWidget(self.breakpoints_button)
        row.addStretch(1)
        g.addLayout(row)
        self.lead_label = QLabel("")
        self.lead_label.setStyleSheet(theme.STATUS_MUTED)
        self.local_label = QLabel("")
        self.local_label.setWordWrap(True)
        self.local_label.setStyleSheet(theme.STATUS_MUTED)
        # Spec 16 item 4: a run sent from the VCF tab says which variant, its genotype, and the VCF tab's scores
        # next to the engine's; a "?" says why they differ where they do.
        self.vcf_header = QWidget()
        vh = QHBoxLayout(self.vcf_header)
        vh.setContentsMargins(0, 0, 0, 0)
        self.vcf_label = QLabel("")
        self.vcf_label.setWordWrap(True)
        self.vcf_label.setTextFormat(Qt.RichText)
        self.vcf_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        vh.addWidget(self.vcf_label, stretch=1)
        self.vcf_why = make_notice_button(self, "caution", t("Why the VCF tab's scores differ"))
        vh.addWidget(self.vcf_why, alignment=Qt.AlignTop)
        self.vcf_header.setVisible(False)
        g.addWidget(self.vcf_header)
        self.results_table = QTableWidget(0, len(ca.COLUMNS))
        self.results_table.setHorizontalHeaderLabels(ca.COLUMNS)
        self.results_table.horizontalHeader().setSectionResizeMode(len(ca.COLUMNS) - 1, QHeaderView.Stretch)
        self.results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.results_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.results_table.setMinimumHeight(260)
        self.results_table.cellClicked.connect(self._on_row_clicked)
        g.addWidget(self.results_table)
        # Spec 16 item 1: one line per applied event under the table; a click opens the viewer there.
        self.events_label = QLabel("")
        self.events_label.setWordWrap(True)
        self.events_label.setTextFormat(Qt.RichText)
        self.events_label.setTextInteractionFlags(Qt.TextBrowserInteraction)
        self.events_label.linkActivated.connect(self._on_event_link)
        g.addWidget(self.events_label)
        g.addWidget(self.lead_label)
        g.addWidget(self.local_label)
        self.specificity_label = QLabel("")
        self.specificity_label.setWordWrap(True)
        g.addWidget(self.specificity_label)
        self.breakpoints_label = QLabel("")
        self.breakpoints_label.setWordWrap(True)
        g.addWidget(self.breakpoints_label)
        self.viewer = None
        return group

    # --- database and genes ---------------------------------------------------------------------------
    def _update_database_view(self):
        m = self.main.refdb
        checking = any(s["state"] == "checking" for s in m.status.values())
        none = not checking and not any(m.usable(b) for b in m.status) and m.proc is None
        self.missing_panel.setVisible(none and not self._skipped_panel)
        self.body.setVisible(not none or self._skipped_panel)
        for build, scores in self.scores.items():
            scores.reader = m.reader(build)     # keeps its loaded models
        self._update_gene_label()
        self._update_status()

    def show_intro(self):
        if self.intro is None:
            self.intro = IntroDialog(self)
        self.intro.show()
        self.intro.raise_()
        return self.intro

    def showEvent(self, event):
        super().showEvent(event)
        if not self._intro_offered:
            self._intro_offered = True
            if config.load()["cdna_show_intro"]:
                QTimer.singleShot(0, self.show_intro)

    def _on_continue_without_database(self):
        self._skipped_panel = True
        self._update_database_view()

    def _fasta(self, build):
        path = self.main.fasta_edits[build].text().strip()
        if not path or not os.path.isfile(path):
            raise ValueError(f"No {build} reference FASTA is set (VCF tab > Settings > Advanced settings).")
        if self.fastas.get(build, (None,))[0] != path:
            from pyfaidx import Fasta

            self.fastas[build] = (path, Fasta(path))
        return self.fastas[build][1]

    def _scores(self, build):
        if build not in self.scores:
            m = self.main.refdb
            self.scores[build] = refdb.ReferenceScores(
                build, self.main.fasta_edits[build].text().strip(), reader=m.reader(build),
                cache_dir=os.path.join(m.cache_dir(), build), refseq_dir=refseq_dir(),
                mane_dir=mane_dir(self.main._mane_dir()))
        return self.scores[build]

    def _table(self, build):
        if build not in self.tables:
            mane = find_cached_summary(mane_dir(self.main._mane_dir()))
            if not mane:
                raise ValueError("The MANE Select summary file isn't there (VCF tab > Settings > Advanced settings).")
            self.tables[build] = ca.GeneTable(build, refseq_dir(), mane)
        return self.tables[build]

    def _start(self, fn, done, failed=None):
        worker = TaskWorker(fn, self)
        self._workers.add(worker)
        worker.finished_ok.connect(done)
        worker.failed.connect(failed or (lambda msg: QMessageBox.warning(self, t("Sequence to cDNA"), msg)))
        worker.finished.connect(lambda: self._workers.discard(worker))
        worker.start()
        return worker

    def _on_build_changed(self, _build):
        self.genes, self.gene = [], None
        self.gene_choice.blockSignals(True)
        self.gene_choice.clear()                 # the other build's genes
        self.gene_choice.blockSignals(False)
        self.gene_choice.setVisible(False)
        self._invalidate_check()
        if not self._refseq_missing(self.build_combo.currentText()) and self.gene_edit.text().strip():
            self._on_gene_entered()
        self._update_enabled()

    def _refseq_missing(self, build):
        """Says so in place of the gene when the build's RefSeq table is missing (no empty gene list)."""
        problem = build and transcripts.missing_table(build, refseq_dir())
        if problem:
            self._gene_error(problem)
        return bool(problem)

    def _gene_error(self, msg):
        self.gene_label.setText(t(str(msg)))
        self.gene_label.setStyleSheet(theme.STATUS_CRIT)

    def _on_generate(self):
        """Spec 14 C6: the sample sequence from events given as coordinates; the events become the differences."""
        build, gene, text = self.build_combo.currentText(), self.gene, self.coords_edit.toPlainText()
        if not build or gene is None:
            QMessageBox.warning(self, t("Generate sequence"), t("Choose the genome build and the gene first."))
            return
        try:
            record = self._fasta(build)["chr" + gene.chrom]
        except ValueError as exc:
            QMessageBox.warning(self, t("Generate sequence"), t(str(exc)))
            return
        fetch = lambda a, b: record[a - 1:b].seq
        self._invalidate_check()                    # the input changes: the last check no longer holds
        self.coords_status.setText(t("Building the sample sequence..."))
        self.coords_status.setStyleSheet(theme.STATUS_MUTED)
        self._start(lambda: ca.from_coordinates(text, gene, build, fetch, len(record)), self._generated,
                    self._generate_failed)

    def _generate_failed(self, message):
        self.coords_status.setText(t(message))
        self.coords_status.setStyleSheet(theme.STATUS_CRIT)
        if self._auto_run:
            self._vcf_stopped(message)

    def _generated(self, compat):
        g, (a, b) = self.gene, compat.region
        self.sequence_edit.blockSignals(True)       # the events are the differences already: no check undone
        self.sequence_edit.setPlainText(fasta_text(f"{g.name}_{self.build_combo.currentText()}_sample "
                                                   f"chr{g.chrom}:{a}-{b} from coordinates", compat.sequence))
        self.sequence_edit.blockSignals(False)
        n = len(compat.alignment.events)
        self.sequence_label.setText(t("Generated from {n} event(s) given as coordinates ({bases:,} bases)", n=n,
                                       bases=len(compat.sequence)))
        self.coords_status.setText("✓ " + t("{n} event(s) applied to {region}; the sequence is in the paste box.",
                                             n=n, region=f"chr{g.chrom}:{a:,}-{b:,}"))
        self.coords_status.setStyleSheet(theme.STATUS_OK)
        self.save_fasta_button.setEnabled(True)
        self._check_done(compat)
        if self._auto_run:
            QTimer.singleShot(0, self._auto_start)

    def _on_save_fasta(self):
        name = self._default_name().replace("_sequence-to-cdna", "") + "_sample.fa"
        path, _ = QFileDialog.getSaveFileName(self, t("Save as FASTA"), name, t("FASTA (*.fa *.fasta)"))
        if path:
            with open(path, "w", encoding="ascii", newline="") as fh:
                fh.write(self.sequence_edit.toPlainText().strip() + "\n")

    def _on_coordinates_entered(self):
        """Spec 14 C5: with no gene entered, the gene(s) at the coordinates are suggested."""
        build, text = self.build_combo.currentText(), self.window_edit.text().strip()
        if not build or not text or self.gene_edit.text().strip() or self._refseq_missing(build):
            return
        try:
            chrom, a, b, _single = align.parse_coordinates(text)
        except ValueError as exc:
            self._gene_error(str(exc))
            return
        self.gene_label.setText(t("Looking for genes there..."))
        self.gene_label.setStyleSheet("")
        self._start(lambda: self._table(build).at(chrom, a, b), lambda genes: self._suggest(build, chrom, a, b, genes),
                    self._gene_error)

    def _suggest(self, build, chrom, a, b, genes):
        where = f"chr{chrom}:{a:,}" + (f"-{b:,}" if b != a else "")
        if build != self.build_combo.currentText() or self.gene_edit.text().strip():
            return
        if not genes:
            self._gene_error(t("No gene of SpliceAI's {build} gene table at {where}.", build=build, where=where))
            return
        self.gene_edit.blockSignals(True)
        self.gene_edit.setText(genes[0].name)
        self.gene_edit.blockSignals(False)
        self._genes_found(build, genes[0].name, genes)
        names = ", ".join(g.name for g in genes)
        self.gene_label.setText((t("Genes at {where}: {names}; {gene} is chosen (change it in the list). ",
                                   where=where, names=names, gene=genes[0].name) if len(genes) > 1 else
                                 t("Gene at {where}: {names}; {gene} is chosen. ", where=where, names=names,
                                   gene=genes[0].name)) + self.gene_label.text())

    def _on_gene_entered(self):
        build, symbol = self.build_combo.currentText(), self.gene_edit.text().strip()
        if not build or not symbol or self._refseq_missing(build):
            return
        self.gene_label.setText(t("Looking the gene up..."))
        self.gene_label.setStyleSheet("")
        self._start(lambda: self._table(build).find(symbol), lambda genes: self._genes_found(build, symbol, genes),
                    self._gene_error)

    def _genes_found(self, build, symbol, genes):
        if build != self.build_combo.currentText() or symbol != self.gene_edit.text().strip():
            return
        self.genes = genes
        self.gene_choice.blockSignals(True)
        self.gene_choice.clear()
        for g in genes:
            self.gene_choice.addItem(f"{g.name} ({g.canonical.accession if g.canonical else t('no NM_')}, chr{g.chrom} "
                                     f"{g.start:,}-{g.end:,})")
        self.gene_choice.blockSignals(False)
        self.gene_choice.setVisible(len(genes) > 1)
        self.gene = genes[0] if genes else None
        self._invalidate_check()
        self._update_gene_label()

    def _on_gene_chosen(self, index):
        if 0 <= index < len(self.genes):
            self.gene = self.genes[index]
            self._invalidate_check()
            self._update_gene_label()

    def _note(self, text):
        self.gene_note.setText(t(text) if text else "")
        self.gene_note.setVisible(bool(text))

    def _update_gene_label(self):
        g, build = self.gene, self.build_combo.currentText()
        self.build_gene_button.setVisible(False)
        self._note(None)
        self._no_prediction = False
        if g is None:
            if self.genes == [] and self.gene_edit.text().strip() and build:
                self.gene_label.setText(t("No gene \"{name}\" in SpliceAI's {build} gene table.",
                                          name=self.gene_edit.text().strip(), build=build))
                self.gene_label.setStyleSheet(theme.STATUS_CRIT)
            return
        c = g.canonical
        text = f"{g.name}" + (f" ({c.symbol})" if c and c.symbol != g.name else "")
        strand = t("+ strand") if g.strand == "+" else t("− strand")     # spec 15 item 1
        text += f", chr{g.chrom}, {strand}: " + (f"{c.accession} ({g.rule})" if c else t("no NM_ transcript"))
        if g.other_mane:
            text += ". " + t("The locus is shared with {genes}.", genes=", ".join(g.other_mane))
        known = wt_concordance.for_gene(build, g.name) if c else None
        if known:
            # Spec 14 C1: the genome-wide run, with the default parameters.
            d = known["differing"]
            listed = ", ".join(map(str, d[:12])) + (t(" and {n} more", n=len(d) - 12) if len(d) > 12 else "")
            text += (". " + t("WT concordance for this gene: {same}/{exons} exons", same=known["identical"],
                              exons=known["exons"])
                     + (t(", the whole transcript exactly") if known["whole_exact"] else "")
                     + ((t("; predicted differently: exons {list}", list=listed) if len(d) > 1 else
                         t("; predicted differently: exon {list}", list=listed)) if d else "")
                     + t(" (default parameters)"))
        m = self.main.refdb
        reader = m.reader(build) if build else None
        in_db = reader is not None and reader.plan(g.name) is not None
        kind, note = self.tables[build].reference_note(g) if build in self.tables else (None, None)
        self._no_prediction = kind == "none"
        if in_db:
            text += ". " + t("Reference values: from the database.")
        elif not self._no_prediction:
            # Before Run: how long the live computation of this gene will take.
            bases = g.mask[1] - g.mask[0] + 1 + 2 * refdb.MARGIN
            text += ". " + t("Reference values: computed live (not in the database), about {time} for its "
                             "{bases:,} bases, once; then cached.", time=format_elapsed(live_seconds(bases)), bases=bases)
        self.gene_label.setText(text)
        self.gene_label.setStyleSheet("")
        self._note(note)
        self.build_gene_button.setVisible(not in_db and not self._no_prediction
                                          and m.status.get(build, {}).get("state") != "refused")
        self._update_enabled()

    def _on_build_gene_now(self):
        if self.gene is not None:
            self.main.refdb.request_gene(self.build_combo.currentText(), self.gene.name)
            self.build_gene_button.setVisible(False)
            self.gene_label.setText(self.gene_label.text() + " " + t("Building it now (Reference database at the top)."))

    # --- the sequence ---------------------------------------------------------------------------------
    def _on_load_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, t("Load the sample's sequence"), "",
                                              t("Sequence files (*.fa *.fasta *.fna *.txt);;All files (*)"))
        if path:
            self.load_sequence_file(path)

    def load_sequence_file(self, path):
        if path.lower().endswith(bundle.EXTENSION):
            self.open_results(path)
            return
        try:
            with open(path, encoding="utf-8-sig", errors="replace") as fh:
                self.sequence_edit.setPlainText(fh.read())
            self.sequence_label.setText(os.path.basename(path))
        except OSError as exc:
            QMessageBox.warning(self, t("Load FASTA"), t("{path} can't be read: {error}", path=path, error=exc))

    def handle_drop(self, path):
        """A file dropped anywhere on the window while this tab is shown."""
        if path.lower().endswith(bundle.EXTENSION):
            self.open_results(path)
        elif path.lower().endswith(SEQUENCE_SUFFIXES):
            self.load_sequence_file(path)
        else:
            QMessageBox.information(self, t("Can't use this file"), t("{name} can't be loaded here.\n\nDrag in a sequence (.fa, .fasta, .txt) or a results file "
                                      "(*{ext}).", name=os.path.basename(path), ext=bundle.EXTENSION))

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        files = [u.toLocalFile() for u in event.mimeData().urls() if u.isLocalFile()]
        if files:
            event.acceptProposedAction()
            self.handle_drop(os.path.normpath(files[0]))

    # --- check compatibility ----------------------------------------------------------------------------
    def _invalidate_check(self, *_):
        self.compat = None
        self.check_status.setText("")
        self.report_label.setText("")
        self.details_label.setText("")
        for widget in (self.report_label, self.details_button, self.details_label):
            widget.setVisible(False)
        self.events_table.setRowCount(0)
        for button in (self.caution_button, self.critical_button):
            button.setVisible(False)
        self.run_anyway.setVisible(False)
        self.run_anyway.setChecked(False)
        self.estimate_label.setText("")
        self._update_enabled()

    def _on_check(self):
        build, gene, text = self.build_combo.currentText(), self.gene, self.sequence_edit.toPlainText()
        window_text = self.window_edit.text().strip()
        try:
            window = None
            position = None
            if window_text:
                chrom, a, b, single = align.parse_coordinates(window_text)
                if chrom != gene.chrom:
                    raise ValueError(f"The coordinates are on chr{chrom}; {gene.name} is on chr{gene.chrom}.")
                if single:
                    position = a
                else:
                    window = (a, b)
            fasta = self._fasta(build)
        except ValueError as exc:
            QMessageBox.warning(self, t("Check compatibility"), t(str(exc)))
            return
        record = fasta["chr" + gene.chrom]
        fetch = lambda a, b: record[a - 1:b].seq
        other = self._other_build(build, gene)
        self.check_button.setEnabled(False)
        self.check_status.setText(t("Reading and aligning the sequence..."))
        self.check_status.setStyleSheet(theme.STATUS_MUTED)

        def work():
            compat = ca.check_compatibility(text, gene, fetch, len(record), window=window, position=position)
            if other is not None:
                compat.report = align.compatibility(compat.alignment, gene.mask[0], gene.mask[1],
                                                    spliceai=(gene.start, gene.end),
                                                    canonical=(gene.canonical.start, gene.canonical.end)
                                                    if gene.canonical else None, other=other(compat.sequence))
                if compat.search:
                    compat.report["lines"].insert(0, compat.search)
            return compat

        def failed(message):
            self.check_button.setEnabled(True)
            self.check_status.setText(message)
            self.check_status.setStyleSheet(theme.STATUS_CRIT)

        self._start(work, self._check_done, failed)

    def _other_build(self, build, gene):
        """A function aligning a sequence to the same gene on the other build, if its FASTA is set; else None."""
        other = "hg19" if build == "hg38" else "hg38"
        path = self.main.fasta_edits[other].text().strip()
        if not path or not os.path.isfile(path) or not os.path.exists(path + ".fai"):
            return None

        def align_other(sequence):
            try:
                table = self._table(other)
                hits = table.find(gene.canonical.symbol if gene.canonical else gene.name) or table.find(gene.name)
                if not hits:
                    return None
                g = hits[0]
                record = self._fasta(other)["chr" + g.chrom]
                start, end = align.analysed_region(g.mask[0], g.mask[1], len(record))
                a = align.align(sequence, record[start - 1:end].seq.upper(), start,
                                lambda c, x, y: record[x - 1:y].seq.upper(), "chr" + g.chrom)
                return other, a
            except (ValueError, KeyError, OSError):
                return None

        return align_other

    def _check_done(self, compat):
        self.check_button.setEnabled(True)
        self.run_anyway.setChecked(False)          # asked again for every new check
        if compat.found:
            # A single position was given: the range found goes into the field (without undoing the check).
            self.window_edit.blockSignals(True)
            self.window_edit.setText(f"chr{self.gene.chrom}:{compat.found[0]:,}-{compat.found[1]:,}")
            self.window_edit.blockSignals(False)
        self.compat = compat
        a, report = compat.alignment, compat.report
        card = ca.compat_card(compat, self.gene)                   # notices: the "?" and "!" icons
        self.report_label.setText("<table cellspacing='0' cellpadding='1'>" + "".join(
            f"<tr><td nowrap style='color: {theme.INK_3}; padding-right: 8px;'>{html.escape(t(k))}</td>"
            f"<td>{html.escape(t(v))}</td></tr>" for k, v in card) + "</table>")
        self.report_label.setVisible(bool(card))
        self.details_label.setText("\n".join(t(line) for line in report["lines"]))
        self.details_button.setVisible(bool(report["lines"]))
        self.details_label.setVisible(self.details_button.isChecked() and bool(report["lines"]))
        self.check_status.setText(t("✗ The sequence can't be used.") if report["problems"] else t("✓ Checked."))
        self.check_status.setStyleSheet(theme.STATUS_CRIT if report["problems"] else theme.STATUS_OK)
        self.events_table.blockSignals(True)
        self.events_table.setRowCount(len(a.events))
        places = ca.difference_places(compat, self.gene)
        # Most important first: the applied differences, then the larger ones, then in order along the gene.
        self._event_rows = sorted(range(len(a.events)), key=lambda k: (not a.events[k].apply, -a.events[k].size,
                                                                        a.events[k].start))
        right = Qt.AlignRight | Qt.AlignVCenter
        for i, k in enumerate(self._event_rows):
            e, (location, origin) = a.events[k], places[k]
            apply_item = QTableWidgetItem()
            apply_item.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            apply_item.setCheckState(Qt.Checked if e.apply else Qt.Unchecked)
            self.events_table.setItem(i, 0, apply_item)
            where = f"{e.start:,}" if e.end <= e.start else f"{e.start:,}..{e.end:,}"
            if e.end < e.start:
                where = f"before {e.start:,}"
            values = [e.label or e.kind, location, where, f"{e.size:,}", ca.short_sequence(e.inserted), origin, e.note]
            for col, value in enumerate(values, 1):
                item = QTableWidgetItem(value)
                item.setToolTip(value)
                if EVENT_COLUMNS[col] in ("Reference position", "Size"):
                    item.setTextAlignment(right)
                self.events_table.setItem(i, col, item)
            if e.inserted:
                self.events_table.item(i, INSERTED_COLUMN).setToolTip(sequence_tooltip(e.inserted))
        self.events_table.blockSignals(False)
        self.events_table.resizeRowsToContents()
        self._update_chosen()

    def _on_details_toggled(self, on):
        self.details_button.setArrowType(Qt.DownArrow if on else Qt.RightArrow)
        self.details_label.setVisible(on and bool(self.details_label.text()))

    def event_of_row(self, row):
        """The difference (align.Event) on the differences table's row (the table is sorted by importance)."""
        return self.compat.alignment.events[self._event_rows[row]]

    def copy_inserted(self, row):
        """Puts the inserted bases of the differences table's row on the clipboard."""
        QApplication.clipboard().setText(self.event_of_row(row).inserted)

    def _on_events_menu(self, pos):
        row = self.events_table.rowAt(pos.y())
        if self.compat is None or row < 0 or not self.event_of_row(row).inserted:
            return
        menu = QMenu(self)
        n = len(self.event_of_row(row).inserted)
        copy_action = menu.addAction(t("Copy inserted sequence ({n:,} bp)", n=n))
        if menu.exec(self.events_table.viewport().mapToGlobal(pos)) == copy_action:
            self.copy_inserted(row)

    def _on_event_toggled(self, item):
        if self.compat is None or item.column() != 0:
            return
        self.event_of_row(item.row()).apply = item.checkState() == Qt.Checked
        self._update_chosen()

    def _reapply_defaults(self, size):
        if self.compat is None:
            return
        for e in self.compat.alignment.events:
            e.apply = align.default_apply(e, size)
        self._check_done(self.compat)

    def _update_chosen(self):
        """The chosen events' check, its notices (spec 15 item 5) and the time estimate."""
        notices = ca.check_notices(self.compat)
        problems = [text for level, text, confirmable in notices if level == ca.CRITICAL and not confirmable]
        confirm = [text for level, text, confirmable in notices if level == ca.CRITICAL and confirmable]
        cautions = [text for level, text, _ in notices if level == ca.CAUTION]
        self._show_notices(cautions, problems + confirm)
        self.run_anyway.setVisible(bool(confirm) and not problems)
        if not problems:
            try:
                edits = [haplotype.Edit.from_event(e) for e in self.compat.alignment.events if e.apply]
                _, cmap = haplotype.build(self.compat.ref_seq, self.compat.region[0], edits)
                reader = self.main.refdb.reader(self.build_combo.currentText())
                in_db = reader is not None and reader.plan(self.gene.name) is not None
                start, end = self.gene.mask[0] - 50, self.gene.mask[1] + 50
                windows = haplotype.recompute_windows(cmap, *(cmap.lift(start, 1), cmap.lift(end, -1)))
                seconds = haplotype.estimate_seconds(windows) + (0 if in_db else (end - start) / 38_000) + 15
                self.estimate_label.setText(
                    t("Estimated run time: about {time} ({n} stretch(es) near the edits predicted again).",
                      time=format_elapsed(seconds), n=len(windows)) if in_db else
                    t("Estimated run time: about {time} ({n} stretch(es) near the edits predicted again; the "
                      "reference computed live first).", time=format_elapsed(seconds), n=len(windows)))
            except ValueError as exc:
                problems = [str(exc)]
                self._show_notices(cautions, problems + confirm)
        self._blocking = bool(problems) or (bool(confirm) and not self.run_anyway.isChecked())
        self._confirm = bool(confirm) and not problems
        if problems:
            self.check_status.setText("\u2717 " + t("Run is blocked: {n} critical problem(s) -- the red ! on the "
                                                      "right says which.", n=len(problems)))
            self.check_status.setStyleSheet(theme.STATUS_CRIT)
        elif confirm:
            self.check_status.setText(t("\u2717 A change touches the ATG or the transcript's start or end (the red !): "
                                      "tick Run anyway to go on."))
            self.check_status.setStyleSheet(theme.STATUS_CRIT)
        else:
            self.check_status.setText("\u2713 " + t("Checked.") + (" " + t("{n} caution(s): the ? on the right.",
                                                                          n=len(cautions)) if cautions else ""))
            self.check_status.setStyleSheet(theme.STATUS_OK)
        self._update_enabled()

    def _show_notices(self, cautions, critical):
        for button, texts in ((self.caution_button, cautions), (self.critical_button, critical)):
            button.setToolTip("\n\n".join(t(x) for x in texts))
            button.setVisible(bool(texts))

    # --- run ---------------------------------------------------------------------------------------------
    def _update_enabled(self, *_):
        blocked = getattr(self, "_blocking", True)
        if getattr(self, "_confirm", False):
            blocked = not self.run_anyway.isChecked()
        ready = self.compat is not None and not blocked
        running = self._run_worker is not None
        self.check_button.setEnabled(bool(self.build_combo.currentText() and self.gene is not None
                                          and not getattr(self, "_no_prediction", False)
                                          and self.sequence_edit.toPlainText().strip()) and not running)
        self.run_button.setEnabled(ready and not running)
        self.cancel_button.setEnabled(running)
        self.events_table.setEnabled(not running)
        has = self.result is not None
        for b in (self.export_csv_button, self.download_button, self.bedgraph_button, self.viewer_button,
                  self.specificity_button, self.breakpoints_button):
            b.setEnabled(has and not running)

    def _settings(self):
        return ca.Settings(self.threshold.value(), self.min_exon.value(), self.min_intron.value(), self.top.value(),
                           self.apply_min_size.value(), config.load()["cdna_low_concordance"], self.weak_delta.value())

    def _save_settings(self):
        settings = config.load()
        settings.update({"cdna_threshold": self.threshold.value(), "cdna_min_exon": self.min_exon.value(),
                         "cdna_min_intron": self.min_intron.value(), "cdna_apply_min_size": self.apply_min_size.value(),
                         "cdna_top": self.top.value(), "cdna_weak_delta": self.weak_delta.value(),
                         "cdna_specificity_n": self.specificity_n.value()})
        config.save(settings)

    def _on_run(self):
        build, gene, compat, settings = self.build_combo.currentText(), self.gene, self.compat, self._settings()
        self._run_request, self._pending_vcf = self._pending_vcf, None     # a run sent from the VCF tab, if this is
        self._save_settings()
        try:
            scores, table, record = self._scores(build), self._table(build), self._fasta(build)["chr" + gene.chrom]
        except ValueError as exc:
            QMessageBox.warning(self, t("Run"), t(str(exc)))
            return
        fetch = lambda a, b: record[a - 1:b].seq
        worker = RunWorker(lambda checkpoint, progress: ca.run_analysis(build, gene, compat, scores, fetch, table,
                                                                        settings, checkpoint, progress), self)
        worker.progress.connect(self._on_progress)
        worker.finished_ok.connect(self._on_run_done)
        worker.failed.connect(self._on_run_failed)
        worker.cancelled.connect(self._on_run_cancelled)
        self._run_worker = worker
        self._run_started, self._fraction, self._step = time.time(), 0.0, N_("Starting")
        self.progress_bar.setValue(0)
        self._timer.start()
        worker.start()
        self._update_enabled()

    def _run_extra(self, step, work, done):
        """An optional run on the results shown (spec 14 C3, C6), with the run's
        progress and Cancel: work(result, scores, checkpoint, on_progress(done,
        total)) in the background, then done(its output)."""
        result = self.result
        self._save_settings()
        try:
            scores = self._scores(result.build)
        except ValueError as exc:
            QMessageBox.warning(self, t(step), t(str(exc)))
            return
        worker = RunWorker(lambda checkpoint, progress: work(result, scores, checkpoint,
                                                             lambda d, t: progress(d / max(t, 1), step)), self)
        worker.progress.connect(self._on_progress)
        worker.finished_ok.connect(lambda out: self._extra_done(step, done, out))
        worker.failed.connect(self._on_run_failed)
        worker.cancelled.connect(self._on_run_cancelled)
        self._run_worker = worker
        self._run_started, self._fraction, self._step = time.time(), 0.0, step
        self.progress_bar.setValue(0)
        self._timer.start()
        worker.start()
        self._update_enabled()

    def _extra_done(self, step, done, out):
        self._finish_run()
        done(out)
        self.progress_label.setText(t("{step} done in {time}.", step=t(step),
                                       time=format_elapsed(time.time() - self._run_started)))

    def _on_specificity(self):
        """Spec 14 C3: the optional sequence-specificity control."""
        n = self.specificity_n.value()
        self._run_extra(N_("Sequence-specificity control"),
                        lambda r, s, c, p: specificity.run(r, s, n, checkpoint=c, on_progress=p),
                        lambda out: (setattr(self.result, "specificity", out), self._show_extras()))

    def _on_breakpoints(self):
        """Spec 14 C6: the optional breakpoint sensitivity."""
        self._run_extra(N_("Breakpoint sensitivity"),
                        lambda r, s, c, p: specificity.breakpoint_sensitivity(r, s, checkpoint=c, on_progress=p),
                        lambda out: (setattr(self.result, "breakpoints", out), self._show_extras()))

    def _show_extras(self):
        r = self.result
        controls, points = (r.specificity, r.breakpoints) if r else ([], [])
        self.specificity_label.setText(t("Sequence-specificity control (an empirical comparison, not a p-value): ")
                                       + " ".join(t(c.text()) + "." for c in controls) if controls else "")
        self.breakpoints_label.setText(t("Breakpoint sensitivity (each breakpoint moved by 10 and 50 bases either "
                                         "way): ") + " ".join(t(b.text()) + "." for b in points) if points else "")

    def _on_progress(self, fraction, step):
        self._fraction, self._step = fraction, step
        self.progress_bar.setValue(int(1000 * fraction))
        self._update_status()

    def _finish_run(self):
        self._run_worker = None
        self._timer.stop()
        self._update_status()
        self._update_enabled()

    def _on_run_done(self, result):
        self._finish_run()
        self.progress_label.setText(t("Done in {time}.", time=format_elapsed(result.seconds))
                                    + ("" if result.wt_source == "database" else
                                       " " + t("The reference values were computed live (this gene isn't in the "
                                               "reference score database yet); the results are the same, only "
                                               "slower.")))
        request, self._run_request = self._run_request, None
        if request is not None:
            result.from_vcf = (request, request.compare(result))
        self.show_result(result)
        if request is not None and self._vcf_folder:
            name = f"{request.gene.name}_{request.build}_{request.row.chrom}-{request.row.pos}-" \
                   f"{request.row.ref[:12]}-{request.row.alt[:12]}"
            self._vcf_report.append(f"{request.name} ({request.gene.name}): "
                                    + os.path.basename(self.save_results(os.path.join(self._vcf_folder, name))))
            QTimer.singleShot(0, self._next_vcf)

    def _on_run_failed(self, message):
        self._finish_run()
        self.progress_label.setText(t("Failed."))
        if self._run_request is not None and self._vcf_folder:
            self._vcf_report.append(f"{self._run_request.name}: failed -- {message}")
            self._run_request = None
            QTimer.singleShot(0, self._next_vcf)
            return
        self._run_request = None
        QMessageBox.critical(self, t("Run failed"), t(message))

    def _on_run_cancelled(self):
        self._finish_run()
        self.progress_bar.setValue(0)
        self.progress_label.setText(t("Cancelled.") + (" " + t("The {n} variant(s) still queued are not run.",
                                                               n=len(self._vcf_queue)) if self._vcf_queue else ""))
        self._vcf_queue, self._run_request = [], None

    # --- from the VCF tab (spec 16 item 4) ---------------------------------------------------------------
    def start_from_vcf(self, requests, folder=None):
        """Runs each vcf_bridge.Request in turn: the build and gene filled in, the variant(s) entered in "From
        coordinates", the sequence generated, checked and run. With a folder (a queue), each result is saved
        there as its own results file."""
        self._vcf_queue, self._vcf_folder, self._vcf_report = list(requests), folder, []
        if self.missing_panel.isVisible():
            self._on_continue_without_database()        # the dialog said the reference is computed live
        self._next_vcf()

    def _next_vcf(self):
        if not self._vcf_queue:
            if self._vcf_folder and self._vcf_report:
                QMessageBox.information(self, t("Analyse transcript effect"),
                                        t("{n} variant(s), results in {folder}:", n=len(self._vcf_report),
                                          folder=self._vcf_folder) + "\n\n" + "\n".join(self._vcf_report))
            self._vcf_folder = None
            return
        req = self._vcf_queue.pop(0)
        for widget, setter, value in ((self.build_combo, self.build_combo.setCurrentText, req.build),
                                      (self.gene_edit, self.gene_edit.setText, req.gene.name),
                                      (self.window_edit, self.window_edit.setText, "")):
            widget.blockSignals(True)           # filled in as a whole below, not looked up again
            setter(value)
            widget.blockSignals(False)
        self._table(req.build)                  # loaded already by the dialog
        self._genes_found(req.build, req.gene.name, [req.gene])
        self.input_tabs.setCurrentIndex(1)
        self.coords_edit.setPlainText("\n".join(req.lines))
        self._pending_vcf, self._auto_run = req, True
        self.progress_label.setText(t("From the VCF tab: {variant} in {gene}", variant=req.name, gene=req.gene.name)
                                    + (" " + t("({i} of {n})", i=len(self._vcf_report) + 1,
                                               n=len(self._vcf_report) + 1 + len(self._vcf_queue))
                                       if self._vcf_folder else "") + ".")
        self._on_generate()

    def _auto_start(self):
        self._auto_run = False
        if self.run_button.isEnabled():
            self._on_run()
            return
        self._vcf_stopped(self.check_status.text() + " " + self.critical_button.toolTip())

    def _vcf_stopped(self, why):
        """A variant sent from the VCF tab that can't be run as it is: in a queue, skipped; else said."""
        self._auto_run = False
        req, self._pending_vcf = self._pending_vcf, None
        if self._vcf_folder:
            self._vcf_report.append(f"{req.name if req else '?'}: not run -- {why.strip()}")
            QTimer.singleShot(0, self._next_vcf)
        else:
            QMessageBox.warning(self, t("Analyse transcript effect"), t("{variant} isn't run: {why}\n\nThe check's notices are next to Check compatibility; "
                                  "once they allow it, Run starts it.", variant=req.name if req else t("The variant"),
                                  why=why.strip()))

    def _show_vcf_header(self, result):
        found = getattr(result, "from_vcf", None)
        self.vcf_header.setVisible(found is not None)
        if found is None:
            return
        req, comparison = found
        vcf_line, engine_line = comparison.lines(req.mode)
        others = ("<br>" + t("Applied together with {n} other variant(s) of the VCF on one haplotype (their phase is "
                             "unknown): ", n=len(req.others))
                  + ", ".join(f"{r.chrom}:{r.pos:,} {r.ref}>{r.alt}" for r in req.others) if req.others else "")
        verdict = ("✓ " + t("within the rounding tolerance (±0.005)") if comparison.agree else
                   t("differs beyond the rounding tolerance (the ? says why): ") + "; ".join(t(p) for p in comparison.problems))
        self.vcf_label.setText(
            t("<b>From the VCF tab</b>: {variant} in {gene}; genotype {genotype}", variant=html.escape(req.name),
              gene=html.escape(req.gene.name), genotype=html.escape(t(req.genotype.text()))) + f"{others}<br>"
            f"{html.escape(t(vcf_line))}<br>"
            f"{html.escape(t(engine_line))} — {html.escape(verdict)}")
        self.vcf_why.setToolTip("\n\n".join(t(r) for r in comparison.reasons))
        self.vcf_why.setVisible(not comparison.agree)

    def _on_cancel(self):
        if self._run_worker is not None:
            self._run_worker.stop()
            self.progress_label.setText(t("Cancelling..."))

    def _update_status(self):
        parts = []
        if self._run_worker is not None and self._run_started:
            elapsed = time.time() - self._run_started
            self.progress_label.setText(f"{t(self._step)} ({100 * self._fraction:.0f}%)")
            eta = t(", about {time} left", time=format_elapsed(elapsed / self._fraction - elapsed)) \
                if self._fraction > 0.05 else ""
            parts.append(t("Elapsed {time}", time=format_elapsed(elapsed)) + eta)
        _cpu, ram = self._meter.sample()
        if ram:
            parts.append(f"RAM {ram / 1e9:.1f} GB")
        build = self.build_combo.currentText()
        m = self.main.refdb
        if build:
            s = m.status.get(build, {})
            db = {"ready": "OK", "partial": N_("partly built"), "building": N_("being built"), "missing": N_("missing"),
                  "refused": N_("not usable"), "checking": N_("checking")}.get(s.get("state"), "?")
            version = (s.get("manifest") or {}).get("software_version")
            parts.append(t("Database {build}: {state}", build=build, state=t(db))
                         + (" " + t("(version {version})", version=version) if version else ""))
        parts.append(t("cache {mb:,.1f} MB", mb=refdb.folder_size(m.cache_dir()) / 1e6))
        self.status_label.setText(" · ".join(parts))

    # --- results -------------------------------------------------------------------------------------
    def show_result(self, result):
        self.result = result
        self._note(result.reference_note[1])
        self._notices = ca.result_notices(result)       # spec 15 item 5: icons next to the rows they concern
        self.lead_label.setText(t(ca.lead_text(result)))
        self._show_vcf_header(result)
        self.local_label.setText("\n".join(t(d["text"]) for d in ca.local_concordance(result)))
        self.event_lines = [d["line"] for d in ca.event_summary(result)]
        colors = [EVENT_COLORS.get(d["kind"], "#555") for d in ca.event_summary(result)]
        self.events_label.setText((t("<b>Applied events</b> (click one to see it in the viewer):") + "<br>" + "<br>".join(
            f"<span style='color: {c};'>■</span> <a href='event:{i}'>{html.escape(line.split(' · ')[0])}</a> · "
            f"{html.escape(line.split(' · ', 1)[1])}" for i, (line, c) in enumerate(zip(self.event_lines, colors))))
            if self.event_lines else "")
        self._show_extras()
        self._fill_table()
        self._update_enabled()

    def visible_rows(self):
        rows = ca.table_rows(self.result) if self.result else []
        if not self.show_alternatives.isChecked():
            rows = [r for r in rows if not (r["Type"] == "WT predicted" and r["Rank"] not in (1, ""))]
        return rows

    def _fill_table(self, *_):
        rows = self.visible_rows()
        self.results_table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            kind = {"Sample": "Patient"}.get(row["Type"], row["Type"])
            mine = [(level, text) for level, text, where in getattr(self, "_notices", [])
                    if (kind, row["Rank"] or 0) in where]
            for j, col in enumerate(ca.COLUMNS):
                item = QTableWidgetItem(str(row[col]).replace(" -> ", " \u2192 "))
                if col in ("Change vs Reference", "Change vs WT predicted", "Evidence vs WT", "RNA check (RT-PCR)"):
                    item.setToolTip(str(row[col]))
                if j == 0 and mine:
                    level = ca.CRITICAL if any(lv == ca.CRITICAL for lv, _ in mine) else ca.CAUTION
                    item.setIcon(notice_icon(level))
                    item.setToolTip("\n\n".join(text for _, text in mine))
                    item.setData(Qt.UserRole, level)
                self.results_table.setItem(i, j, item)
        self.results_table.resizeColumnsToContents()

    def _on_row_clicked(self, row, _column):
        rows = self.visible_rows()
        if 0 <= row < len(rows):
            r = rows[row]
            self.open_viewer(r["Rank"] if r["Type"] == "Sample" else 1)

    def _on_event_link(self, link):
        viewer = self.open_viewer()
        if viewer is not None:
            viewer._on_event_index(int(link.split(":")[1]))

    def open_viewer(self, candidate=1):
        """The viewer of the results shown, the given sample candidate chosen in its sequence view."""
        from .viewer_window import ViewerWindow

        if self.result is None:
            return None
        if self.viewer is None or getattr(self.viewer.data, "result", None) is not self.result:
            self.viewer = ViewerWindow(self.result, self, candidate=candidate)
        else:
            self.viewer.set_candidate(candidate)
        self.viewer.show()
        self.viewer.raise_()
        return self.viewer

    def _on_export_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, t("Export CSV"), self._default_name() + ".csv",
                                              t("CSV (*.csv);;Tab-separated (*.tsv)"))
        if not path:
            return
        keys = [f"c{i}" for i in range(len(ca.COLUMNS))]
        rows = [SimpleNamespace(**{k: r[c] for k, c in zip(keys, ca.COLUMNS)}) for r in self.visible_rows()]
        # Spec 16: the applied events (with where inserted bases come from) above the table.
        notes = [ca.RESEARCH_NOTE] + [f"Event {i}: {d['line']}" for i, d in enumerate(ca.event_summary(self.result), 1)]
        write_rows_ordered(rows, path, list(zip(ca.COLUMNS, keys)), delimiter=delimiter_for_path(path), note=notes)

    def _default_name(self):
        r = self.result
        return f"{r.gene.name}_{r.build}_sequence-to-cdna" if r else "sequence-to-cdna"

    def _on_download(self):
        path, _ = QFileDialog.getSaveFileName(self, t("Download results"), self._default_name() + bundle.EXTENSION,
                                              t("Sequence-to-cDNA results (*{ext})", ext=bundle.EXTENSION))
        if path:
            self.save_results(path)

    def save_results(self, path):
        from spliceai_gui.spliceai_setup import SPLICEAI_VERSION

        r = self.result
        status = self.main.refdb.status.get(r.build, {})
        manifest = status.get("manifest") or {}
        versions = {"spliceai": SPLICEAI_VERSION, "models": manifest.get("model_sha256"),
                    "database": {k: manifest.get(k) for k in ("software_version", "updated", "scope", "precision")}
                    if manifest else None, "reference_values": r.wt_source}
        if getattr(r, "from_vcf", None):
            versions["from_vcf"] = vcf_bridge.summary(*r.from_vcf)        # spec 16 item 4
        return bundle.save(r, path, versions=versions)

    def _on_open_results(self):
        path, _ = QFileDialog.getOpenFileName(self, t("Open calculated results"), "",
                                              t("Sequence-to-cDNA results (*{ext})", ext=bundle.EXTENSION))
        if path:
            self.open_results(path)

    def open_results(self, path):
        try:
            result, manifest, notes = bundle.load(path)
        except bundle.BundleError as exc:
            QMessageBox.warning(self, t("Open calculated results"), str(exc))
            return None
        self.show_result(result)
        self.progress_label.setText(t("Opened {name}: {gene} on {build}, made {created}.", name=os.path.basename(path),
                                      gene=manifest["gene"], build=manifest["build"], created=manifest["created"])
                                    + (f" {' '.join(notes)}" if notes else ""))
        return result

    def _on_export_bedgraph(self):
        folder = QFileDialog.getExistingDirectory(self, t("Folder for the bedGraph tracks"))
        if folder:
            self.export_bedgraph(os.path.join(folder, self._default_name()))

    def export_bedgraph(self, base):
        """Writes reference, sample and delta tracks (acceptor and donor) as
        <base>.<track>.bedgraph.gz with their .tbi; returns the paths."""
        r = self.result
        plan, chrom = r.plan, "chr" + r.gene.chrom
        patient = haplotype.ref_aligned(r.patient, r.cmap, plan["start"], plan["end"])
        tracks = {"reference_acceptor": r.wt_values[:, 0], "reference_donor": r.wt_values[:, 1],
                  "sample_acceptor": patient[:, 0], "sample_donor": patient[:, 1],
                  "delta_acceptor": patient[:, 0] - r.wt_values[:, 0], "delta_donor": patient[:, 1] - r.wt_values[:, 1]}
        paths = []
        for name, values in tracks.items():
            path = f"{base}.{name}.bedgraph.gz"
            bedgraph.write(path, chrom, plan["start"], values,
                           description=f"SpliceAI {name.replace('_', ' ')} probability, {r.gene.name} {r.build}")
            paths.append(path)
        return paths
