"""Main window: wires the input area, options, background run, and results table
around the existing spliceai_pipeline package. No pipeline logic lives here --
only UI state and calls into spliceai_pipeline.cli.run_pipeline_core (via
worker.PipelineWorker) and spliceai_pipeline.writer.write_rows_ordered (for export).
"""
import concurrent.futures
import gzip
import os
import tempfile
import time

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog, QFileDialog, QFrame,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMenu, QMessageBox, QPlainTextEdit,
    QProgressBar, QPushButton, QRadioButton, QScrollArea, QSplitter, QStyle, QTableView, QTabWidget, QToolButton,
    QVBoxLayout, QWidget,
)

from spliceai_pipeline import haplotype, refdb, vcf_bridge
from spliceai_pipeline.cli import BUILDS, MODES, format_elapsed
from spliceai_pipeline.mane import (
    DEFAULT_MANE_DIR, find_cached_summary, is_summary_filename, sibling_mane_dir,
)
from spliceai_pipeline.normalize import normalize_variant
from spliceai_pipeline.snpeff import (
    DEFAULT_SNPEFF_DIR, JavaNotFoundError, SnpEffNotSetUpError, check_snpeff_setup, find_java, snpeff_jar_path,
)
from spliceai_pipeline.vcfio import strip_vcf_extension
from spliceai_pipeline.writer import ResultsFileError, read_rows, write_rows_ordered

from . import config
from . import spliceai_setup
from . import theme
from .assets_paths import LOGO_PNG
from .cdna_tab import CdnaTab
from .cpu_meter import CpuMeter
from .download_dialog import ReferenceDownloadDialog
from .gene_import_dialog import GeneImportDialog
from .help_button import make_help_button
from .mane_download_dialog import ManeDownloadDialog
from . import reference_locator
from .i18n import LANGUAGES, language, t
from .refdb_wizard import RefdbDialog, RefdbManager, database_folder_of, describe as describe_refdb, live_seconds
from .snpeff_download_dialog import SnpEffDownloadDialog
from .spliceai_setup_dialog import SpliceAISetupDialog
from .score_delegate import THRESHOLDS as SCORE_THRESHOLDS
from .score_delegate import ScoreBarDelegate, SeverityStripeDelegate
from .table_model import COLUMNS, MAX_SCORE_COLUMN, ScoreFilterProxyModel, ScoreTableModel, sort_key_value
from .transcript_effect_dialog import TranscriptEffectDialog
from .vcf_check import check_vcf_text
from .vcf_info import detect_vcf_build, summarize_vcf_text
from .worker import PipelineWorker

# Shown in front of the progress text while a run is paused.
PAUSED_PREFIX = t("⏸ Paused -- ")

# The filter buttons above the results and the colour bands in the table are
# the same three SpliceAI cutoffs, taken from one place so they can't drift
# apart (see score_delegate.THRESHOLDS for what each one means).
THRESHOLDS = tuple(sorted(lower for lower, _colour in SCORE_THRESHOLDS))


# Rows per page in the results table (see _refresh_page, _on_prev_page and
# _on_next_page).
#
RESULTS_PAGE_SIZE = 100








MODE_HELP_TEXT = (
    t("raw: the model's unfiltered delta scores.\n"
      "\n"
      "masked: predictions that do not match an annotated splice boundary in the expected direction are set to "
      "zero, to suppress biologically implausible predictions. For the same variant, masked scores are never "
      "higher than raw scores.")
)

SCORE_HELP_TEXT = (
    t("Column abbreviations:\n"
      "  • DS = Delta Score: the predicted probability (0-1) that the variant creates or disrupts a splice site of "
      "the given type\n"
      "  • DP = Delta Position: the distance in base pairs between the variant and the predicted splicing change "
      "(negative or positive, indicating the direction)\n"
      "  • AG = Acceptor Gain: a new splice acceptor site is predicted\n"
      "  • AL = Acceptor Loss: the loss of an existing splice acceptor site is predicted\n"
      "  • DG = Donor Gain: a new splice donor site is predicted\n"
      "  • DL = Donor Loss: the loss of an existing splice donor site is predicted\n"
      "\n"
      "Together these form four score/position column pairs: DS_AG / DP_AG, DS_AL / DP_AL, DS_DG / DP_DG and DS_DL "
      "/ DP_DL.\n"
      "\n"
      "max_score is the largest of the four scores and summarises the overall likelihood that the variant affects "
      "splicing.\n"
      "\n"
      "Commonly used interpretation thresholds (conventions, not diagnostic criteria): the original SpliceAI "
      "publication (Jaganathan et al., Cell 2019) and many laboratories use about 0.2 as a sensitive cut-off, "
      "about 0.5 as a standard threshold and about 0.8 as a high-confidence threshold.\n"
      "\n"
      "A high score indicates a predicted splicing effect; it does not by itself establish clinical pathogenicity. "
      "Variant interpretation also requires the ACMG criteria, clinical correlation and other supporting "
      "evidence.\n"
      "\n"
      "Allele Fraction is read from the VCF's AF field (the per-sample FORMAT value is preferred to the INFO "
      "value). If there is no AF, it is calculated from AD as alt / (ref + alt); if neither is present, the column "
      "is left blank. It is not needed for scoring.")
)

PRECOMPUTED_HELP_TEXT = (
    t('The SpliceAI authors (Jaganathan et al., Cell 2019) have precomputed scores for every possible singl'
      'e-nucleotide variant and for common small indels (1 bp insertions, 1-4 bp deletions). These files ar'
      'e optional and can be downloaded from Illumina BaseSpace (basespace.illumina.com/s/otSPW8hnhaZR) und'
      "er Illumina's terms. Scores are available for hg19 and hg38, in raw and masked mode; for one build, "
      'the masked SNV file is about 29 GB and the masked indel file about 69 GB.\n\nWhen the program can read'
      ' these files, variants already in them are not scored again: their scores are taken from the files, '
      'which shortens the run. Without them ("I don\'t have precomputed data") every variant is scored live;'
      ' this takes longer, but the results are the same.\n\nFile names: spliceai_scores.{raw,masked}.{snv,ind'
      'el}.hg{19,38}.vcf.gz, each with its .tbi index. Keep the files compressed as downloaded (bgzip, with'
      ' the .tbi index); do not decompress them. Then select "Use precomputed data" and choose the folder t'
      'hat contains them.')
)

SCORE_GENES_HELP_TEXT = (
    t("Scores only the variants within the selected genes, so a whole genome or exome takes minutes instead of "
      "hours. A variant is considered to lie within a gene when it falls in SpliceAI's interval for that gene (the "
      "same rule that assigns a score to a gene), and only rows for those genes are kept.\n"
      "\n"
      "Run is enabled once one of the two is chosen: All genes scores every variant; Choose gene list... accepts "
      "a pasted list (for example an HPO phenotype list or a gene panel) and can save it under a name for later "
      "use. The last list is offered again.\n"
      "\n"
      "This differs from the gene-set filter above the results table, which only hides rows of a run in which all "
      "variants were scored.")
)

COMPAT_HELP_TEXT = (
    t("\"Check compatibility\" performs a quick structural check of the VCF text currently loaded or pasted, "
      "without running the scoring, so that problems can be found before a potentially long run.\n"
      "\n"
      "It verifies that:\n"
      "  • a #CHROM header line with the required columns (#CHROM, POS, ID, REF, ALT) is present\n"
      "  • at least one variant record follows the header\n"
      "  • a sample of records can be parsed: POS is a positive integer, REF contains valid bases and ALT is not "
      "empty\n"
      "\n"
      "It also reports the number of variant records and whether the chromosome names carry the \"chr\" prefix.\n"
      "\n"
      "Typical variant counts (these vary considerably between laboratories, capture kits and filtering):\n"
      "  • WES (whole exome): tens of thousands\n"
      "  • WGS (whole genome): several million\n"
      "  • CES / clinical exome or gene panels: from hundreds to a few thousand")
)

FASTA_HELP_TEXT = t(
    "There is one field per genome build. A run uses the field for the build selected above, which is detected "
    "automatically from the VCF header when the header states it (##reference, or the length of chr1 in its "
    "##contig lines).\n"
    "\n"
    "The simplest option is \"Download...\" next to a field: it downloads the build's FASTA from UCSC (about 1 "
    "GB compressed, about 3 GB decompressed), decompresses it and creates its .fai index.\n"
    "\n"
    "Alternatively, download hg19.fa.gz or hg38.fa.gz from hgdownload.soe.ucsc.edu/goldenPath/{build}/bigZips/, "
    "decompress it, create a .fai index (for example with `samtools faidx`), and select the .fa file with "
    "Browse.\n"
    "\n"
    "The reference FASTA is required even when precomputed data are selected, because variants not covered by "
    "them are scored live using the reference sequence."
)

SNPEFF_HELP_TEXT = (
    t("SnpEff adds transcript and coding-position annotation (see \"About SnpEff\" (?)). SnpEff "
      "is a separate Java program, not part of this Python program.\n"
      "\n"
      "One-time setup, if needed:\n"
      "  \"Download SnpEff...\" below downloads SnpEff and the RefSeq transcript database for the selected build "
      "(GRCh37.p13/GRCh38.p14, which provide NM_ identifiers), with a progress bar as for the reference FASTA.\n"
      "  SnpEff requires Java 21 or later. If this computer does not have it, the same download also obtains a "
      "private copy of Java (Eclipse Temurin, about 50 MB) in the program's \"java\" folder, used only by this "
      "program; any Java already installed is left unchanged.\n"
      "\n"
      "Once set up, Run fills in two additional columns.")
)

SNPEFF_ABOUT_HELP_TEXT = (
    t("SnpEff is an open-source variant annotation tool that determines the gene and transcript in which a variant "
      "lies and converts its genomic position into coding-sequence notation (e.g. NM_018136.5, c.875A>T).\n"
      "\n"
      "When a variant overlaps more than one transcript, the transcript shown is chosen in this order: the gene's "
      "NCBI MANE Select transcript, if available (matched by RefSeq accession, independent of the genome build); "
      "otherwise a curated RefSeq transcript (NM_/NR_) in preference to a computationally predicted one (XM_/XR_); "
      "otherwise the first transcript listed by SnpEff.\n"
      "\n"
      "SpliceAI predicts effects on splicing but does not provide transcript or coding notation; SnpEff "
      "complements it, so that each row shows both the splicing prediction and the standard transcript and coding "
      "notation.\n"
      "\n"
      "SnpEff is open source under the MIT license and was developed by Pablo Cingolani.")
)

MANE_HELP_TEXT = (
    t("NCBI's MANE Select is a single, expert-curated representative transcript per gene, agreed by NCBI and "
      "Ensembl. When a variant overlaps more than one transcript, SnpEff's transcript and coding-position columns "
      "report this transcript in preference to the others.\n"
      "\n"
      "It is optional: without it, a curated RefSeq transcript (NM_/NR_) is preferred to a predicted one "
      "(XM_/XR_), and otherwise the first transcript listed by SnpEff is used.\n"
      "\n"
      "The file is small (about 1 MB) and is downloaded from NCBI's FTP site; the latest release is located "
      "automatically.")
)


def _is_gzip(path):
    """Sniff the actual magic bytes rather than trusting the filename -- real-world
    VCFs are sometimes gzip/BGZF-compressed without a .gz/.bgz suffix."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _read_vcf_text(path):
    """Read a VCF file as text, auto-detecting gzip/BGZF by content (not filename).

    Raises RuntimeError with a message identifying which step of the
    load/decompress/decode chain failed, instead of a bare low-level exception.
    """
    try:
        is_gz = _is_gzip(path)
    except OSError as exc:
        raise RuntimeError(f"could not read file to detect compression: {exc}") from exc

    try:
        opener = gzip.open if is_gz else open
        # utf-8-sig: drops a byte-order mark, which would otherwise sit in front
        # of the first header line (see spliceai_pipeline/vcfio.py).
        with opener(path, "rt", encoding="utf-8-sig") as fh:
            return fh.read()
    except gzip.BadGzipFile as exc:
        raise RuntimeError(f"file looked gzip-compressed (magic bytes 1f 8b) but failed to decompress: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise RuntimeError(f"decompressed to bytes that aren't valid UTF-8 text: {exc}") from exc
    except OSError as exc:
        raise RuntimeError(f"could not open/read file: {exc}") from exc


# File names treated as VCFs (by "Load file..." and drag and drop) and as
# reference FASTAs when dropped onto the window.
VCF_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz", ".bgz")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna")

# Starting heights of the three panes (results + run controls, Input VCF,
# Settings) the first time the window opens; after that whatever the user
# dragged them to is remembered. Results get the room, the VCF box only as
# much as it needs to show which file was loaded.
DEFAULT_PANE_SIZES = [520, 150, 260]


def dropped_local_files(mime_data):
    """Local paths in a drag (e.g. files from Explorer), in order; [] for a
    text drag."""
    if not mime_data.hasUrls():
        return []
    # normpath: Qt gives C:/Users/..., shown in the fields as C:\Users\...
    return [os.path.normpath(url.toLocalFile()) for url in mime_data.urls() if url.isLocalFile()]


def dropped_vcf_path(mime_data):
    """The first dropped file with a VCF name, or None."""
    return next((p for p in dropped_local_files(mime_data) if p.lower().endswith(VCF_SUFFIXES)), None)


class VcfDropTextEdit(QPlainTextEdit):
    """QPlainTextEdit that also accepts a VCF file dropped onto it -- normal
    text drag-drop (e.g. dragging a text selection) still falls through to the
    base class untouched. Other files are left to the main window (see
    MainWindow.dropEvent) instead of being pasted in as a file:/// URL."""

    fileDropped = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptDrops(True)
        self._base_style = ""
        self._drag_style = "QPlainTextEdit { border: 2px dashed #1a5fb4; }"

    def set_drop_highlight(self, on):
        self.setStyleSheet(self._drag_style if on else self._base_style)

    def dragEnterEvent(self, event):
        if dropped_vcf_path(event.mimeData()):
            event.acceptProposedAction()
            self.set_drop_highlight(True)
        elif event.mimeData().hasUrls():
            event.ignore()
        else:
            super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if dropped_vcf_path(event.mimeData()):
            event.acceptProposedAction()
        elif event.mimeData().hasUrls():
            event.ignore()
        else:
            super().dragMoveEvent(event)

    def dragLeaveEvent(self, event):
        self.set_drop_highlight(False)
        super().dragLeaveEvent(event)

    def dropEvent(self, event):
        self.set_drop_highlight(False)
        path = dropped_vcf_path(event.mimeData())
        if path:
            event.acceptProposedAction()
            self.fileDropped.emit(path)
        elif event.mimeData().hasUrls():
            event.ignore()
        else:
            super().dropEvent(event)


class PathLineEdit(QLineEdit):
    """QLineEdit for a file or folder path that also takes one dragged in from
    Explorer: it emits pathDropped(local_path) -- MainWindow decides what the
    file is for -- instead of inserting its file:/// URL as text, which is
    what a plain QLineEdit does. Typing and pasting work as usual."""

    pathDropped = Signal(str)

    def dragEnterEvent(self, event):
        if dropped_local_files(event.mimeData()):
            event.acceptProposedAction()
        else:
            super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if dropped_local_files(event.mimeData()):
            event.acceptProposedAction()
        else:
            super().dragMoveEvent(event)

    def dropEvent(self, event):
        files = dropped_local_files(event.mimeData())
        if files:
            event.acceptProposedAction()
            self.pathDropped.emit(files[0])
        else:
            super().dropEvent(event)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(t("SpliceAI Variant Scoring"))
        self.resize(1200, 800)
        # Start no bigger than the screen; on a small one everything scrolls
        # (see _build_ui).
        screen = QApplication.primaryScreen()
        if screen is not None:
            available = screen.availableGeometry()
            self.resize(min(1200, available.width() - 40), min(800, available.height() - 60))
        if LOGO_PNG.exists():
            self.setWindowIcon(QIcon(str(LOGO_PNG)))

        # The look is one stylesheet on the QApplication, so it reaches every
        # widget including the dialogs (theme.py). main() applies it before
        # anything is built; this is the safety net for a window constructed
        # directly, as the checks in tests/ do -- and it must not overwrite a
        # stylesheet that is already there.
        app = QApplication.instance()
        if app is not None and not app.styleSheet():
            theme.apply_theme(app)

        self._worker = None
        self._temp_vcf_path = None
        self._loaded_vcf_basename = None
        self._all_rows = []
        self._gene_set_filter = None
        self._quick_genes = []
        self._quick_gene_chip_widgets = {}
        self._current_page = 0

        self._run_start_time = None
        self._last_run_settings = None
        self._last_run_duration = None

        # The reference score database and its builds (refdb_wizard.py): its
        # status shows in the top bar; a build left unfinished continues now.
        self.refdb = RefdbManager(lambda build: self.fasta_edits[build].text().strip(), self._mane_dir, self)
        self._refdb_notices = []

        self._build_ui()
        # Files can be dragged onto the whole window, not just the input box
        # (see dropEvent).
        self.setAcceptDrops(True)
        self._load_settings()
        self._init_advanced_section()
        self._update_precomputed_enabled()
        self._update_run_enabled()
        self.refdb.changed.connect(self._update_refdb_indicator)
        self.refdb.notice.connect(self._show_refdb_notice)
        self.refdb.refresh()
        self.refdb.start_next()

    def _build_ui(self):
        # Everything sits in one scroll area, so the window still works on a
        # small screen (or a small window): it scrolls instead of squeezing.
        central = QWidget()
        root = QVBoxLayout(central)

        # Top to bottom: what you came for (the results, with the run controls
        # under them), then the VCF, then the settings. The order a user's eye
        # needs them in -- the settings are read once and then left alone, so
        # they sit at the bottom. Input has to be built before the settings:
        # _update_run_enabled, which the settings' widgets call as they are set
        # up, reads the VCF box.
        #
        # The three are panes of a splitter, so the line between any two can be
        # dragged to give one of them more room, the way Explorer's panes work.
        # Where they are left is remembered (config key "pane_sizes").
        root.addLayout(self._build_top_bar())

        results_pane = QWidget()
        results_layout = QVBoxLayout(results_pane)
        results_layout.setContentsMargins(0, 0, 0, 0)
        results_layout.addWidget(self._build_results_group(), stretch=1)
        # The run controls belong with the results and are never resized on
        # their own, so they ride along at the bottom of this pane.
        results_layout.addLayout(self._build_run_row())

        self.panes = QSplitter(Qt.Vertical)
        self.panes.setChildrenCollapsible(False)
        self.panes.setHandleWidth(7)
        self.panes.addWidget(results_pane)
        self.panes.addWidget(self._build_input_group())
        self.panes.addWidget(self._build_options_group())
        # Extra height goes to the results; the other two keep what they have.
        self.panes.setStretchFactor(0, 1)
        self.panes.setStretchFactor(1, 0)
        self.panes.setStretchFactor(2, 0)
        self.panes.setSizes(DEFAULT_PANE_SIZES)
        # Two tabs: VCF scoring, exactly the page it always was, and Sequence to
        # cDNA (cdna_tab.py). The top bar above them belongs to both.
        self.tabs = QTabWidget()
        self.tabs.addTab(self.panes, t("VCF scoring"))
        self.cdna_tab = CdnaTab(self)
        self.tabs.addTab(self.cdna_tab, t("Sequence to cDNA"))
        root.addWidget(self.tabs, stretch=1)

        # The results table keeps a usable height; when space runs out the
        # page scrolls rather than shrinking it to a few rows.
        self.table_view.setMinimumHeight(260)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setWidget(central)
        self.setCentralWidget(scroll)

    def _make_help_button(self, dialog_title, help_text, font_size="10pt", logo=False):
        """A "(?)" icon that pops up an explanation (help_button.py, shared with the other tab)."""
        return make_help_button(self, dialog_title, help_text, font_size, logo)

    def _build_top_bar(self):
        row = QHBoxLayout()
        title = QLabel(t("SpliceAI Variant Scoring"))
        title.setStyleSheet(f"font-size: 16pt; font-weight: 600; color: {theme.INK};")
        row.addWidget(title)
        row.addStretch(1)

        # Status button for the SpliceAI install (see _update_spliceai_status_indicator
        # and _on_spliceai_status_clicked).
        self.spliceai_status_button = QPushButton()
        self.spliceai_status_button.setStyleSheet(theme.pill("crit"))
        self.spliceai_status_button.setCursor(Qt.PointingHandCursor)
        self.spliceai_status_button.clicked.connect(self._on_spliceai_status_clicked)
        row.addWidget(self.spliceai_status_button)
        self._update_spliceai_status_indicator()
        self.refdb_status_button = QPushButton()
        self.refdb_status_button.setCursor(Qt.PointingHandCursor)
        self.refdb_status_button.clicked.connect(lambda: self._open_refdb_dialog())
        row.addWidget(self.refdb_status_button)
        self._update_refdb_indicator()
        # English / Türkçe (i18n.py): takes effect at the next start, so the notice says it in both.
        self.language_combo = QComboBox()
        for code, name in LANGUAGES.items():
            self.language_combo.addItem(name, code)
        self.language_combo.setCurrentIndex(list(LANGUAGES).index(language()))
        self.language_combo.setToolTip("Language / Dil")
        self.language_combo.currentIndexChanged.connect(self._on_language_changed)
        row.addWidget(self.language_combo)
        return row

    def _on_language_changed(self):
        settings = config.load()
        settings["language"] = self.language_combo.currentData()
        config.save(settings)
        QMessageBox.information(self, "Language / Dil", "The language changes the next time the program starts.\n\n"
                                "Dil, program bir sonraki açılışında değişir.")

    def _update_spliceai_status_indicator(self):
        """Always shown: green when SpliceAI is found (click to check it), red
        when it isn't (click to install it)."""
        if spliceai_setup.is_spliceai_installed():
            text, kind = t("SpliceAI ready"), "ok"
            tip = t("SpliceAI is installed. Click to see its version and location, and whether all its files are there.")
        else:
            text, kind = t("SpliceAI not installed"), "crit"
            tip = t("Live scoring can't run without it. Click to install SpliceAI.")
        self.spliceai_status_button.setText(text)
        self.spliceai_status_button.setToolTip(tip)
        self.spliceai_status_button.setStyleSheet(theme.pill(kind))
        self.spliceai_status_button.setVisible(True)

    def _update_refdb_indicator(self):
        """The reference score database's state, for the Sequence-to-cDNA work
        (the VCF tab doesn't need it): building, ready for which builds, a
        problem, or none yet. Click for details, the build and import."""
        m = self.refdb
        usable = [b for b in m.status if m.usable(b)]
        if m.proc is not None:
            p = m.progress
            pct = f" {100 * p['bases_done'] / p['bases_total']:.0f}%" if p.get("bases_total") else ""
            text, kind = t("Building reference DB") + pct + (t(" (paused)") if p.get("state") == "paused" else ""), "warn"
        elif any(s["state"] == "refused" for s in m.status.values()) or m.errors:
            text, kind = t("Reference DB: problem"), "warn"
        elif usable:
            text, kind = t("Reference DB: {builds}", builds=", ".join(usable)), "ok"
        elif any(s["state"] == "checking" for s in m.status.values()):
            text, kind = t("Reference DB: checking..."), "muted"
        else:
            text, kind = t("No reference DB"), "muted"
        self.refdb_status_button.setText(text)
        self.refdb_status_button.setStyleSheet(theme.pill(kind))
        self.refdb_status_button.setToolTip(
            "\n".join(describe_refdb(s) for s in m.status.values()) + "\n\n"
            + t("Click for details, to build it or import one."))

    def _open_refdb_dialog(self, import_folder=None):
        dialog = RefdbDialog(self.refdb, self)
        if import_folder:
            dialog.import_folder(import_folder)
        dialog.exec()

    def _show_refdb_notice(self, title, text):
        # Not modal: a build can end while the user is busy with something else.
        box = QMessageBox(QMessageBox.Information, title, text, QMessageBox.Ok, self)
        box.setModal(False)
        box.show()
        self._refdb_notices.append(box)
        box.finished.connect(lambda _r, b=box: self._refdb_notices.remove(b))

    def _on_spliceai_status_clicked(self):
        if spliceai_setup.is_spliceai_installed():
            self._show_spliceai_check()
        else:
            SpliceAISetupDialog(self).exec()
        self._update_spliceai_status_indicator()

    def _show_spliceai_check(self):
        info = spliceai_setup.check_installation()
        if info is None:
            SpliceAISetupDialog(self).exec()
            return
        models_found = info["n_models"] - len(info["missing_models"])
        annotations_found = info["n_annotations"] - len(info["missing_annotations"])
        lines = [
            f"SpliceAI {info['version'] or t('(version unknown)')} -- {t(info['origin'])}.",
            "",
            t("Location: {path}", path=info["location"]),
            t("Model files (spliceai1-5.h5): {found} of {n} found", found=models_found, n=info["n_models"]),
            t("Gene annotation files (grch37, grch38): {found} of {n} found", found=annotations_found,
              n=info["n_annotations"]),
            "",
        ]
        box = QMessageBox(self)
        box.setWindowTitle(t("SpliceAI check"))
        missing = info["missing_models"] + info["missing_annotations"]
        if missing:
            box.setIcon(QMessageBox.Warning)
            lines.append("⚠ " + t("Missing: {files}. Live scoring won't work until SpliceAI is reinstalled.",
                                  files=", ".join(missing)))
            reinstall = box.addButton(t("Open SpliceAI setup..."), QMessageBox.ActionRole)
        else:
            box.setIcon(QMessageBox.Information)
            lines.append("✓ " + t("Everything live scoring needs is in place."))
            if info["version"] and info["version"] != spliceai_setup.SPLICEAI_VERSION:
                lines.append(t("(This program is tested with SpliceAI {version}.)",
                               version=spliceai_setup.SPLICEAI_VERSION))
            reinstall = None
        box.setText("\n".join(lines))
        box.addButton(QMessageBox.Ok)
        box.exec()
        if reinstall is not None and box.clickedButton() is reinstall:
            SpliceAISetupDialog(self).exec()

    def _build_input_group(self):
        group = QGroupBox(t("Input VCF"))
        layout = QVBoxLayout(group)

        self.vcf_text = VcfDropTextEdit()
        self.vcf_text.setPlaceholderText(
            t('Paste VCF content here, use "Load file..." below, or drag a .vcf / .vcf.gz file anywhere onto this window.')
        )
        # A few lines is enough: a loaded VCF is read from its file, not from
        # what is on show here, so this box only has to prove the right file
        # arrived. Drag the line below the group to make it taller.
        self.vcf_text.setMinimumHeight(58)
        self.vcf_text.textChanged.connect(self._update_run_enabled)
        self.vcf_text.textChanged.connect(self._on_vcf_text_changed)
        self.vcf_text.fileDropped.connect(self._on_file_dropped)
        layout.addWidget(self.vcf_text)

        row = QHBoxLayout()
        load_btn = QPushButton(t("Load file..."))
        load_btn.setObjectName("stepButton")
        load_btn.setToolTip(t("Open a VCF file -- or drag one from Explorer anywhere onto this window."))
        load_btn.clicked.connect(self._on_load_file)
        row.addWidget(load_btn)
        self.load_status_label = QLabel("")
        row.addWidget(self.load_status_label)

        self.about_vcf_button = QToolButton()
        self.about_vcf_button.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_FileIcon))
        self.about_vcf_button.setAutoRaise(True)
        self.about_vcf_button.setCursor(Qt.PointingHandCursor)
        self.about_vcf_button.setToolTip(
            t("About this VCF -- shows the currently loaded file's reference build, date, sample(s), and "
            "variant count, parsed from its own header.")
        )
        self.about_vcf_button.clicked.connect(self._on_about_vcf)
        row.addWidget(self.about_vcf_button)

        self.check_compat_button = QPushButton(t("Check compatibility"))
        self.check_compat_button.setObjectName("stepButton")
        self.check_compat_button.clicked.connect(self._on_check_compatibility)
        row.addWidget(self.check_compat_button)
        self.compat_help_button = self._make_help_button(t("Check compatibility"), COMPAT_HELP_TEXT)
        row.addWidget(self.compat_help_button)
        self.compat_status_label = QLabel("")
        row.addWidget(self.compat_status_label)

        self.clear_vcf_button = QPushButton(t("Clear"))
        self.clear_vcf_button.clicked.connect(self._on_clear_vcf)
        row.addWidget(self.clear_vcf_button)

        row.addStretch(1)
        layout.addLayout(row)

        # Limits the run itself to some genes, so a whole genome doesn't take
        # hours (the results table's gene-set filter only hides rows). A choice
        # made on purpose: Run stays off until one of the two is picked.
        score_genes_row = QHBoxLayout()
        score_genes_row.addWidget(QLabel(t("Genes to score:")))
        self._score_genes, self._score_genes_chosen = [], False
        self.all_genes_button = QPushButton(t("All genes"))
        self.all_genes_button.setObjectName("allGenesButton")
        self.all_genes_button.setCheckable(True)
        self.all_genes_button.clicked.connect(lambda: self._set_score_genes([], chosen=True))
        score_genes_row.addWidget(self.all_genes_button)
        self.choose_score_genes_button = QPushButton(t("Choose gene list..."))
        self.choose_score_genes_button.setObjectName("geneListButton")
        self.choose_score_genes_button.setCheckable(True)
        self.choose_score_genes_button.clicked.connect(self._on_choose_score_genes)
        score_genes_row.addWidget(self.choose_score_genes_button)
        score_genes_row.addWidget(self._make_help_button(t("Genes to score"), SCORE_GENES_HELP_TEXT))
        score_genes_row.addSpacing(12)
        score_genes_row.addWidget(self.run_button)
        self.score_genes_label = QLabel("")
        self.score_genes_label.setWordWrap(True)     # wraps rather than widening the window on a small screen
        score_genes_row.addWidget(self.score_genes_label, stretch=1)
        score_genes_row.addStretch(1)
        layout.addLayout(score_genes_row)
        self._set_score_genes([], chosen=False)

        return group

    def _build_options_group(self):
        """Settings. Only the two things that really are a per-run decision are
        on show -- raw/masked scoring and whether there is precomputed data.

        Everything else (the genome build, the reference FASTAs, SnpEff, MANE)
        either comes with Setup or is worked out from the VCF, so it lives
        behind "Advanced settings" and is only there for when that goes wrong.
        The build and its evidence are still reported in the open, read-only:
        scoring against the wrong build gives wrong answers quietly, so that
        one must never be out of sight.
        """
        group = QGroupBox(t("Settings"))
        layout = QVBoxLayout(group)

        self.build_combo = QComboBox()
        self.build_combo.addItems(list(BUILDS))
        self.build_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.build_combo.setMinimumWidth(100)
        self.build_combo.currentTextChanged.connect(self._update_run_enabled)
        self.build_combo.currentTextChanged.connect(self._on_build_changed)
        # Where the selected build came from -- the VCF's own header, when it
        # says (see _detect_vcf_build).
        self.build_source_label = QLabel("")
        self._detected_build = (None, None)
        self._vcf_present = False
        self._build_detect_timer = QTimer(self)
        self._build_detect_timer.setSingleShot(True)
        self._build_detect_timer.setInterval(300)
        self._build_detect_timer.timeout.connect(self._detect_vcf_build)

        top_row = QHBoxLayout()
        top_row.addWidget(QLabel(t("Mode:")))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems(list(MODES))
        self.mode_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.mode_combo.setMinimumWidth(100)
        top_row.addWidget(self.mode_combo)

        self.mode_help_button = self._make_help_button(t("Raw vs Masked scoring"), MODE_HELP_TEXT)
        top_row.addWidget(self.mode_help_button)

        top_row.addSpacing(24)
        # Read-only twin of the build row inside Advanced settings: which build
        # this run will use, and whether the VCF agrees.
        self.build_summary_label = QLabel("")
        top_row.addWidget(self.build_summary_label)

        top_row.addStretch(1)
        layout.addLayout(top_row)

        precomp_row = QHBoxLayout()
        self.use_precomputed_radio = QRadioButton(t("Use precomputed data"))
        self.no_precomputed_radio = QRadioButton(t("I don't have precomputed data"))
        self.precomputed_group = QButtonGroup(self)
        self.precomputed_group.addButton(self.use_precomputed_radio)
        self.precomputed_group.addButton(self.no_precomputed_radio)
        self.use_precomputed_radio.toggled.connect(self._update_precomputed_enabled)
        self.no_precomputed_radio.toggled.connect(self._update_precomputed_enabled)
        precomp_row.addWidget(self.use_precomputed_radio)
        precomp_row.addWidget(self.no_precomputed_radio)
        self.precomputed_help_button = self._make_help_button(t("Precomputed data"), PRECOMPUTED_HELP_TEXT)
        precomp_row.addWidget(self.precomputed_help_button)
        precomp_row.addStretch(1)
        layout.addLayout(precomp_row)

        precomp_dir_row = QHBoxLayout()
        precomp_dir_row.addWidget(QLabel(t("Precomputed data folder:")))
        self.precomputed_dir_edit = PathLineEdit()
        self.precomputed_dir_edit.pathDropped.connect(
            lambda p: self._handle_dropped_file(p, folder_edit=self.precomputed_dir_edit))
        self.precomputed_dir_edit.textChanged.connect(self._update_run_enabled)
        precomp_dir_row.addWidget(self.precomputed_dir_edit, stretch=1)
        self.precomputed_dir_browse = QPushButton(t("Browse..."))
        self.precomputed_dir_browse.clicked.connect(self._on_browse_precomputed_dir)
        precomp_dir_row.addWidget(self.precomputed_dir_browse)
        layout.addLayout(precomp_dir_row)

        self.skip_precomputed_checkbox = QCheckBox(t("Skip variants already in precomputed data"))
        layout.addWidget(self.skip_precomputed_checkbox)

        # --- advanced settings --------------------------------------------
        # The reference genomes, SnpEff and MANE all come with Setup and are
        # found automatically, so none of this is a question anyone should have
        # to answer. It is folded away, with one line saying whether everything
        # is in place; it opens by itself when something isn't.
        advanced_row = QHBoxLayout()
        self.advanced_button = QPushButton()
        self.advanced_button.setCheckable(True)
        self.advanced_button.toggled.connect(self._on_advanced_toggled)
        advanced_row.addWidget(self.advanced_button)
        self.setup_status_label = QLabel()
        self.setup_status_label.setWordWrap(True)
        advanced_row.addWidget(self.setup_status_label, stretch=1)
        layout.addLayout(advanced_row)

        self.advanced_widget = QWidget()
        advanced_layout = QVBoxLayout(self.advanced_widget)
        advanced_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.advanced_widget)

        build_row = QHBoxLayout()
        build_row.addWidget(QLabel(t("Genome build:")))
        build_row.addWidget(self.build_combo)
        build_row.addWidget(self.build_source_label)
        build_row.addStretch(1)
        advanced_layout.addLayout(build_row)

        # One reference FASTA per build; a run uses the selected build's
        # (_selected_fasta), and that row is marked (_update_fasta_rows).
        self.fasta_edits = {}
        self.fasta_row_labels = {}
        self.fasta_rows_widget = QWidget()
        fasta_rows_layout = QVBoxLayout(self.fasta_rows_widget)
        fasta_rows_layout.setContentsMargins(0, 0, 0, 0)
        advanced_layout.addWidget(self.fasta_rows_widget)
        label_width = self.fontMetrics().horizontalAdvance(t("{build} reference FASTA (used for this run):", build="hg38")) + 12
        for build in BUILDS:
            fasta_row = QHBoxLayout()
            label = QLabel()
            label.setMinimumWidth(label_width)
            fasta_row.addWidget(label)
            edit = PathLineEdit()
            edit.setPlaceholderText(t("{build}.fa -- keep its .fa.fai next to it", build=build))
            edit.setToolTip(t("The {build} reference FASTA -- Browse..., or drag the .fa file here.", build=build))
            edit.pathDropped.connect(lambda p, b=build: self._handle_dropped_file(p, fasta_build=b))
            edit.textChanged.connect(self._update_run_enabled)
            fasta_row.addWidget(edit, stretch=1)
            browse = QPushButton(t("Browse..."))
            browse.clicked.connect(lambda _checked=False, b=build: self._on_browse_fasta(b))
            fasta_row.addWidget(browse)
            download = QPushButton(t("Download..."))
            download.setToolTip(t("Download {build}.fa from UCSC, decompress it and build its index", build=build))
            download.clicked.connect(lambda _checked=False, b=build: self._on_download_reference(b))
            fasta_row.addWidget(download)
            # A "?" on every row also keeps both rows' fields exactly aligned.
            help_button = self._make_help_button(t("Reference FASTA"), FASTA_HELP_TEXT)
            fasta_row.addWidget(help_button)
            if build == BUILDS[0]:
                self.fasta_help_button = help_button
            self.fasta_edits[build] = edit
            self.fasta_row_labels[build] = label
            edit.textChanged.connect(self._update_setup_status)
            fasta_rows_layout.addLayout(fasta_row)
        self._update_fasta_rows()

        snpeff_opt_row = QHBoxLayout()
        self._snpeff_unavailable = False    # the box is off only because SnpEff isn't there, not by choice
        self.use_snpeff_checkbox = QCheckBox(t("Annotate with SnpEff"))
        self.use_snpeff_checkbox.setToolTip(
            t("On by default: SnpEff comes with the program. Turning it off drops the "
            "coding position and transcript (NM) columns, and makes a run a little faster.")
        )
        self.use_snpeff_checkbox.toggled.connect(self._on_use_snpeff_toggled)
        # Always on: SnpEff comes with Setup, so there is no setting. The box
        # stays as the run's internal switch (off only when SnpEff isn't there).
        self.use_snpeff_checkbox.setVisible(False)
        snpeff_opt_row.addWidget(self.use_snpeff_checkbox)
        snpeff_opt_row.addWidget(QLabel(t("SnpEff (always on): gives the coding position and transcript (NM) columns.")))
        self.use_mane_checkbox = QCheckBox(t("MANE Select"))
        self.use_mane_checkbox.setToolTip(
            t("On by default: where a gene has a MANE Select transcript, that NM_ is the one reported. Off: SnpEff's "
            "transcript order decides, curated (NM_) ones first.")
        )
        snpeff_opt_row.addWidget(self.use_mane_checkbox)
        self.snpeff_about_help_button = self._make_help_button(t("About SnpEff"), SNPEFF_ABOUT_HELP_TEXT)
        snpeff_opt_row.addWidget(self.snpeff_about_help_button)
        self.snpeff_help_button = self._make_help_button(t("SnpEff setup"), SNPEFF_HELP_TEXT)
        snpeff_opt_row.addWidget(self.snpeff_help_button)
        snpeff_opt_row.addStretch(1)
        advanced_layout.addLayout(snpeff_opt_row)


        # The SnpEff install-folder row and the MANE Select row share one
        # container widget, snpeff_config_widget, which is added to the
        # options layout below but starts out hidden
        # (setVisible(False)).
        self.snpeff_config_widget = QWidget()
        snpeff_config_layout = QVBoxLayout(self.snpeff_config_widget)
        snpeff_config_layout.setContentsMargins(0, 0, 0, 0)

        snpeff_row = QHBoxLayout()
        snpeff_row.addWidget(QLabel(t("SnpEff install folder:")))
        self.snpeff_dir_edit = PathLineEdit()
        self.snpeff_dir_edit.pathDropped.connect(lambda p: self._handle_dropped_file(p, folder_edit=self.snpeff_dir_edit))
        self.snpeff_dir_edit.setPlaceholderText(DEFAULT_SNPEFF_DIR)
        self.snpeff_dir_edit.textChanged.connect(self._update_run_enabled)
        snpeff_row.addWidget(self.snpeff_dir_edit, stretch=1)
        snpeff_dir_browse = QPushButton(t("Browse..."))
        snpeff_dir_browse.clicked.connect(self._on_browse_snpeff_dir)
        snpeff_row.addWidget(snpeff_dir_browse)
        self.download_snpeff_button = QPushButton(t("Download SnpEff..."))
        self.download_snpeff_button.clicked.connect(self._on_download_snpeff)
        snpeff_row.addWidget(self.download_snpeff_button)
        self.snpeff_status_label = QLabel("")
        snpeff_row.addWidget(self.snpeff_status_label)
        snpeff_config_layout.addLayout(snpeff_row)

        mane_row = QHBoxLayout()
        mane_label = QLabel(t("MANE Select folder:"))
        mane_label.setToolTip(t("Optional: improves the choice of transcript."))
        mane_row.addWidget(mane_label)
        # Where the MANE summary file is -- empty means the program's own
        # "mane" folder (DEFAULT_MANE_DIR). Passed to the run (worker.py).
        self.mane_dir_edit = PathLineEdit()
        self.mane_dir_edit.pathDropped.connect(lambda p: self._handle_dropped_file(p, folder_edit=self.mane_dir_edit))
        self.mane_dir_edit.setPlaceholderText(DEFAULT_MANE_DIR)
        self.mane_dir_edit.textChanged.connect(self._update_mane_status_label)
        mane_row.addWidget(self.mane_dir_edit, stretch=1)
        mane_browse = QPushButton(t("Browse..."))
        mane_browse.setToolTip(t("Select a MANE summary file you already have (MANE.GRCh38.vX.X.summary.txt.gz)"))
        mane_browse.clicked.connect(self._on_browse_mane_file)
        mane_row.addWidget(mane_browse)
        self.download_mane_button = QPushButton(t("Download MANE Select..."))
        self.download_mane_button.clicked.connect(self._on_download_mane)
        mane_row.addWidget(self.download_mane_button)
        self.mane_status_label = QLabel("")
        mane_row.addWidget(self.mane_status_label)
        self.mane_help_button = self._make_help_button(t("MANE Select data"), MANE_HELP_TEXT)
        mane_row.addWidget(self.mane_help_button)
        snpeff_config_layout.addLayout(mane_row)

        self.snpeff_config_widget.setVisible(False)
        advanced_layout.addWidget(self.snpeff_config_widget)

        # Turning SnpEff off hides its folder and the MANE row with it, because
        # MANE only refines the transcript SnpEff picked -- on its own it has
        # nothing to act on. Rather than leaving a gap where those rows were,
        # say what was lost and how to get it back.
        self.snpeff_off_label = QLabel(
            t("SnpEff wasn't found, so results have no transcript (NM) or coding position. It comes with the "
            "program: reinstalling restores it, or use \"Download SnpEff...\" below.")
        )
        self.snpeff_off_label.setWordWrap(True)
        self.snpeff_off_label.setStyleSheet(theme.STATUS_WARN)
        advanced_layout.addWidget(self.snpeff_off_label)

        self.advanced_widget.setVisible(False)
        return group

    def _build_run_row(self):
        row = QHBoxLayout()
        self.run_button = QPushButton(t("Run"))
        # The one button that starts the work; styled as the accent (theme.py).
        self.run_button.setObjectName("runButton")
        self.run_button.clicked.connect(self._on_run_clicked)
        # Placed in the Input VCF section, right of the Genes-to-score choice it waits for.

        # Only enabled while a run is going (see _set_running).
        self.pause_button = QPushButton(t("Pause"))
        self.pause_button.setToolTip(t("Pause the calculation after the current variant; click again to resume."))
        self.pause_button.setEnabled(False)
        self.pause_button.clicked.connect(self._on_pause_clicked)
        row.addWidget(self.pause_button)
        self.end_button = QPushButton(t("End"))
        self.end_button.setToolTip(t("End the calculation now (its results are discarded)."))
        self.end_button.setEnabled(False)
        self.end_button.clicked.connect(self._on_end_clicked)
        row.addWidget(self.end_button)

        self.summary_button = QPushButton(t("Summary"))
        self.summary_button.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_FileIcon))
        self.summary_button.setToolTip(t("Summary statistics for the most recently completed run (all results, not the filtered or paginated view)."))
        self.summary_button.setEnabled(False)
        self.summary_button.clicked.connect(self._on_show_summary)
        row.addWidget(self.summary_button)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        row.addWidget(self.progress_bar, stretch=1)

        self.progress_label = QLabel("")
        row.addWidget(self.progress_label, stretch=2)

        # CPU use while a run is going (see cpu_meter.py and _update_cpu_label).
        self.cpu_label = QLabel("")
        self.cpu_label.setToolTip(
            t("CPU use: this program's share of the whole processor, the whole computer's "
            "(including SnpEff's Java), and the time since the run started.")
        )
        self.cpu_label.setVisible(False)
        row.addWidget(self.cpu_label)
        self._cpu_meter = CpuMeter()
        self._cpu_timer = QTimer(self)
        self._cpu_timer.setInterval(1000)
        self._cpu_timer.timeout.connect(self._update_cpu_label)
        return row

    def _build_results_group(self):
        group = QGroupBox(t("Results"))
        layout = QVBoxLayout(group)

        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("max_score >="))
        self.threshold_buttons = {}
        for threshold in THRESHOLDS:
            btn = QPushButton(str(threshold))
            btn.setCheckable(True)
            btn.clicked.connect(lambda checked, val=threshold: self._on_threshold_clicked(val, checked))
            self.threshold_buttons[threshold] = btn
            filter_row.addWidget(btn)

        self.score_help_button = self._make_help_button(t("How to interpret scores"), SCORE_HELP_TEXT)
        filter_row.addWidget(self.score_help_button)
        filter_row.addSpacing(12)

        filter_row.addWidget(QLabel(t("Quick gene filter:")))
        self.quick_gene_edit = QLineEdit()
        self.quick_gene_edit.setPlaceholderText(t("gene symbol, e.g. BRCA1"))
        self.quick_gene_edit.returnPressed.connect(self._on_add_quick_gene)
        filter_row.addWidget(self.quick_gene_edit)
        self.add_quick_gene_button = QPushButton(t("Add"))
        self.add_quick_gene_button.clicked.connect(self._on_add_quick_gene)
        filter_row.addWidget(self.add_quick_gene_button)

        self.quick_gene_chips_layout = QHBoxLayout()
        self.quick_gene_chips_layout.setSpacing(4)
        filter_row.addLayout(self.quick_gene_chips_layout)
        filter_row.addStretch(1)
        layout.addLayout(filter_row)

        files_row = QHBoxLayout()           # its own row: spec 17 item 5, the tab fits 731 logical px
        self.export_button = QPushButton(t("Export visible rows..."))
        self.export_button.setToolTip(t("Exports only the rows currently shown, following the active gene and threshold filters and the sort order."))
        self.export_button.clicked.connect(self._on_export_clicked)
        files_row.addWidget(self.export_button)

        self.download_full_button = QPushButton(t("Download full results..."))
        self.download_full_button.setToolTip(t("Exports the complete result set from the most recent run, ignoring any active filters."))
        self.download_full_button.clicked.connect(self._on_download_full_results)
        files_row.addWidget(self.download_full_button)

        self.open_calculated_button = QPushButton(t("Open calculated data..."))
        self.open_calculated_button.setToolTip(t("Loads a previously saved results file (CLI or GUI output) directly into the table, skipping scoring."))
        self.open_calculated_button.clicked.connect(self._on_open_calculated_data)
        files_row.addWidget(self.open_calculated_button)

        self.clear_results_button = QPushButton(t("Clear results"))
        self.clear_results_button.clicked.connect(self._on_clear_results)
        files_row.addWidget(self.clear_results_button)
        files_row.addStretch(1)
        layout.addLayout(files_row)

        gene_set_row = QHBoxLayout()
        self.import_gene_list_button = QPushButton(t("Import gene list..."))
        self.import_gene_list_button.clicked.connect(self._on_import_gene_list)
        gene_set_row.addWidget(self.import_gene_list_button)

        self.gene_set_status_label = QLabel("")
        gene_set_row.addWidget(self.gene_set_status_label)

        self.clear_gene_set_button = QPushButton(t("Clear gene-set filter"))
        self.clear_gene_set_button.clicked.connect(self._on_clear_gene_set)
        self.clear_gene_set_button.setVisible(False)
        gene_set_row.addWidget(self.clear_gene_set_button)

        gene_set_row.addStretch(1)
        layout.addLayout(gene_set_row)

        count_row = QHBoxLayout()
        self.row_count_label = QLabel(t("0 rows"))
        count_row.addWidget(self.row_count_label)
        count_row.addStretch(1)
        # Spec 16 item 4: the selected variant(s) to the Sequence-to-cDNA tab.
        self.transcript_effect_button = QPushButton(t("Analyse transcript effect..."))
        self.transcript_effect_button.setToolTip(
            t("Sends the selected variant to the Sequence to cDNA tab: the sample is built from the reference with "
              "this variant only, and its predicted transcripts are compared with the reference and the WT "
              "prediction. Up to {n} selected variants are run one after another.", n=vcf_bridge.MAX_QUEUE))
        self.transcript_effect_button.clicked.connect(lambda: self.analyse_transcript_effect())
        count_row.addWidget(self.transcript_effect_button)
        layout.addLayout(count_row)

        self.no_results_label = QLabel("")
        self.no_results_label.setStyleSheet(f"color: {theme.CRIT}; font-style: italic;")
        self.no_results_label.setVisible(False)
        layout.addWidget(self.no_results_label)

        self.table_model = ScoreTableModel()
        self.proxy_model = ScoreFilterProxyModel(self.table_model, on_changed=self._on_proxy_changed)

        # page_table_model is a second ScoreTableModel; the table view below is
        # attached to it rather than to proxy_model (see _refresh_page and
        # RESULTS_PAGE_SIZE).

        self.page_table_model = ScoreTableModel()

        self.table_view = QTableView()
        self.table_view.setModel(self.page_table_model)
        self.table_view.setSortingEnabled(True)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table_view.customContextMenuRequested.connect(self._on_results_context_menu)
        self.table_view.horizontalHeader().setSectionsMovable(True)

        # Max Score is drawn as a number with a bar behind it, and the first
        # column carries a severity stripe, so a long result set can be
        # scanned. Presentation only -- see score_delegate.py.
        self._score_delegate = ScoreBarDelegate(self.table_view)
        self.table_view.setItemDelegateForColumn(MAX_SCORE_COLUMN, self._score_delegate)
        self._stripe_delegate = SeverityStripeDelegate(MAX_SCORE_COLUMN, self.table_view)
        self.table_view.setItemDelegateForColumn(0, self._stripe_delegate)

        self.table_view.horizontalHeader().sortIndicatorChanged.connect(self._on_sort_indicator_changed)


        self.table_view.horizontalHeader().setHighlightSections(False)
        self.table_view.verticalHeader().setHighlightSections(False)

        # The "ref" and "alt" columns are hidden in the view: each one's index
        # is looked up in COLUMNS (3-tuples whose middle item is the row
        # attribute name).

        for hidden_attr in ("ref", "alt"):
            hidden_col = next(i for i, (_, attr, _) in enumerate(COLUMNS) if attr == hidden_attr)
            self.table_view.setColumnHidden(hidden_col, True)
        layout.addWidget(self.table_view, stretch=1)

        pagination_row = QHBoxLayout()
        self.prev_page_button = QPushButton(t("< Previous"))
        self.prev_page_button.clicked.connect(self._on_prev_page)
        pagination_row.addWidget(self.prev_page_button)
        self.page_indicator_label = QLabel(t("Page 1 of 1"))
        pagination_row.addWidget(self.page_indicator_label)
        self.next_page_button = QPushButton(t("Next >"))
        self.next_page_button.clicked.connect(self._on_next_page)
        pagination_row.addWidget(self.next_page_button)
        pagination_row.addStretch(1)
        layout.addLayout(pagination_row)

        return group

    def _load_settings(self):
        settings = config.load()
        # Fills in any reference FASTA that isn't set (or whose saved file has
        # gone) from the usual places -- see reference_locator. self.
        # detected_references names the builds that came from there, for the
        # summary line and so _save_settings can write them back.
        paths, self.detected_references = reference_locator.resolve(settings)
        for build, edit in self.fasta_edits.items():
            edit.setText(paths[build])
        self.precomputed_dir_edit.setText(settings["precomputed_dir"])
        self.snpeff_dir_edit.setText(settings["snpeff_dir"])
        self.mane_dir_edit.setText(settings["mane_dir"])
        # On unless SnpEff genuinely isn't there: it ships with Setup, but a
        # run from source (or a half-finished install) has no database, and an
        # unticked box the user can see beats a disabled Run button.
        self.use_snpeff_checkbox.setChecked(True)
        self.use_mane_checkbox.setChecked(settings["use_mane"])
        self.use_mane_checkbox.setEnabled(True)
        if not self._snpeff_ready():
            self.use_snpeff_checkbox.setChecked(False)
            # Saved as it was, so annotating comes back by itself once SnpEff is there.
            self._snpeff_unavailable = True
        self.skip_precomputed_checkbox.setChecked(settings["skip_precomputed"])
        if settings["build"] in BUILDS:
            self.build_combo.setCurrentText(settings["build"])
        if settings["mode"] in MODES:
            self.mode_combo.setCurrentText(settings["mode"])
        if settings["precomputed_mode"] == "use":
            self.use_precomputed_radio.setChecked(True)
        elif settings["precomputed_mode"] == "none":
            self.no_precomputed_radio.setChecked(True)

        if settings["column_order"]:
            self._apply_column_order(settings["column_order"])
        if isinstance(settings["score_genes"], list):
            self._set_score_genes([g for g in settings["score_genes"] if isinstance(g, str)], chosen=False)

        # Only a list of the right shape, so an edited or older settings file
        # can't leave a pane at zero height with no way to get it back.
        sizes = settings["pane_sizes"]
        if (isinstance(sizes, list) and len(sizes) == self.panes.count()
                and all(isinstance(n, int) and n > 0 for n in sizes)):
            self.panes.setSizes(sizes)
            # Saved while Advanced settings were open (or at another window size): they start closed, and the
            # pane is fitted to its contents once the window is shown (showEvent).

    def _save_settings(self):
        if self.use_precomputed_radio.isChecked():
            precomputed_mode = "use"
        elif self.no_precomputed_radio.isChecked():
            precomputed_mode = "none"
        else:
            precomputed_mode = None
        # Starts from config.load() and update()s only the keys below, so any
        # other keys already saved are kept.
        settings = config.load()
        settings.update({
            **{config.fasta_key(build): edit.text().strip() for build, edit in self.fasta_edits.items()},
            "fasta_path": "",  # the old single field, migrated by config.load()
            "precomputed_mode": precomputed_mode,
            "precomputed_dir": self.precomputed_dir_edit.text().strip(),
            "build": self.build_combo.currentText(),
            "mode": self.mode_combo.currentText(),
            "skip_precomputed": self.skip_precomputed_checkbox.isChecked(),
            "column_order": self._current_column_order(),
            "column_order_v2_applied": True,
            "pane_sizes": self.panes.sizes(),
            "score_genes": self._score_genes,
            "snpeff_dir": self.snpeff_dir_edit.text().strip(),
            "mane_dir": self.mane_dir_edit.text().strip(),
            "use_snpeff": True,
            "use_mane": self.use_mane_checkbox.isChecked(),
            "snpeff_default_applied": True,
        })
        config.save(settings)

    def _current_column_order(self):
        """ScoreRow attr names in the exact left-to-right order currently shown in the
        table -- reflects any drag-to-reorder the user has done via the movable
        section headers. Used both to persist column order across launches and to
        make CSV export match what's actually on screen."""
        header = self.table_view.horizontalHeader()
        return [COLUMNS[header.logicalIndex(visual)][1] for visual in range(header.count())]

    def _apply_column_order(self, order_attrs):
        """Reorders the header's visual sections to match a saved attr-name order.
        Unknown/missing names (e.g. a column added in a later version, or one the
        user never rearranged) are handled gracefully: unrecognized saved names are
        ignored, and any current column not mentioned in the saved order is placed
        after the column it follows in COLUMNS (at the end if that one is missing too),
        rather than dropped."""
        header = self.table_view.horizontalHeader()
        attr_to_logical = {attr: i for i, (_, attr, _) in enumerate(COLUMNS)}
        target = [attr_to_logical[a] for a in order_attrs if a in attr_to_logical]
        for logical in range(len(COLUMNS)):
            if logical not in target:
                # A column added since: next to the one it follows by default.
                after = target.index(logical - 1) + 1 if logical - 1 in target else len(target)
                target.insert(after, logical)
        for visual_pos, logical_index in enumerate(target):
            current_visual = header.visualIndex(logical_index)
            if current_visual != visual_pos:
                header.moveSection(current_visual, visual_pos)

    def closeEvent(self, event):
        if self._worker is not None and self._worker.isRunning():
            reply = QMessageBox.warning(
                self, t("Calculation still running"),
                t("The calculation is not finished yet. Are you sure you want to exit?\n\n"
                "The run will be stopped and its results discarded."),
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                event.ignore()
                return
            # Stops at the next variant (SnpEff is killed at once); a thread
            # still running when the window goes would crash the program.
            self.progress_label.setText(t("Stopping the calculation..."))
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self._worker.stop()
                if not self._worker.wait(60000):
                    self._worker.terminate()
                    self._worker.wait(5000)
            finally:
                QApplication.restoreOverrideCursor()
        if self.refdb.proc is not None:
            reply = QMessageBox.question(
                self, t("Reference database build running"),
                t("A reference database build is running. Close the program anyway?\n\n"
                "The build stops now and continues from the last finished gene the next time you start the program."),
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                event.ignore()
                return
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.refdb.shutdown()
            finally:
                QApplication.restoreOverrideCursor()
        # A Sequence-to-cDNA run stops at its next model call; its thread must
        # not outlive the window.
        if self.cdna_tab._run_worker is not None:
            self.cdna_tab._run_worker.stop()
            self.cdna_tab._run_worker.wait(30000)
        self.cdna_tab._save_settings()
        self._save_settings()
        self._cleanup_temp_vcf()
        super().closeEvent(event)

    def _on_load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, t("Load VCF file"), "", t("VCF files (*.vcf *.vcf.gz *.vcf.bgz *.bgz);;All files (*)")
        )
        if not path:
            return
        self._load_vcf_from_path(path)

    def _on_file_dropped(self, path):
        self._load_vcf_from_path(path)

    # Drag and drop anywhere on the window (the input box and the path fields
    # handle drops onto themselves, and pass other files on to here).
    def dragEnterEvent(self, event):
        files = dropped_local_files(event.mimeData())
        if files:
            event.acceptProposedAction()
            self.vcf_text.set_drop_highlight(bool(dropped_vcf_path(event.mimeData())))
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if dropped_local_files(event.mimeData()):
            event.acceptProposedAction()
        else:
            event.ignore()

    def dragLeaveEvent(self, event):
        self.vcf_text.set_drop_highlight(False)

    def dropEvent(self, event):
        self.vcf_text.set_drop_highlight(False)
        files = dropped_local_files(event.mimeData())
        if not files:
            event.ignore()
            return
        event.acceptProposedAction()
        self._handle_dropped_file(dropped_vcf_path(event.mimeData()) or files[0])

    def _handle_dropped_file(self, path, fasta_build=None, folder_edit=None):
        """What a dragged-in file is for, by its name: a VCF loads (like "Load
        file..."), a FASTA goes into a reference row (fasta_build's when dropped
        onto that row, else the build its name says, else the selected one), a
        MANE summary sets the MANE folder, a reference score database (its
        folder or manifest.json) is offered for import, and anything dropped
        onto a folder field sets that field."""
        name = path.lower()
        database = database_folder_of(path) if folder_edit is None else None
        if database:
            self._open_refdb_dialog(import_folder=database)
        elif self.tabs.currentWidget() is self.cdna_tab and folder_edit is None and fasta_build is None:
            self.cdna_tab.handle_drop(path)
        elif name.endswith(VCF_SUFFIXES):
            self._load_vcf_from_path(path)
        elif name.endswith(FASTA_SUFFIXES):
            build = fasta_build or config.guess_build_from_name(path) or self.build_combo.currentText()
            self._set_fasta(build, path)
        elif is_summary_filename(path):
            self.mane_dir_edit.setText(os.path.dirname(path))
        elif folder_edit is not None:
            folder_edit.setText(path if os.path.isdir(path) else os.path.dirname(path))
        else:
            QMessageBox.information(
                self, t("Can't use this file"),
                t("{name} can't be loaded here.\n\n"
                  "Drag in a VCF (.vcf, .vcf.gz) to score it, a reference FASTA (.fa, .fasta, .fna) "
                  "for its build's row, or a MANE summary file (MANE.GRCh38.vX.X.summary.txt.gz).",
                  name=os.path.basename(path)),
            )

    def _load_vcf_from_path(self, path):
        """Shared by "Load file..." and drag-and-drop onto the input area --
        both go through identical read/decompress/decode + validation."""
        try:
            self.setCursor(Qt.WaitCursor)
            text = _read_vcf_text(path)
        except RuntimeError as exc:
            QMessageBox.critical(self, t("Failed to load file"), f"{path}\n\n{exc}")
            return
        finally:
            self.unsetCursor()
        # A newly loaded file always gets its build detected afresh, even if it
        # names the same build as the previous one (see _detect_vcf_build).
        self._detected_build = (None, None)
        self.vcf_text.setPlainText(text)
        self.load_status_label.setText(t("Loaded: {name}", name=os.path.basename(path)))
        self._loaded_vcf_basename = strip_vcf_extension(os.path.basename(path))

    def _on_clear_vcf(self):
        """Resets the Input VCF section back to its initial empty state. Independent
        of the Results section on purpose -- clearing what's about to be run shouldn't
        discard results still being viewed from a prior run."""
        self.vcf_text.clear()
        self.load_status_label.setText("")
        self._set_score_genes(self._score_genes, chosen=False)
        self.compat_status_label.setStyleSheet("")
        # The compat label's text and tooltip are cleared by _on_vcf_text_changed,
        # which vcf_text.clear() above triggers through the textChanged signal.

    def _on_vcf_text_changed(self):
        self.compat_status_label.setText("")
        self.compat_status_label.setToolTip("")
        # Re-detect the build shortly after the text stops changing (typing or
        # pasting fires this for every change).
        self._build_detect_timer.start()

    def _detect_vcf_build(self):
        """Selects the build the VCF's own header names (vcf_info.detect_vcf_build)
        -- but only when that result changes, i.e. for a different VCF, so a
        build the user picked by hand isn't overridden while they edit the same
        one. A header that doesn't say leaves the selection alone."""
        text = self.vcf_text.toPlainText()
        self._vcf_present = bool(text.strip())
        detected = detect_vcf_build(text) if self._vcf_present else (None, None)
        if detected != self._detected_build:
            self._detected_build = detected
            if detected[0] in BUILDS and detected[0] != self.build_combo.currentText():
                self.build_combo.setCurrentText(detected[0])
        self._update_build_source_label()

    def _update_build_source_label(self):
        build, evidence = self._detected_build
        if not self._vcf_present:
            text, style = "", ""
        elif build is None:
            text, style = "⚠ " + t("build not stated in the VCF -- make sure this is right"), theme.STATUS_WARN
        elif build == self.build_combo.currentText():
            text, style = "✓ " + t("from the VCF ({evidence})", evidence=t(evidence)), theme.STATUS_OK
        else:
            text, style = "⚠ " + t("the VCF says {build} ({evidence})", build=build, evidence=t(evidence)), \
                theme.STATUS_CRIT
        self.build_source_label.setText(text)
        self.build_source_label.setStyleSheet(style)
        # The same thing again in the open, where the build combo no longer is.
        summary = getattr(self, "build_summary_label", None)
        if summary is not None:
            selected = self.build_combo.currentText()
            summary.setText(t("Genome build: {build}", build=f"{selected} {text}".rstrip() if text else selected))
            summary.setStyleSheet(style or theme.STATUS_MUTED)
            summary.setToolTip(t("Chosen from the VCF's own header. To set it by hand, "
                               "open Advanced settings."))

    def _on_build_changed(self, _build=None):
        self._update_fasta_rows()
        if hasattr(self, "build_source_label"):
            self._update_build_source_label()

    def _update_fasta_rows(self):
        """Marks the FASTA row the next run will use (the selected build's)."""
        labels = getattr(self, "fasta_row_labels", None)
        if not labels:
            return
        selected = self.build_combo.currentText()
        for build, label in labels.items():
            active = build == selected
            label.setText(t("{build} reference FASTA (used for this run):", build=build) if active
                          else t("{build} reference FASTA:", build=build))
            label.setStyleSheet("font-weight: 600;" if active else theme.STATUS_MUTED)
        self._update_setup_status()

    def _on_advanced_toggled(self, checked):
        self.advanced_widget.setVisible(checked)
        self.advanced_button.setText(t("Hide advanced settings") if checked else t("Advanced settings..."))
        if not checked:
            self._fit_settings_pane()

    def showEvent(self, event):
        super().showEvent(event)
        # The Settings pane at its contents' height, whatever was saved or the default gave it -- once the window
        # is laid out, as the splitter's sizes mean nothing before.
        if not getattr(self, "_fitted_once", False):
            self._fitted_once = True
            QTimer.singleShot(0, self._fit_settings_pane)

    def _fit_settings_pane(self):
        """Opening Advanced settings makes the Settings pane taller; closing them
        must give that height back (to the results), not leave it empty."""
        pane = self.panes.widget(2)
        pane.layout().activate()
        need = pane.sizeHint().height()
        sizes = self.panes.sizes()
        if sizes[2] > need:
            self.panes.setSizes([sizes[0] + sizes[2] - need, sizes[1], need])

    def setup_problems(self):
        """What, if anything, stops this being a normal run -- in the order a
        user would have to deal with it. Everything here comes with Setup, so
        an empty list is the expected case and the only one most people see."""
        problems = []
        if not self._selected_fasta():
            build = self.build_combo.currentText()
            problems.append(t("the {build} reference genome hasn't been found", build=build))
        # getattr: the FASTA rows are built (and ask for this) before the
        # SnpEff row below them exists.
        if getattr(self, "use_snpeff_checkbox", None) and self.use_snpeff_checkbox.isChecked():
            ready, short_message, _detail = self._snpeff_status()
            if not ready:
                problems.append(short_message[0].lower() + short_message[1:] if language() == "en" else short_message)
        return problems

    def _update_setup_status(self, _text=None):
        """The one line next to "Advanced settings...": either everything is in
        place (and where the genomes came from), or what isn't."""
        label = getattr(self, "setup_status_label", None)
        if label is None or not getattr(self, "fasta_edits", None):
            return
        problems = self.setup_problems()
        if problems:
            text = "⚠ " + "; ".join(problems) + t(" -- open Advanced settings.")
            style = theme.STATUS_CRIT
        else:
            # Only the paths still holding what was detected -- once the user
            # browses to a file of their own, the window stops claiming it
            # found it.
            detected = getattr(self, "detected_references", {})
            set_builds = [b for b in BUILDS if self.fasta_edits[b].text().strip()]
            folders = {os.path.dirname(detected[b]) for b in set_builds
                       if detected.get(b) == self.fasta_edits[b].text().strip()}
            where = t(", found in {folder}", folder=sorted(folders)[0]) if len(folders) == 1 else ""
            parts = [t("reference genome {builds}", builds=t(" and ").join(set_builds))]
            if self.use_snpeff_checkbox.isChecked():
                parts.append("SnpEff")
            text = t("Ready: {what}{where}.", what=", ".join(parts), where=where)
            style = theme.STATUS_MUTED
        label.setText(text)
        label.setStyleSheet(style)

    def _init_advanced_section(self):
        """Called once, after the settings are loaded: advanced settings stay
        folded away when everything is in place, and start open when something
        has to be sorted out."""
        self.advanced_button.setChecked(bool(self.setup_problems()))
        # setChecked() emits nothing when the value doesn't change, so set the
        # button's text and the section's visibility here either way.
        self._on_advanced_toggled(self.advanced_button.isChecked())

    def _reveal_reference_section(self):
        """Opens advanced settings, so a file that has just been put into one
        of its rows (dropped, browsed to, downloaded) is somewhere visible."""
        self.advanced_button.setChecked(True)
        self._on_advanced_toggled(True)

    def reference_paths(self):
        """{build: FASTA path} as the window currently has them, including the
        ones found automatically. Used by the first-launch reference note, so
        it isn't shown about files the window already has."""
        return {build: edit.text().strip() for build, edit in self.fasta_edits.items()}

    def _selected_fasta(self):
        return self.fasta_edits[self.build_combo.currentText()].text().strip()

    def _on_about_vcf(self):
        text = self.vcf_text.toPlainText()
        if not text.strip():
            QMessageBox.information(self, t("About this VCF"), t("No VCF is currently loaded."))
            return

        info = summarize_vcf_text(text)
        lines = []

        if info.detected_build_raw:
            lines.append(t("Reference (from header): {ref}", ref=info.detected_build_raw))
        build, evidence = detect_vcf_build(text)
        selected_build = self.build_combo.currentText()
        if build:
            lines.append(t("Genome build: {build}", build=f"{build} ({t(evidence)})"))
            if build != selected_build:
                lines.append("⚠ " + t("The file is {build}, but {selected} is currently selected.", build=build,
                                      selected=selected_build))
        else:
            lines.append(t("Genome build: not stated in the file (no ##reference line or ##contig lengths)"))

        if info.file_date:
            lines.append(t("File date: {date}", date=info.file_date))
        if info.samples:
            label = t("Sample") if len(info.samples) == 1 else t("Samples")
            lines.append(f"{label}: {', '.join(info.samples)}")
        lines.append(t("Variant records: {n:,}", n=info.variant_count))
        if info.source:
            lines.append(t("Source: {source}", source=info.source))
        if info.contig_count:
            lines.append(t("Contigs: {n}", n=info.contig_count))

        title = t("About this VCF")
        if self._loaded_vcf_basename:
            title += f" -- {self._loaded_vcf_basename}"
        QMessageBox.information(self, title, "\n".join(lines))

    def _on_check_compatibility(self):
        result = check_vcf_text(self.vcf_text.toPlainText())
        prefix = "✓ " + t("PASS") if result.ok else "✗ " + t("FAIL")
        color = "#1a7f37" if result.ok else "#b00020"
        self.compat_status_label.setText(f"{prefix}: {result.summary}")
        self.compat_status_label.setStyleSheet(f"color: {color}; font-weight: bold;")
        if result.details:
            self.compat_status_label.setToolTip("\n".join(result.details))
        else:
            self.compat_status_label.setToolTip("")

    def _on_import_gene_list(self):
        dialog = GeneImportDialog(self)
        if dialog.exec() == QDialog.Accepted:
            genes = dialog.validated_genes()
            self._gene_set_filter = set(genes)
            self.proxy_model.set_gene_set(self._gene_set_filter)
            self.gene_set_status_label.setText(t("Gene-set filter: {n} gene(s) active", n=len(genes)))
            self.clear_gene_set_button.setVisible(True)

    def _on_choose_score_genes(self):
        # The last list (saved in the settings) is offered again.
        dialog = GeneImportDialog(self, title=t("Genes to score"), apply_text=t("Score only these genes"),
                                  genes=self._score_genes)
        genes = dialog.validated_genes() if dialog.exec() == QDialog.Accepted else []
        if genes:
            self._set_score_genes(genes, chosen=True)
        else:   # cancelled: the click must not leave the button looking chosen
            self._set_score_genes(self._score_genes, chosen=self._score_genes_chosen)

    def _set_score_genes(self, genes, chosen=True):
        """genes [] = all genes. chosen=False: nothing picked yet for this VCF,
        so Run waits (the list is kept, to be offered again)."""
        self._score_genes, self._score_genes_chosen = list(genes), chosen
        all_genes, a_list = chosen and not genes, chosen and bool(genes)
        self.all_genes_button.setChecked(all_genes)
        self.choose_score_genes_button.setChecked(a_list)
        self.all_genes_button.setText(("✓ " if all_genes else "") + t("All genes"))
        self.choose_score_genes_button.setText(("✓ " if a_list else "") + t("Choose gene list..."))
        if a_list:
            shown = (t("{genes} and {n} more", genes=", ".join(genes[:6]), n=len(genes) - 6) if len(genes) > 6
                     else ", ".join(genes))
            self.score_genes_label.setText(t("{n} gene(s): {genes}", n=len(genes), genes=shown))
            self.score_genes_label.setStyleSheet(theme.STATUS_OK)
            self.score_genes_label.setToolTip(", ".join(genes))
        else:
            self.score_genes_label.setText(t("all genes") if all_genes else t("Choose one to enable Run."))
            self.score_genes_label.setStyleSheet(theme.STATUS_MUTED if all_genes else theme.STATUS_WARN)
            self.score_genes_label.setToolTip("")
        if hasattr(self, "fasta_edits"):    # not while the window is still being built
            self._update_run_enabled()

    def _on_clear_gene_set(self):
        self._gene_set_filter = None
        self.proxy_model.set_gene_set(None)
        self.gene_set_status_label.setText("")
        self.clear_gene_set_button.setVisible(False)

    def _on_browse_precomputed_dir(self):
        path = QFileDialog.getExistingDirectory(self, t("Select precomputed-scores folder"))
        if path:
            self.precomputed_dir_edit.setText(path)

    def _on_browse_snpeff_dir(self):
        path = QFileDialog.getExistingDirectory(self, t("Select SnpEff install folder"), self.snpeff_dir_edit.text() or DEFAULT_SNPEFF_DIR)
        if path:
            self.snpeff_dir_edit.setText(path)

    def _on_browse_fasta(self, build):
        path, _ = QFileDialog.getOpenFileName(
            self, t("Select the {build} reference FASTA", build=build), os.path.dirname(self.fasta_edits[build].text()),
            t("FASTA files (*.fa *.fasta *.fna);;All files (*)"),
        )
        if path:
            self._set_fasta(build, path)

    def _set_fasta(self, build, path):
        """Puts a FASTA into a build's row -- from Browse... or a drop -- asking
        first when its file name names the other build."""
        named = config.guess_build_from_name(path)
        if named and named != build:
            reply = QMessageBox.question(
                self, t("Different build?"),
                t("{name} looks like an {named} file, but this is the {build} field.\n\nUse it for {build} anyway?",
                  name=os.path.basename(path), named=named, build=build),
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return
        self.fasta_edits[build].setText(path)
        self._reveal_reference_section()

    def _on_download_reference(self, build):
        dialog = ReferenceDownloadDialog(self, default_build=build)
        if dialog.exec() == QDialog.Accepted:
            path = dialog.result_path()
            if path:
                self.fasta_edits[dialog.result_build()].setText(path)
                self._reveal_reference_section()

    def _update_precomputed_enabled(self):
        use_precomputed = self.use_precomputed_radio.isChecked()
        self.precomputed_dir_edit.setEnabled(use_precomputed)
        self.precomputed_dir_browse.setEnabled(use_precomputed)
        self.skip_precomputed_checkbox.setEnabled(use_precomputed)
        self._update_run_enabled()

    def _update_run_enabled(self):
        has_vcf = bool(self.vcf_text.toPlainText().strip())
        has_fasta = bool(self._selected_fasta())
        precomputed_chosen = self.use_precomputed_radio.isChecked() or self.no_precomputed_radio.isChecked()
        precomputed_dir_ok = True
        if self.use_precomputed_radio.isChecked():
            precomputed_dir_ok = bool(self.precomputed_dir_edit.text().strip())

        self._update_snpeff_status_label()
        snpeff_ready, snpeff_short_message, _ = self._snpeff_status()
        if self.use_snpeff_checkbox.isChecked():
            self._update_mane_status_label()

        # Everything still missing is listed in the Run button's tooltip
        # ("Before you can run: ..."); Run is enabled only when nothing is
        # missing and no worker is running.
        missing = []
        if not has_vcf:
            missing.append(t("load a VCF"))
        if not has_fasta:
            missing.append(t("set the {build} reference FASTA", build=self.build_combo.currentText()))
        if not precomputed_chosen:
            missing.append(t("choose a precomputed-data option"))
        elif not precomputed_dir_ok:
            missing.append(t("set the Precomputed data folder"))
        if not snpeff_ready:
            missing.append(snpeff_short_message)
        if not self._score_genes_chosen:
            missing.append(t("choose All genes or a gene list"))
        self.run_button.setToolTip(
            t("Before you can run: {what}.", what="; ".join(missing)) if missing else ""
        )

        enabled = (has_vcf and has_fasta and precomputed_chosen and precomputed_dir_ok and snpeff_ready
                   and self._score_genes_chosen)
        self.run_button.setEnabled(enabled and self._worker is None)
        self._update_setup_status()

    def _snpeff_ready(self):
        return self._snpeff_status()[0]

    def _snpeff_status(self):
        """Returns (ready, short_message, detail). ready/short_message/detail are
        all "" and True when SnpEff isn't required (checkbox unchecked -- the
        FASTA is the only always-required reference file) or is fully set up
        for the CURRENTLY SELECTED build. This is re-evaluated fresh every call
        against the real filesystem (see check_snpeff_setup) -- there is no
        cached/shared "is SnpEff installed" flag, so hg19 and hg38 (each with
        their own database directory, e.g. data/GRCh37.p13 vs data/GRCh38.p14)
        are always checked independently, keyed to build_combo's current value.

        short_message names the specific build when that's what's missing
        (e.g. "SnpEff database for hg38 not installed"), rather than an
        ambiguous "not installed" that reads as if SnpEff itself is gone --
        distinct from the case where SnpEff/its jar isn't installed at all.
        detail is the fuller path/command-level message, for a tooltip.
        """
        if not self.use_snpeff_checkbox.isChecked():
            return True, "", ""
        snpeff_dir = self.snpeff_dir_edit.text().strip() or DEFAULT_SNPEFF_DIR
        build = self.build_combo.currentText()
        try:
            check_snpeff_setup(snpeff_dir, build)
            find_java(snpeff_dir)
            return True, "", ""
        except JavaNotFoundError as exc:
            # SnpEff is there but no usable Java: "Download SnpEff..." now fetches
            # just the Java part (installed parts are skipped).
            return False, t("Java for SnpEff isn't installed"), str(exc)
        except SnpEffNotSetUpError as exc:
            detail = str(exc)
            if not os.path.isfile(snpeff_jar_path(snpeff_dir)):
                short_message = t("SnpEff isn't installed")
            else:
                short_message = t("SnpEff database for {build} not installed", build=build)
            return False, short_message, detail

    def _update_snpeff_status_label(self):
        ready, short_message, detail = self._snpeff_status()
        if not self.use_snpeff_checkbox.isChecked():
            self.snpeff_status_label.setText("")
            self.snpeff_status_label.setToolTip("")
            return
        if ready:
            self.snpeff_status_label.setText(t("✓ Installed"))
            self.snpeff_status_label.setStyleSheet(theme.STATUS_OK)
            self.snpeff_status_label.setToolTip("")

        #
        # Not ready: the short message plus a pointer to "Download SnpEff...",
        # shown in gray (#5f5f5f), with the fuller detail message (which names
        # the missing path/command) as the tooltip.
        else:
            self.snpeff_status_label.setText(short_message + " -- " + t('click "Download SnpEff..."'))
            self.snpeff_status_label.setStyleSheet(theme.STATUS_MUTED)
            self.snpeff_status_label.setToolTip(detail)

    def _mane_dir(self):
        """Where to look for the MANE summary: the field if it is set,
        otherwise the "mane" folder next to the app, and failing that the one
        beside the chosen SnpEff install -- Setup puts snpeff, java and mane
        side by side, so that sibling is where it is when the app itself is
        running from somewhere else."""
        chosen = self.mane_dir_edit.text().strip()
        if chosen:
            return chosen
        if find_cached_summary(DEFAULT_MANE_DIR):
            return DEFAULT_MANE_DIR
        beside = sibling_mane_dir(self.snpeff_dir_edit.text().strip() or DEFAULT_SNPEFF_DIR)
        return beside or DEFAULT_MANE_DIR

    def _update_mane_status_label(self):
        """Whether the MANE folder (field, or the default) has a summary file
        the run will find -- the same lookup the pipeline does."""
        mane_dir = self._mane_dir()
        path = find_cached_summary(mane_dir)
        if path:
            # Short, so it fits a small window; the full path is the tooltip.
            version = os.path.basename(path).split(".summary")[0].rsplit(".v", 1)[-1]
            self.mane_status_label.setText("✓ " + t("Found (MANE v{version})", version=version))
            self.mane_status_label.setStyleSheet(theme.STATUS_OK)
            self.mane_status_label.setToolTip(path)
        else:
            self.mane_status_label.setText(t("Not found"))
            self.mane_status_label.setStyleSheet(theme.STATUS_MUTED)
            self.mane_status_label.setToolTip(t("No MANE.GRCh38.vX.X.summary.txt.gz in {folder}", folder=mane_dir))

    def _on_browse_mane_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, t("Select the MANE Select summary file"), self._mane_dir(),
            t("MANE summary (MANE.GRCh38.v*.summary.txt.gz);;All files (*)"),
        )
        if not path:
            return
        if not is_summary_filename(path):
            QMessageBox.warning(
                self, t("Not a MANE summary file"),
                t("{name} isn't a MANE summary file.\n\nExpected a name like "
                  "MANE.GRCh38.v1.5.summary.txt.gz (from ftp.ncbi.nlm.nih.gov/refseq/MANE).",
                  name=os.path.basename(path)),
            )
            return
        self.mane_dir_edit.setText(os.path.dirname(path))

    def _on_use_snpeff_toggled(self, checked):
        self._snpeff_unavailable = False
        self.use_mane_checkbox.setEnabled(checked)
        self.snpeff_config_widget.setVisible(checked)
        self.snpeff_off_label.setVisible(not checked)
        if checked:
            self._update_mane_status_label()
        self._update_run_enabled()

    def _on_download_snpeff(self):
        dest_dir = self.snpeff_dir_edit.text().strip() or DEFAULT_SNPEFF_DIR
        dialog = SnpEffDownloadDialog(self, build=self.build_combo.currentText(), dest_dir=dest_dir)
        if dialog.exec() == QDialog.Accepted:
            self.snpeff_dir_edit.setText(dialog.dir_edit.text().strip())
        self._update_run_enabled()

    def _on_download_mane(self):
        dialog = ManeDownloadDialog(self, dest_dir=self._mane_dir())
        dialog.exec()
        self._update_mane_status_label()

    def _on_threshold_clicked(self, value, checked):
        for t, btn in self.threshold_buttons.items():
            if t != value:
                btn.setChecked(False)
        self.proxy_model.set_threshold(value if checked else None)

    def _on_add_quick_gene(self):
        text = self.quick_gene_edit.text().strip()
        if not text:
            return
        self.quick_gene_edit.clear()
        if text.lower() in (g.lower() for g in self._quick_genes):
            return
        self._quick_genes.append(text)
        self._add_gene_chip(text)
        self._apply_quick_gene_filter()

    def _add_gene_chip(self, gene_name):
        chip = QWidget()
        chip.setStyleSheet("QWidget { background: palette(midlight); border-radius: 8px; }")
        chip_layout = QHBoxLayout(chip)
        chip_layout.setContentsMargins(8, 2, 4, 2)
        chip_layout.setSpacing(2)
        chip_layout.addWidget(QLabel(gene_name))
        remove_btn = QToolButton()
        remove_btn.setText("×")
        remove_btn.setAutoRaise(True)
        remove_btn.setCursor(Qt.PointingHandCursor)
        remove_btn.setToolTip(t("Remove {gene} from the filter", gene=gene_name))
        remove_btn.clicked.connect(lambda: self._on_remove_quick_gene(gene_name))
        chip_layout.addWidget(remove_btn)
        self.quick_gene_chips_layout.addWidget(chip)
        self._quick_gene_chip_widgets[gene_name] = chip

    def _on_remove_quick_gene(self, gene_name):
        self._quick_genes = [g for g in self._quick_genes if g != gene_name]
        widget = self._quick_gene_chip_widgets.pop(gene_name, None)
        if widget is not None:
            self.quick_gene_chips_layout.removeWidget(widget)
            # hide() takes the chip off screen immediately; deleteLater() frees
            # it once control returns to the Qt event loop.
            widget.hide()
            widget.deleteLater()
        self._apply_quick_gene_filter()

    def _apply_quick_gene_filter(self):
        self.proxy_model.set_quick_gene_filters(self._quick_genes)

    def _on_proxy_changed(self, *_):
        """Fires whenever the FULL filtered/sorted set changes shape (a new run,
        Clear results, a gene/threshold filter added or removed, gene-set import/
        clear) -- NOT on plain page navigation, which never touches proxy_model.
        Always jumps back to page 1: "page 3 of the old filter" has no meaning
        once the filter itself changed."""
        self._refresh_page(reset_to_first_page=True)
        self._update_no_results_message()

    def _refresh_page(self, reset_to_first_page):
        """Slices exactly one page's worth of rows out of the full filtered/
        sorted proxy_model and copies just those into page_table_model, which is
        what the QTableView actually renders. This is the crux of the fix for
        the reported freeze: no matter how large the filtered result set is
        (e.g. all ~120k rows after clearing a gene filter), the view is never
        asked to render more than RESULTS_PAGE_SIZE rows at once."""
        if reset_to_first_page:
            self._current_page = 0

        total = self.proxy_model.rowCount()
        page_count = max(1, -(-total // RESULTS_PAGE_SIZE))
        self._current_page = max(0, min(self._current_page, page_count - 1))

        start = self._current_page * RESULTS_PAGE_SIZE
        end = min(start + RESULTS_PAGE_SIZE, total)
        page_rows = [self.proxy_model.row_at(r) for r in range(start, end)]
        self.page_table_model.set_rows(page_rows)

        if total == 0:
            self.row_count_label.setText(t("Showing 0 of 0 rows"))
        else:
            self.row_count_label.setText(t("Showing {first:,}–{last:,} of {total:,} rows", first=start + 1, last=end,
                                           total=total))
        self.page_indicator_label.setText(t("Page {page} of {pages}", page=self._current_page + 1, pages=page_count))
        self.prev_page_button.setEnabled(self._current_page > 0)
        self.next_page_button.setEnabled(self._current_page < page_count - 1)

    def _on_prev_page(self):
        if self._current_page > 0:
            self._current_page -= 1
            self._refresh_page(reset_to_first_page=False)

    def _on_next_page(self):
        self._current_page += 1
        self._refresh_page(reset_to_first_page=False)

    def _on_sort_indicator_changed(self, column, order):
        """Sorts the full result set -- never just whatever ~100 rows happen to
        be in page_table_model at click time -- in pure Python (a single
        native Timsort pass over self._all_rows) rather than through
        proxy_model's own sort()/lessThan(), which is far too slow at 100k+
        rows (see ScoreFilterProxyModel's class docstring). table_model.set_rows()
        triggers proxy_model's modelReset, which _on_proxy_changed is already
        connected to, so that alone re-pages from page 1 -- no separate
        _refresh_page() call needed here."""
        _, attr, kind = COLUMNS[column]
        self._all_rows = sorted(self._all_rows, key=lambda r: sort_key_value(r, attr, kind), reverse=order == Qt.DescendingOrder)
        self.table_model.set_rows(self._all_rows)

    def _update_no_results_message(self):
        total = self.table_model.rowCount()
        visible = self.proxy_model.rowCount()
        if total > 0 and visible == 0 and self._quick_genes:
            genes_str = ", ".join(self._quick_genes)
            self.no_results_label.setText(t("No variants found for: {genes}", genes=genes_str))
            self.no_results_label.setVisible(True)
        else:
            self.no_results_label.setText("")
            self.no_results_label.setVisible(False)

    def _on_results_context_menu(self, pos):
        index = self.table_view.indexAt(pos)
        if not index.isValid():
            return
        row = self.page_table_model.row_at(index.row())
        menu = QMenu(self)
        copy_action = menu.addAction(t("Copy variant ID ({id})", id=row.variant_id))
        effect_action = menu.addAction(t("Analyse transcript effect..."))
        chosen = menu.exec(self.table_view.viewport().mapToGlobal(pos))
        if chosen == copy_action:
            QApplication.clipboard().setText(row.variant_id)
        elif chosen == effect_action:
            selected = self.selected_result_rows()
            self.analyse_transcript_effect(selected if row.variant_id in {r.variant_id for r in selected} else [row])

    def selected_result_rows(self):
        """The variants of the rows selected in the results table, one row each (a variant scored in two genes
        is one variant), in table order; rows that weren't scored are left out."""
        rows, seen = [], set()
        for index in sorted(self.table_view.selectionModel().selectedRows(), key=lambda i: i.row()):
            r = self.page_table_model.row_at(index.row())
            if r.variant_id not in seen and not (r.source or "").startswith("skipped"):
                seen.add(r.variant_id)
                rows.append(r)
        return rows

    def analyse_transcript_effect(self, rows=None, dialog_hook=None):
        """Spec 16 item 4: the confirmation (TranscriptEffectDialog), then the Sequence-to-cDNA tab runs the
        variant(s). dialog_hook(dialog) answers the dialog instead of a person (tests); it returns True for
        Compute."""
        rows = self.selected_result_rows() if rows is None else rows
        if not rows:
            QMessageBox.information(self, t("Analyse transcript effect"), t("Select a variant in the results first."))
            return None
        if len(rows) > vcf_bridge.MAX_QUEUE:
            QMessageBox.warning(self, t("Analyse transcript effect"),
                                t("{n} variants are selected; at most {most} can be queued.", n=len(rows),
                                  most=vcf_bridge.MAX_QUEUE))
            return None
        last = getattr(self, "_last_run_settings", None) or {}
        build, mode = last.get("build") or self.build_combo.currentText(), last.get("mode") or self.mode_combo.currentText()
        tab = self.cdna_tab
        try:
            table, fasta = tab._table(build), tab._fasta(build)
        except ValueError as exc:
            QMessageBox.warning(self, t("Analyse transcript effect"), str(exc))
            return None
        bare = lambda chrom: chrom[3:] if chrom.lower().startswith("chr") else chrom

        def genes_of(row):
            names = sorted({r.gene for r in self._all_rows if r.variant_id == row.variant_id and r.gene} | {row.gene})
            found, missing = [], []
            for name in names:
                hits = [g for g in table.find(name) if g.chrom == bare(row.chrom).upper()]
                if hits:
                    found += [g for g in hits if all(g is not f for f in found)]
                else:
                    missing.append(name)
            return found, missing

        def others_of(row, gene):
            # By position as well as by name: a row can carry another symbol of the gene (precomputed data,
            # older annotations) or a neighbouring gene's name only.
            out, seen = [], {row.variant_id}
            for r in self._all_rows:
                inside = bare(r.chrom).upper() == gene.chrom and gene.mask[0] <= r.pos <= gene.mask[1]
                if (r.gene == gene.name or inside) and r.variant_id not in seen \
                        and not (r.source or "").startswith("skipped"):
                    seen.add(r.variant_id)
                    out.append(r)
            return sorted(out, key=lambda r: r.pos)

        vcf = self.vcf_text.toPlainText()
        fetch = lambda chrom, a, b: fasta["chr" + bare(chrom)][a - 1:b].seq.upper()
        normalize = lambda chrom, pos, ref, alt: normalize_variant(fetch, chrom, pos, ref, alt)
        reader = self.refdb.reader(build)

        def estimate_of(gene, edits):
            seconds = haplotype.estimate_seconds([(0, 2 * haplotype.CONTEXT)] * edits) + 15
            if reader is not None and reader.plan(gene.name) is not None:
                return seconds, t("reference values from the database")
            bases = gene.mask[1] - gene.mask[0] + 1 + 2 * refdb.MARGIN
            return seconds + live_seconds(bases), t("the reference computed live first, {bases:,} bases", bases=bases)

        dialog = TranscriptEffectDialog(self, build, rows, genes_of, others_of,
                                        lambda r: vcf_bridge.genotype(vcf, r.chrom, r.pos, r.ref, r.alt, normalize),
                                        estimate_of, mode)
        accepted = dialog_hook(dialog) if dialog_hook else dialog.exec() == QDialog.Accepted
        if not accepted:
            return dialog
        self.tabs.setCurrentWidget(tab)
        tab.start_from_vcf(dialog.requests(), dialog.folder())
        return dialog

    def _on_clear_results(self):
        """Empties the results table and resets every filter (quick gene chips,
        gene-set import, score threshold). Independent of the VCF input section --
        does not touch whatever's currently loaded there. Confirms first since this
        can discard a run that took a while to produce."""
        if self.table_model.rowCount() > 0:
            reply = QMessageBox.question(
                self, t("Clear results?"),
                t("This will discard the current results table. Continue?"),
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return

        self._all_rows = []
        self.table_model.set_rows([])
        self.summary_button.setEnabled(False)
        self._last_run_settings = None
        self._last_run_duration = None

        self._quick_genes = []
        for widget in self._quick_gene_chip_widgets.values():
            self.quick_gene_chips_layout.removeWidget(widget)
            widget.hide()
            widget.deleteLater()
        self._quick_gene_chip_widgets = {}
        self.quick_gene_edit.clear()
        self.proxy_model.set_quick_gene_filters([])

        self._on_clear_gene_set()

        for btn in self.threshold_buttons.values():
            btn.setChecked(False)
        self.proxy_model.set_threshold(None)

    def _cleanup_temp_vcf(self):
        if self._temp_vcf_path and os.path.exists(self._temp_vcf_path):
            try:
                os.remove(self._temp_vcf_path)
            except OSError:
                pass
        self._temp_vcf_path = None

    def _on_run_clicked(self):
        # The build is normally selected from the VCF automatically; this only
        # triggers when it was changed by hand to the other one. (The pipeline
        # also stops by itself when most REF alleles don't match the FASTA.)
        detected_build, evidence = detect_vcf_build(self.vcf_text.toPlainText())
        selected_build = self.build_combo.currentText()
        if detected_build and detected_build != selected_build:
            reply = QMessageBox.warning(
                self, t("Genome build mismatch"),
                t("This VCF is {detected} ({evidence}), but {selected} is selected.\n\n"
                  "Scoring it as {selected} gives wrong or empty results. Run anyway?", detected=detected_build,
                  evidence=t(evidence), selected=selected_build),
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return

        self._save_settings()
        self._cleanup_temp_vcf()

        fd, temp_path = tempfile.mkstemp(suffix=".vcf", prefix="spliceai_gui_")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(self.vcf_text.toPlainText())
        self._temp_vcf_path = temp_path

        precomputed_dir = self.precomputed_dir_edit.text().strip() if self.use_precomputed_radio.isChecked() else None

        self.run_button.setEnabled(False)
        self.summary_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_label.setText(t("Starting..."))
        self.table_model.set_rows([])

        self._run_start_time = time.time()
        self._last_run_settings = {
            "build": self.build_combo.currentText(),
            "mode": self.mode_combo.currentText(),
            "used_precomputed": self.use_precomputed_radio.isChecked(),
            "used_snpeff": self.use_snpeff_checkbox.isChecked(),
            "score_genes": list(self._score_genes),
        }

        self._worker = PipelineWorker(
            vcf_path=temp_path,
            build=self.build_combo.currentText(),
            mode=self.mode_combo.currentText(),
            fasta_path=self._selected_fasta(),
            precomputed_dir=precomputed_dir,
            skip_precomputed=self.skip_precomputed_checkbox.isChecked(),
            use_snpeff=self.use_snpeff_checkbox.isChecked(),
            snpeff_dir=self.snpeff_dir_edit.text().strip(),
            mane_dir=self.mane_dir_edit.text().strip(),
            use_mane=self.use_mane_checkbox.isChecked(),
            genes=self._score_genes or None,
        )
        self._worker.progress.connect(self._on_worker_progress)
        self._worker.finished_ok.connect(self._on_worker_finished)
        self._worker.failed.connect(self._on_worker_failed)
        self._worker.cancelled.connect(self._on_worker_cancelled)
        self._worker.start()
        self._set_running(True)

    def _on_worker_progress(self, message, current, total):
        if total:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(current or 0)
        else:
            self.progress_bar.setRange(0, 0)
        if self._worker is not None and self._worker.control.paused:
            message = PAUSED_PREFIX + message
        self.progress_label.setText(message)

    def _on_worker_finished(self, rows):
        self._all_rows = rows
        self.table_model.set_rows(rows)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        self.progress_label.setText(t("Done -- {n} rows.", n=len(rows)))
        self._worker = None
        self._set_running(False)
        if self._run_start_time is not None:
            self._last_run_duration = time.time() - self._run_start_time
            self._run_start_time = None
        self.summary_button.setEnabled(bool(rows))
        self._update_run_enabled()

    def _on_worker_failed(self, message):
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.progress_label.setText(t("Failed."))
        self._worker = None
        self._set_running(False)
        self._update_spliceai_status_indicator()
        self._run_start_time = None
        self._update_run_enabled()
        QMessageBox.critical(self, t("Pipeline error"), message)

    def _on_worker_cancelled(self):
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.progress_label.setText(t("Run ended -- no results."))
        self._worker = None
        self._set_running(False)
        self._run_start_time = None
        self._update_run_enabled()

    def _set_running(self, running):
        """Pause/End and the CPU display are only live during a run."""
        self.pause_button.setText(t("Pause"))
        self.pause_button.setEnabled(running)
        self.end_button.setEnabled(running)
        self.cpu_label.setVisible(running)
        if running:
            self.cpu_label.setText("")
            self._cpu_meter.reset()
            self._cpu_timer.start()
        else:
            self._cpu_timer.stop()

    def _on_pause_clicked(self):
        if self._worker is None:
            return
        control = self._worker.control
        if control.paused:
            control.resume()
            self.pause_button.setText(t("Pause"))
            self.progress_label.setText(self.progress_label.text().replace(PAUSED_PREFIX, ""))
        else:
            control.pause()
            self.pause_button.setText(t("Resume"))
            self.progress_label.setText(PAUSED_PREFIX + self.progress_label.text())

    def _on_end_clicked(self):
        if self._worker is None:
            return
        reply = QMessageBox.question(
            self, t("End the calculation?"),
            t("End this calculation now? The variants scored so far are discarded."),
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply != QMessageBox.Yes or self._worker is None:
            return
        self._worker.stop()
        self.pause_button.setEnabled(False)
        self.end_button.setEnabled(False)
        self.progress_label.setText(t("Ending the calculation..."))

    def _update_cpu_label(self):
        app_pct, pc_pct = self._cpu_meter.sample()
        if app_pct is None:
            return
        text = t("CPU: this program {pct:.0f}%", pct=app_pct)
        if pc_pct is not None:
            text += " · " + t("whole PC {pct:.0f}%", pct=pc_pct)
        if self._run_start_time is not None:
            text += f" · {format_elapsed(time.time() - self._run_start_time)}"
        self.cpu_label.setText(text)

    def _build_summary_text(self):
        """Statistics for the full result set of the most recently completed run
        -- self._all_rows, deliberately NOT the filtered/paginated proxy_model/
        page_table_model, same "full data" principle as Download full results."""
        rows = self._all_rows
        total_rows = len(rows)

        # total_variants counts distinct (chrom, pos, ref, alt) tuples; the row
        # count can be higher (one row per matching gene/source, as the "Result
        # rows" line below says).
        total_variants = len({(r.chrom, r.pos, r.ref, r.alt) for r in rows})

        n_precomputed = sum(1 for r in rows if r.source == "precomputed")
        n_live = sum(1 for r in rows if r.source == "live")
        n_none = sum(1 for r in rows if r.source == "none")
        n_skipped = sum(1 for r in rows if (r.source or "").startswith("skipped"))

        def pct(n):
            return f"{n / total_rows * 100:.1f}%" if total_rows else "0.0%"

        threshold_counts = {
            x: sum(1 for r in rows if r.max_score is not None and r.max_score >= x)
            for x in THRESHOLDS
        }

        settings = self._last_run_settings or {}
        duration_text = (
            format_elapsed(self._last_run_duration)
            if self._last_run_duration is not None else t("unknown")
        )

        lines = [
            t("Total variants processed: {n:,}", n=total_variants),
            t("Result rows: {n:,} (can exceed variant count -- one row per matching gene/source)", n=total_rows),
            "",
            t("Source breakdown (result rows):"),
            "  " + t("Precomputed lookup: {n:,} ({pct})", n=n_precomputed, pct=pct(n_precomputed)),
            "  " + t("Live-scored: {n:,} ({pct})", n=n_live, pct=pct(n_live)),
            "  " + t("No gene overlap: {n:,} ({pct})", n=n_none, pct=pct(n_none)),
            "  " + t("Skipped, not scorable (reason in Source column): {n:,} ({pct})", n=n_skipped, pct=pct(n_skipped)),
            "",
            t("Score thresholds (max_score >=):"),
        ]
        lines += [f"  >= {x}: {threshold_counts[x]:,}" for x in THRESHOLDS]
        lines += [
            "",
            t("Settings used for this run:"),
            "  " + t("Build: {build}", build=settings.get("build", t("unknown"))),
            "  " + t("Mode: {mode}", mode=settings.get("mode", t("unknown"))),
            "  " + (t("Precomputed data: used") if settings.get("used_precomputed") else t("Precomputed data: not used")),
            "  " + (t("SnpEff annotation: enabled") if settings.get("used_snpeff") else t("SnpEff annotation: disabled")),
            "",
            t("Run duration: {time}", time=duration_text),
        ]
        return "\n".join(lines)

    def _on_show_summary(self):
        if not self._all_rows:
            QMessageBox.information(self, t("No results yet"), t("Run the pipeline first."))
            return
        QMessageBox.information(self, t("Run Summary"), self._build_summary_text())

    def _current_column_specs(self):
        """CSV export column headers deliberately use the raw ScoreRow attribute
        names (variant_id, gene, DS_AG, ...), NOT the table's friendly display
        labels ("Variant Id", "Gene", "DS AG", ...) -- the display labels are a
        presentation-only concern, and changing them must never change what a
        script parsing this CSV's header row sees."""
        return [(attr, attr) for attr in self._current_column_order()]

    def _default_export_filename(self):
        """<original input VCF's filename, extension stripped>-calculated.csv,
        e.g. loading my_sample.vcf.gz -> "my_sample-calculated.csv". Falls back
        to a generic name when no VCF file was ever loaded (pasted text, or a
        session that only used "Open calculated data...")."""
        if self._loaded_vcf_basename:
            return f"{self._loaded_vcf_basename}-calculated.csv"
        return "spliceai_results-calculated.csv"

    def _on_export_clicked(self):
        if self.table_model.rowCount() == 0:
            QMessageBox.information(self, t("Nothing to export"), t("Run the pipeline first."))
            return
        path, _ = QFileDialog.getSaveFileName(
            self, t("Export visible rows"), self._default_export_filename(), t("CSV files (*.csv);;TSV files (*.tsv)")
        )
        if not path:
            return
        visible_rows = [self.proxy_model.row_at(r) for r in range(self.proxy_model.rowCount())]
        delimiter = "\t" if path.lower().endswith(".tsv") else ","
        write_rows_ordered(visible_rows, path, self._current_column_specs(), delimiter=delimiter)
        QMessageBox.information(self, t("Export complete"), t("Wrote {n} rows to {path}", n=len(visible_rows), path=path))

    def _on_download_full_results(self):
        """Exports the complete result set from the most recent run, ignoring any
        active gene/threshold filters currently applied to the table view -- unlike
        "Export visible rows...", which intentionally respects them."""
        if not self._all_rows:
            QMessageBox.information(self, t("Nothing to download"), t("Run the pipeline first."))
            return
        path, _ = QFileDialog.getSaveFileName(
            self, t("Download full results"), self._default_export_filename(), t("CSV files (*.csv);;TSV files (*.tsv)")
        )
        if not path:
            return
        delimiter = "\t" if path.lower().endswith(".tsv") else ","
        write_rows_ordered(self._all_rows, path, self._current_column_specs(), delimiter=delimiter)
        QMessageBox.information(self, t("Download complete"), t("Wrote {n} rows (full, unfiltered) to {path}",
                                                                 n=len(self._all_rows), path=path))

    def _on_open_calculated_data(self):
        """Loads a previously saved results file (CLI or GUI output) directly into
        the results table, skipping the VCF/scoring pipeline entirely. Independent
        of the Input VCF section -- doesn't touch it."""
        path, _ = QFileDialog.getOpenFileName(
            self, t("Open calculated data"), "", t("Results files (*.csv *.tsv);;All files (*)")
        )
        if not path:
            return
        # A WGS result takes ~20 s to parse: done off the GUI thread so the
        # window keeps responding, polled from here.
        self.open_calculated_button.setEnabled(False)
        self.progress_bar.setRange(0, 0)
        self.progress_label.setText(t("Loading {name}...", name=os.path.basename(path)))
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        future = executor.submit(read_rows, path)
        executor.shutdown(wait=False)  # the thread ends once read_rows returns
        timer = QTimer(self)

        def poll():
            if not future.done():
                return
            timer.stop()
            timer.deleteLater()
            self.open_calculated_button.setEnabled(True)
            try:
                rows = future.result()
            except ResultsFileError as exc:
                self.progress_bar.setRange(0, 1)
                self.progress_label.setText("")
                QMessageBox.critical(self, t("Could not load file"), f"{path}\n\n{exc}")
                return
            self._show_loaded_rows(path, rows)

        timer.timeout.connect(poll)
        timer.start(200)

    def _show_loaded_rows(self, path, rows):
        self._all_rows = rows
        self.table_model.set_rows(rows)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        self.progress_label.setText(t("Loaded {n} rows from {name}.", n=len(rows), name=os.path.basename(path)))
        # Summary is disabled: a loaded file has no run settings or duration.
        self.summary_button.setEnabled(False)
        self._last_run_settings = None
        self._last_run_duration = None
