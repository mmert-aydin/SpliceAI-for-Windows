"""The viewer (spec 8.2): a window with the exon schematic, the stacked
base-level view and the sequence comparison of one analysis, opened from the
Sequence-to-cDNA tab's results, from a results file, or (per-base view only)
from a CSV the viewer exported. Results files and CSVs can be dropped onto it.

Spec 15: it opens on the exon schematic, which is quick; the stacked view
draws its tiles and the sequence comparison is built only when first shown.
Each view has its "Transcript orientation (5'->3')" switch (on for the
schematic, off for the stacked view). Tooltips are short; a click (or Shift
+ hover in the stacked view) puts the full details in the side panel.

Spec 16: the applied events are summed up at the top (one line each, a click
zooms there); the stacked view opens on the first event with its flanking
exons, Focus on. The toolbar is two rows that hand what doesn't fit to an
overflow menu (») at any width. The side panel is a key/value table in
sections; its width is kept in the settings.
"""
import html
import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton,
    QScrollArea, QSizePolicy, QSplitter, QStyle, QTabWidget, QTextBrowser, QToolBar, QToolButton, QVBoxLayout,
    QWidget,
)

from spliceai_pipeline import bundle, cdna, cdna_analysis
from spliceai_pipeline.view_data import ViewData

from . import config, theme
from .exon_view import ExonView
from .help_button import make_help_button
from .stacked_view import EVENT_COLORS, StackedView, orientation_text, rightwards
from .i18n import t

SUMMARY_LINES = 6       # events listed in the summary card; the rest are in the event list

SCORE_TRACKS_HELP_TEXT = (
    t("Each column is one base; reference and sample share the columns like a multiple alignment: an insertion in "
    "the sample opens columns with a gap in the reference rows, a deletion columns with a gap in the sample's, "
    "and an inverted stretch is shown in reference order with the sample's bases reverse-complemented and marked.\n\n"
    "Score rows show SpliceAI's probability, at every base, of being an acceptor (blue) or a donor (orange). Zoomed "
    "out they are heat strips -- the stronger the colour, the higher the probability, each pixel showing the highest "
    "value of the bases it covers, so no strong site is lost --; zoomed in they become bars, with the site threshold "
    "as a dashed line; at base level the values are written out. The delta row is the sample's value minus the "
    "reference's (up: higher in the sample).\n\n"
    "Transcripts are drawn as exon blocks on a thin intron line whose chevrons point the way the gene is "
    "transcribed, coloured by the reference exon each exon matches; grey hatched blocks are exons with no "
    "reference counterpart. Beneath each, its own scores at its exonic bases. Sample bases taken from the "
    "reference (where the pasted sequence didn't reach) are grey; inserted bases have a yellow and inverted ones "
    "a lilac background. Hatched score columns have no score (or were scored with N as context near a "
    "chromosome end).\n\n"
    "\"Transcript orientation (5'→3')\" draws a minus-strand gene reverse-complemented, exon 1 on the left; "
    "the line above the view says which orientation is shown. Positions in tooltips are always genomic.\n\n"
    "Applied events are bands through every row in their type's colour (see the legend), labelled above the ruler "
    "with type, size, coordinates and place in the canonical transcript; an insertion is a line in the reference "
    "rows and the inserted bases in the sample's. \"Show origin\" outlines, dashed, where inserted bases come "
    "from; \"Go to origin\" goes there. Focus shows only the reference, WT predicted 1 and the chosen sample rows; "
    "a click on a sample's name lists the other candidates. The delta row is scaled to the largest |delta| in "
    "view; red ticks mark |delta| ≥ 0.1.\n\n"
    "Mouse wheel: zoom (from the whole gene to single bases); drag: pan; double-click: zoom in; the button at the "
    "bottom left: the whole gene and back; hover: a short tooltip; click (or Shift + hover): everything about the "
    "base or exon, in the panel on the right. "
    "Annotated transcripts are drawn with their UTRs thinner than the coding part.")
)


class _CornerScroll(QScrollArea):
    """A scroll area with a button kept at its bottom left, over the view."""
    corner = None

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.corner is not None:
            self.corner.move(8, self.viewport().height() - self.corner.sizeHint().height() - 8)
            self.corner.raise_()


class ViewerWindow(QMainWindow):
    def __init__(self, source, parent=None, candidate=1):
        """source: a cdna_analysis.Analysis, or the path of a results file or viewer CSV."""
        super().__init__(parent)
        self.resize(1300, 850)
        self.setAcceptDrops(True)
        self.load(source, candidate)

    # --- what is shown ---------------------------------------------------------------------
    def load(self, source, candidate=1):
        if isinstance(source, str):
            if source.lower().endswith(bundle.EXTENSION):
                result, _manifest, _notes = bundle.load(source)
                data = ViewData.from_analysis(result)
            else:
                data = ViewData.from_csv(source)
        else:
            data = ViewData.from_analysis(source)
        self.data = data
        self._event_index = None
        self._candidate = candidate
        strand = t("+ strand") if data.strand == "+" else t("− strand")
        self.setWindowTitle(t("Viewer") + f" -- {data.gene} ({data.build}, {data.chrom}, {strand})")
        central = QWidget()
        layout = QVBoxLayout(central)
        if data.limited:
            banner = QLabel(data.limited)
            banner.setWordWrap(True)
            banner.setStyleSheet(theme.STATUS_WARN)
            layout.addWidget(banner)
        self.summary = self._summary_card()
        layout.addWidget(self.summary)
        for bar in self._toolbars():
            layout.addWidget(bar)
        self.tabs = QTabWidget()
        result = getattr(data, "result", None)
        # 1. The exon schematic (spec 15 item 8: first, and where the viewer opens).
        self.exons = ExonView(data)
        self.exons.exon_clicked.connect(self._jump_to_columns)
        self.exons.event_clicked.connect(self._on_event_index)
        self.exons.details.connect(self.show_details)
        page = QWidget()
        pl = QVBoxLayout(page)
        pl.setContentsMargins(0, 0, 0, 0)
        self.exon_orientation = QCheckBox(t("Transcript orientation (5'→3')"))
        self.exon_orientation.setChecked(True)
        self.exon_orientation.toggled.connect(self._on_exon_orientation)
        pl.addWidget(self.exon_orientation)
        exon_scroll = QScrollArea()
        exon_scroll.setWidgetResizable(True)
        exon_scroll.setWidget(self.exons)
        pl.addWidget(exon_scroll)
        self.tabs.addTab(page, t("Exon schematic"))
        # 2. The stacked view: tiles are drawn only once it is shown.
        self.stacked = StackedView(data)
        self.stacked.details.connect(self.show_details)
        self.stacked.focus = self.focus.isChecked()
        self.stacked.sample_rank = self._candidate or 1
        self.stacked.start_view = self.start_view()
        page = QWidget()
        pl = QVBoxLayout(page)
        pl.setContentsMargins(0, 0, 0, 0)
        header = QHBoxLayout()
        pl.addLayout(header)
        self.stacked_orientation = QCheckBox(t("Transcript orientation (5'→3')"))
        self.stacked_orientation.setChecked(False)
        self.stacked_orientation.setEnabled(data.strand == "-")
        self.stacked_orientation.setToolTip(t("A minus-strand gene drawn reverse-complemented, exon 1 on the left.")
                                            if data.strand == "-" else t("A plus-strand gene reads 5'→3' already."))
        self.stacked_orientation.toggled.connect(self._on_stacked_orientation)
        header.addWidget(self.stacked_orientation)
        # Spec 16: the orientation is said here, not on the ruler (where it covered the coordinates).
        self.orientation_label = QLabel()
        self.orientation_label.setStyleSheet(f"color: {theme.INK_2};")
        self.orientation_label.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)   # never widens the view
        header.addWidget(self.orientation_label, stretch=1)
        self._update_orientation_label()
        hint = QLabel(t("Wheel: zoom · drag: pan · double-click: zoom in · Shift + hover: details on the right"))
        hint.setStyleSheet(theme.STATUS_MUTED)
        hint.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        pl.addWidget(hint)
        scroll = _CornerScroll()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.stacked)
        self.whole_button = QToolButton(scroll)
        self.whole_button.setIcon(self.style().standardIcon(QStyle.SP_TitleBarMaxButton))
        self.whole_button.setToolTip(t("Whole gene, or back to the previous view"))
        self.whole_button.setAutoRaise(False)
        self.whole_button.clicked.connect(self._toggle_whole)
        scroll.corner = self.whole_button
        pl.addWidget(scroll)
        self.tabs.addTab(page, t("Stacked view"))
        self.stacked.sample_changed.connect(self.set_candidate)
        self.stacked.moved.connect(self._update_whole_icon)
        # 3. The sequence comparison: built when first shown.
        self.sequences = None
        if result is not None and result.transcripts:
            self._sequence_page = QWidget()
            QVBoxLayout(self._sequence_page).setContentsMargins(0, 0, 0, 0)
            self.tabs.addTab(self._sequence_page, t("Sequence comparison"))
        self.tabs.setTabEnabled(0, not data.limited)
        self.tabs.setCurrentIndex(1 if data.limited else 0)
        self.tabs.currentChanged.connect(self._on_tab)
        self.details_panel = QTextBrowser()
        self.details_panel.setPlaceholderText(t("Click a base, an exon or an event."))      # one line: Qt shows no more
        self.details_panel.setToolTip(t("Everything about what you clicked shows here (Shift + hover in the stacked "
                                      "view does it too)."))
        self.splitter = QSplitter(Qt.Horizontal)
        self.splitter.addWidget(self.tabs)
        side = QWidget()                          # the panel, headed by the gene's name
        sl = QVBoxLayout(side)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(4)
        self.gene_heading = QLabel(data.gene)
        self.gene_heading.setStyleSheet(f"font-size: 14pt; font-weight: 600; font-style: italic; color: {theme.INK}; "
                                        f"padding: 2px 4px;")
        self.gene_heading.setTextInteractionFlags(Qt.TextSelectableByMouse)
        sl.addWidget(self.gene_heading)
        sl.addWidget(self.details_panel, stretch=1)
        self.splitter.addWidget(side)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setChildrenCollapsible(False)
        self.details_panel.setMinimumWidth(160)
        layout.addWidget(self.splitter, stretch=1)
        if self.isVisible():
            self._size_panel()
        self.setCentralWidget(central)
        self.stacked.setMinimumHeight(self.stacked.total_height())
        self._update_status()
        self.stacked.hovered.connect(lambda _c: self._update_status())
        self.stacked.moved.connect(self._update_status)

    def _size_panel(self):
        """The side panel at the width it was last left at (sizes set before the window is laid out are lost)."""
        panel = min(int(config.load()["viewer_details_width"]), self.splitter.width() // 3)     # never at the views' cost
        self.splitter.setSizes([self.splitter.width() - panel, panel])

    def showEvent(self, event):
        super().showEvent(event)
        if not getattr(self, "_panel_sized", False):
            self._panel_sized = True
            self._size_panel()

    def _on_tab(self, index):
        if self.tabs.widget(index) is getattr(self, "_sequence_page", None) and self.sequences is None:
            self.build_sequences()
        self._update_status()

    def build_sequences(self):
        from .sequence_view import SequenceView

        self.sequences = SequenceView(self.data.result)
        self.sequences.set_candidate(self._candidate)
        self.sequences.candidate.currentIndexChanged.connect(lambda i: self.set_candidate(i + 1))
        self._sequence_page.layout().addWidget(self.sequences)
        return self.sequences

    def set_candidate(self, candidate):
        """The sample candidate shown, the same in the stacked view and the sequence comparison."""
        self._candidate = candidate
        self.stacked.set_sample(candidate)
        if self.sequences is not None:
            self.sequences.set_candidate(candidate)

    def _toggle_whole(self):
        self.stacked.toggle_whole()
        self._update_whole_icon()

    def _update_whole_icon(self):
        back = self.stacked.is_whole() and self.stacked._before_whole is not None
        self.whole_button.setIcon(self.style().standardIcon(QStyle.SP_TitleBarNormalButton if back
                                                            else QStyle.SP_TitleBarMaxButton))

    def show_details(self, content):
        """Plain text, or [(section, [(key, value)])] (the stacked view's) as a compact two-column table."""
        if isinstance(content, str):
            self.details_panel.setPlainText(content)
            return
        cells = []
        for title, pairs in content:
            cells.append(f"<tr><td colspan='2' style='padding-top: 8px;'><b>{html.escape(t(title))}</b></td></tr>")
            cells += [f"<tr><td nowrap valign='top' style='color: {theme.INK_3}; padding-right: 8px;'>"
                      f"{html.escape(t(str(k)))}</td><td valign='top'>{html.escape(t(str(v)))}</td></tr>" for k, v in pairs]
        self.details_panel.setHtml("<table width='100%' cellspacing='0' cellpadding='2'>" + "".join(cells) + "</table>")

    def _summary_card(self):
        """Spec 16 item 1: one line per applied event -- type, size, coordinates, place in the canonical
        transcript, nearest sites, origin, what sample 1 does there; a click zooms there."""
        card = QLabel()
        card.setWordWrap(True)
        card.setTextFormat(Qt.RichText)
        card.setTextInteractionFlags(Qt.TextBrowserInteraction)
        card.linkActivated.connect(lambda link: self._on_event_index(int(link.split(":")[1])))
        info = getattr(self.data, "event_info", [])
        lines = [f"<span style='color: {EVENT_COLORS.get(x['kind'], '#555')};'>■</span> "
                 f"<a href='event:{i}'>{html.escape(x['label'].split(' · ')[0])}</a> · "
                 f"{html.escape(x['line'].split(' · ', 1)[1])}" for i, x in enumerate(info[:SUMMARY_LINES])]
        if len(info) > SUMMARY_LINES:
            lines.append("... " + t("and {n} more (in the event list below).", n=len(info) - SUMMARY_LINES))
        card.setText(t("<b>Applied events</b> (click one to zoom there):") + "<br>" + "<br>".join(lines))
        card.setStyleSheet(f"background: {theme.SURFACE}; border: 1px solid #d5dde0; padding: 4px;")
        card.setVisible(bool(info))
        self.summary_lines = [x["line"] for x in info]
        return card

    def _toolbars(self):
        """Two rows (spec 16 item 3); a QToolBar puts what doesn't fit in its overflow menu (»)."""
        d = self.data
        info = getattr(d, "event_info", [])
        nav, view = QToolBar(t("Navigate")), QToolBar(t("Show and export"))
        for bar in (nav, view):
            bar.setMovable(False)
        self.event_combo = QComboBox()
        self.event_combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.event_combo.setMinimumContentsLength(20)
        names = [f"{i}. {x['label'].split(' · ')[0]} · {x['where']}" for i, x in enumerate(info, 1)] or \
            [f"{i}. {text}" for i, (_, _, text) in enumerate(d.events, 1)]
        self.event_combo.addItems(names or [t("No applied events")])
        self.event_combo.setEnabled(bool(d.events))
        self.event_combo.setToolTip(t("The applied events; choosing one zooms to it."))
        self.event_combo.activated.connect(self._on_event)
        nav.addWidget(self.event_combo)
        self.zoom_event_button = QPushButton(t("Zoom to event"))
        self.zoom_event_button.clicked.connect(lambda: self._on_event(self.event_combo.currentIndex()))
        self.zoom_event_button.setEnabled(bool(d.events))
        nav.addWidget(self.zoom_event_button)
        self.origin_button = QPushButton(t("Show origin"))
        self.origin_button.setCheckable(True)
        self.origin_button.setToolTip(t("Outlines (dashed) where the inserted bases of insertions, duplications, "
                                      "inversions and replacements come from."))
        self.origin_button.toggled.connect(self._on_show_origin)
        self.goto_origin_button = QPushButton(t("Go to origin"))
        self.goto_origin_button.setToolTip(t("Zooms to where the chosen event's inserted bases come from."))
        self.goto_origin_button.clicked.connect(self.go_to_origin)
        for b in (self.origin_button, self.goto_origin_button):
            b.setEnabled(any(x["sources"] for x in info))
            nav.addWidget(b)
        nav.addSeparator()
        prev_button = QPushButton(t("◄ Previous"))
        prev_button.setToolTip(t("The previous difference"))
        prev_button.clicked.connect(lambda: self.step_difference(-1))
        nav.addWidget(prev_button)
        next_button = QPushButton(t("Next ►"))
        next_button.setToolTip(t("The next difference"))
        next_button.clicked.connect(lambda: self.step_difference(1))
        nav.addWidget(next_button)
        nav.addSeparator()
        self.exon_combo = QComboBox()
        self.exon_combo.addItem(t("Zoom to exon..."))
        ref_row = next((r for r in d.rows if r.group == "reference"), None)
        self._ref_blocks = sorted(ref_row.blocks, key=lambda b: b.number) if ref_row else []
        self.exon_combo.addItems([f"Exon {b.number}" for b in self._ref_blocks])
        self.exon_combo.activated.connect(self._on_exon)
        nav.addWidget(self.exon_combo)
        whole = QPushButton(t("Whole gene"))
        whole.clicked.connect(lambda: self.stacked.fit())
        nav.addWidget(whole)
        self.jump_edit = QLineEdit()
        self.jump_edit.setPlaceholderText(t("Go to: chr7:117,639,961 or c.3718-2477"))
        self.jump_edit.setMinimumWidth(170)
        self.jump_edit.returnPressed.connect(lambda: self.jump(self.jump_edit.text()))
        view.addWidget(self.jump_edit)
        view.addSeparator()
        self.focus = QCheckBox(t("Focus"))
        self.focus.setToolTip(t("Stacked view: only Reference DNA and scores, the reference transcript, WT predicted 1, "
                              "the sample's sequence, scores and delta, and sample 1."))
        self.focus.setChecked(bool(d.events))
        self.focus.toggled.connect(self._on_focus)
        view.addWidget(self.focus)
        self.alternatives = QCheckBox(t("Annotated alternatives"))
        self.alternatives.setToolTip(t("The gene's other curated transcripts (NM_), from the annotation."))
        self.alternatives.toggled.connect(self._on_alternatives)
        view.addWidget(self.alternatives)
        self.spliceai_alternatives = QCheckBox(t("SpliceAI alternatives"))
        self.spliceai_alternatives.setToolTip(t("The WT and sample candidates 2 to 5."))
        self.spliceai_alternatives.toggled.connect(self._on_spliceai_alternatives)
        view.addWidget(self.spliceai_alternatives)
        view.addSeparator()
        for text, fn in (("PNG...", self._on_png), ("SVG...", self._on_svg), ("CSV (visible)...", self._on_csv)):
            button = QPushButton(text)
            button.clicked.connect(fn)
            view.addWidget(button)
        view.addWidget(make_help_button(self, t("Score tracks in the viewer"), SCORE_TRACKS_HELP_TEXT))
        self.toolbars = (nav, view)
        return nav, view

    def start_view(self):
        """Spec 16 item 1: genomic columns to open the stacked view on -- centred on the first event, its
        flanking canonical exons in view; None (the whole gene) without events."""
        d = self.data
        info = getattr(d, "event_info", [])
        if not d.events or not info or not hasattr(d, "result"):
            return None
        first, last, _ = d.events[0]
        e = info[0]["event"]
        lo, hi = (e.end, e.start) if e.end < e.start else (e.start, e.end)
        k1, k2 = cdna_analysis.flanking_exons(d.result, lo, hi)
        blocks = [b for b in self._ref_blocks if k1 <= b.number <= k2]
        a = min([first] + [b.first for b in blocks])
        b = max([last] + [b.last for b in blocks])
        centre = (first + last) / 2
        half = max(centre - a, b - centre, 30)
        return max(0, int(centre - half)), min(d.n - 1, int(centre + half + 1))

    def _update_status(self):
        d = self.stacked.data
        a, b = self.stacked.visible()
        bad = self.data.unscored_transcribed()
        text = t("{n:,} columns; showing {first} .. {last}. Reference values: {source}.", n=d.n, first=d.label(a),
                 last=d.label(b), source=t(d.source))
        if bad:
            text += " " + t("{n:,} transcribed position(s) not scored the usual way: ",
                            n=sum(y - x + 1 for x, y, _ in bad)) + "; ".join(sorted({why for _, _, why in bad})) + "."
        self.statusBar().showMessage(text)

    # --- orientation ------------------------------------------------------------------------------
    def _on_exon_orientation(self, on):
        self.exons.transcript_orientation = on
        self.exons.update()

    def _on_stacked_orientation(self, on):
        self.stacked.set_transcript_orientation(on)
        self._update_orientation_label()
        self._update_status()

    def _update_orientation_label(self):
        d = self.stacked.data
        arrow = "5' ──▶ 3'" if rightwards(d) else "3' ◀── 5'"
        self.orientation_label.setText(f"{orientation_text(d, t)}    {arrow}")
        self.orientation_label.setToolTip(self.orientation_label.text())

    # --- navigation -----------------------------------------------------------------------------
    def _jump_to_columns(self, first, last):
        """Shows genomic columns first..last in the stacked view."""
        self.tabs.setCurrentIndex(1)
        self.stacked.show_genomic(first, last)
        self._update_status()

    def _on_event(self, index):
        if 0 <= index < len(self.data.events):
            self._on_event_index(index)

    def _on_event_index(self, i):
        first, last, _ = self.data.events[i]
        self._event_index = sorted(self.data.events).index(self.data.events[i])
        self.event_combo.setCurrentIndex(i)
        self._jump_to_columns(first, last)

    def _on_show_origin(self, on):
        self.stacked.show_origin = on
        self.stacked.update()

    def go_to_origin(self):
        """Zooms to where the chosen event's inserted bases come from (or the first event that has an origin),
        outlined."""
        info = getattr(self.data, "event_info", [])
        i = self.event_combo.currentIndex()
        chosen = info[i] if 0 <= i < len(info) and info[i]["sources"] else next((x for x in info if x["sources"]),
                                                                               None)
        if chosen is None:
            return False
        spans = [self.data.columns_of(r0, r1) for r0, r1, _ in chosen["sources"]]
        self.origin_button.setChecked(True)
        self._jump_to_columns(min(a for a, _ in spans), max(b for _, b in spans))
        return True

    def _on_focus(self, on):
        self.stacked.focus = on
        self.stacked.setMinimumHeight(self.stacked.total_height())
        self.stacked.update()

    def _on_exon(self, index):
        if index > 0:
            b = self._ref_blocks[index - 1]
            self._jump_to_columns(b.first, b.last)

    def step_difference(self, direction):
        """Zooms to the next (direction 1) or previous (-1) difference, in order
        along the gene, starting from the first; it goes round at the ends."""
        events = sorted(self.data.events)
        if not events:
            return
        current = getattr(self, "_event_index", None)
        if current is None:
            current = -1 if direction > 0 else 0
        self._event_index = (current + direction) % len(events)
        first, last, _ = events[self._event_index]
        self._jump_to_columns(first, last)

    def jump(self, text):
        """Centres the view on a genomic position (chr:pos or a plain number) or a c. position."""
        text = text.strip().replace(",", "")
        try:
            if text.lower().startswith("c."):
                result = getattr(self.data, "result", None)
                if result is None or not result.transcripts:
                    raise ValueError("c. positions need a results file, not a CSV")
                pos = cdna.c_to_genomic(result.transcripts[0].cdna, text)
            else:
                pos = int(text.split(":")[-1].split("-")[0])
        except ValueError as exc:
            QMessageBox.warning(self, t("Go to"), t("Can't go to {text}: {error}", text=repr(text), error=t(str(exc))))
            return False
        c = self.data.column_of_ref(pos)
        self._jump_to_columns(max(0, c - 40), c + 40)
        self.stacked.cursor_col = self.stacked.to_display(c)
        return True

    def _on_alternatives(self, on):
        if on:
            self.focus.setChecked(False)        # what was asked for is shown
        self.stacked.show_alternatives = self.exons.show_alternatives = on
        self.stacked.setMinimumHeight(self.stacked.total_height())
        self.stacked.update()
        self.exons.update()

    def _on_spliceai_alternatives(self, on):
        if on:
            self.focus.setChecked(False)
        self.stacked.show_spliceai_alternatives = self.exons.show_spliceai_alternatives = on
        self.stacked.setMinimumHeight(self.stacked.total_height())
        self.stacked.update()
        self.exons.update()

    # --- exports --------------------------------------------------------------------------------------
    def _name(self):
        return f"{self.data.gene}_{self.data.build}_viewer"

    def _on_png(self):
        """The view shown (the exon schematic or the stacked view), at twice the screen's resolution."""
        view = self.exons if self.tabs.currentWidget() is self.tabs.widget(0) else self.stacked
        path, _ = QFileDialog.getSaveFileName(self, t("Save as PNG"), self._name() + ".png", t("PNG (*.png)"))
        if path:
            view.render_png(path)

    def _on_svg(self):
        path, _ = QFileDialog.getSaveFileName(self, t("Save as SVG"), self._name() + ".svg", t("SVG (*.svg)"))
        if path:
            self.stacked.render_svg(path)

    def _on_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, t("Export the visible region as CSV"), self._name() + ".csv",
                                              t("CSV (*.csv)"))
        if path:
            self.data.write_csv(path, *self.stacked.visible_genomic())

    # --- drag and drop --------------------------------------------------------------------------------
    def dragEnterEvent(self, event):
        if any(u.isLocalFile() and u.toLocalFile().lower().endswith((bundle.EXTENSION, ".csv"))
               for u in event.mimeData().urls()):
            event.acceptProposedAction()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile() if url.isLocalFile() else ""
            if path.lower().endswith((bundle.EXTENSION, ".csv")):
                event.acceptProposedAction()
                self.open_file(os.path.normpath(path))
                return

    def open_file(self, path):
        try:
            self.load(path)
        except (bundle.BundleError, ValueError, OSError) as exc:
            QMessageBox.warning(self, t("Open in the viewer"), str(exc))
            return False
        return True

    def closeEvent(self, event):
        width = self.splitter.sizes()[1]
        settings = config.load()
        if width > 0 and width != settings["viewer_details_width"]:
            settings["viewer_details_width"] = width          # spec 16 item 3: kept for the next viewer
            config.save(settings)
        self.stacked.close()
        super().closeEvent(event)
