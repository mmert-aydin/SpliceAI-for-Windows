"""Phase 7: the Sequence-to-cDNA tab in the real window (Qt offscreen,
throw-away home folder, no database, so the reference is computed live): the
window's two tabs, the "no database" panel, choosing build and gene, a
synthetic patient's sequence checked and run, the differences table, the
results table and its warnings and status, CSV / results bundle / bedGraph
exports, the bundle opened again by dropping it onto the window, a gene
without an NM_ (hg19), Cancel, and the help texts.

    .venv\\Scripts\\python tests\\test_cdna_tab.py
"""
import csv
import glob
import json
import os
import random
import sys
import tempfile
import time

import helpers
from helpers import check

helpers.watchdog(1200)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg19",
               "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True,
               "mane_dir": str(helpers.MANE_DIR)}, fh)

from PySide6.QtCore import QMimeData, QPoint, QPointF, Qt, QUrl  # noqa: E402
from PySide6.QtGui import QDragEnterEvent, QDropEvent  # noqa: E402
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QMessageBox, QToolButton  # noqa: E402

from spliceai_pipeline import bundle, cdna_analysis as ca, haplotype  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.tabix import TabixFile  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui import config, theme  # noqa: E402
from spliceai_gui.main_window import MainWindow  # noqa: E402

shown = []
QMessageBox.warning = staticmethod(lambda _p, title, text, *a, **k: shown.append((title, text)) or QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda _p, title, text, *a, **k: shown.append((title, text)) or QMessageBox.Ok)
QMessageBox.critical = staticmethod(lambda _p, title, text, *a, **k: shown.append((title, text)) or QMessageBox.Ok)
QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Yes)
out_dir = str(home / "out")
os.makedirs(out_dir)


def wait(what, until, seconds=240):
    deadline = time.time() + seconds
    while time.time() < deadline:
        app.processEvents()
        if until():
            return True
        time.sleep(0.05)
    check(f"waiting for {what}", False, f"not within {seconds} s")
    return False


w = MainWindow()
tab = w.cdna_tab
check("the window has two tabs: VCF scoring (as before, shown first) and Sequence to cDNA",
      w.tabs.count() == 2 and w.tabs.tabText(0) == "VCF scoring" and w.tabs.tabText(1) == "Sequence to cDNA"
      and w.tabs.widget(0) is w.panes and w.tabs.currentIndex() == 0)
w.tabs.setCurrentWidget(tab)
w.show()
wait("the introduction", lambda: tab.intro is not None and tab.intro.isVisible(), 10)
from spliceai_gui import cdna_help  # noqa: E402

intro = tab.intro.text.toPlainText()
check("no database on this computer (an install from GitHub): the introduction first says how to get one -- "
      "the USB drive's folder, computing it here, or neither", intro.index("Before the First Analysis")
      < intro.index("What It Does") and "From a folder (USB drive)" in intro and "Compute it on this computer" in intro
      and "Continue without it" in intro, intro[:300])
check("introduction (spec 15 item 6): shown the first time the tab is opened, with what the tab does, every "
      "parameter (default, raising, lowering), every result term, the general notes and how to read the result",
      "three kinds of transcript" in intro and "Comparing the WT prediction with the reference" in intro
      and all(name in intro and f"default {default}" in intro for name, default, _, _ in cdna_help.PARAMETERS)
      and all(term in intro for term, _ in cdna_help.RESULT_TERMS)
      and all(note[:60] in intro for note in cdna_help.GENERAL_NOTES)
      and "How to Read the Results" in intro and "at least 0.5" in intro and "confirm the change with RNA" in intro,
      intro[:200])
tab.intro.dont_show.setChecked(True)
tab.intro.close()
check("'Don't show again' is kept in the settings; the '?' in the tab's header opens the introduction again",
      config.load()["cdna_show_intro"] is False and tab.tab_help_button.click() is None and tab.intro.isVisible())
tab.intro.close()
import html  # noqa: E402

results_help = cdna_help.SCORES_HELP_TEXT + cdna_help.CHECKS_HELP_TEXT + cdna_help.HOW_TO_READ_HELP_TEXT
check("the help dialogs use the same texts (written once): Reading the results in three, terms in bold",
      cdna_help.SCORES_HELP_TEXT.startswith("<p><b>Relative score</b>: " + html.escape(cdna_help.RESULT_TERMS[0][1]))
      and "<b>Sequence-specificity control</b>" in cdna_help.CHECKS_HELP_TEXT
      and html.escape(cdna_help.HOW_TO_READ.split(chr(10) * 2)[0]) in cdna_help.HOW_TO_READ_HELP_TEXT
      and cdna_help.GENERAL_NOTES[1] in cdna_help.TAB_HELP_TEXT
      and not any(term in (t for t, _ in cdna_help.RESULT_TERMS) for term in ("In frame", "PTC", "NMD")))
labels = [b.text() for b in tab.findChildren(QToolButton) if b.text().endswith(" (?)")]
check("...three labelled helps next to WT alternatives; the tab's own '?' at the left",
      all(f"{name} (?)" in labels for name in ("Scores and evidence", "Optional checks", "How to read the result"))
      and tab.tab_help_button.x() < 40, f"{labels} x={tab.tab_help_button.x()}")
# items 4 and 7
check("input area (item 4): a short Gene placeholder, the long explanation in its tooltip; the input tabs are "
      "'Paste sequence' and 'From coordinates'", tab.gene_edit.placeholderText() == "Gene symbol, e.g. BRCA2"
      and "RefSeq" in tab.gene_edit.toolTip() and "SpliceAI" in tab.gene_edit.toolTip()
      and [tab.input_tabs.tabText(i) for i in range(2)] == ["Paste sequence", "From coordinates"])


def contrast(a, b):
    lum = lambda c: sum(k * ((v / 255 / 12.92) if v / 255 <= 0.03928 else ((v / 255 + 0.055) / 1.055) ** 2.4)
                        for k, v in zip((0.2126, 0.7152, 0.0722), (int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16))))
    hi, lo = sorted((lum(a), lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


check("'Open viewer' (item 7) in its own accent colour, white text at a contrast of at least 4.5:1, distinct from "
      "the Run button", tab.viewer_button.objectName() == "viewerButton"
      and f"QPushButton#viewerButton {{\n    background: {theme.VIEWER}" in theme.STYLESHEET
      and contrast(theme.VIEWER, "#ffffff") >= 4.5 and theme.VIEWER != theme.TEAL,
      f"{contrast(theme.VIEWER, '#ffffff'):.1f}:1")
wait("the database check", lambda: all(s["state"] != "checking" for s in w.refdb.status.values()), 60)
tab._update_database_view()
check("no database: the tab shows the 'no database yet' panel instead of its inputs",
      not tab.missing_panel.isHidden() and tab.body.isHidden())
next(b for b in tab.missing_panel.findChildren(type(tab.run_button)) if b.text().startswith("Continue")).click()
check("...and 'Continue without it' shows the inputs", tab.missing_panel.isHidden() and not tab.body.isHidden())
check("the genome build isn't preset (chosen anew)", tab.build_combo.currentIndex() == -1 and not tab.check_button.isEnabled())

# --- a patient of MRPL39 (minus strand): a donor +1 SNV of exon 3 and a 70 bp insertion ---------
tab.build_combo.setCurrentText("hg38")
tab.gene_edit.setText("MRPL39")
tab.gene_edit.editingFinished.emit()
wait("the gene", lambda: tab.gene is not None, 60)
check("the gene: its transcript and rule are shown, the reference values said to be computed live with the "
      "time that takes, and 'Build this gene now' offered", "NM_017446.4 (MANE Select)" in tab.gene_label.text()
      and "computed live" in tab.gene_label.text() and "about 0:00:" in tab.gene_label.text()
      and "bases, once" in tab.gene_label.text() and not tab.build_gene_button.isHidden()
      and tab.gene_note.isHidden(), tab.gene_label.text())
check("the gene line gives its strand (spec 15 item 1): MRPL39 is on the − strand",
      "chr21, − strand:" in tab.gene_label.text(), tab.gene_label.text()[:80])
check("the gene: its WT concordance from the genome-wide run is shown (spec 14 C1)",
      "WT concordance for this gene: 10/10 exons, the whole transcript exactly (default parameters)"
      in tab.gene_label.text(), tab.gene_label.text())
from spliceai_pipeline import wt_concordance  # noqa: E402

odd = next(r for r in wt_concordance.read_table(os.path.join(wt_concordance.DATA_DIR, "wt_concordance_hg38.tsv.gz"))
           .values() if r["differing"])
check("the table names the exons predicted differently, one per exon not identical",
      len(odd["differing"]) == odd["exons"] - odd["identical"], str(odd)[:200])
gene = tab.gene
from pyfaidx import Fasta  # noqa: E402

record = Fasta(str(helpers.HG38))["chr" + gene.chrom]
region = (gene.mask[0] - 5000, gene.mask[1] + 5000)
ref = record[region[0] - 1:region[1]].seq.upper()
exons = list(reversed(gene.canonical.exons))
site = exons[2][0] - 1
rng = random.Random(9)
edits = [Edit(site, site, "A" if ref[site - region[0]] != "A" else "C", "SNV"),
         Edit(exons[6][1] + 700, exons[6][1] + 699, "".join(rng.choice("ACGT") for _ in range(70)), "insertion")]
patient, _ = haplotype.build(ref, region[0], edits)
tab.sequence_edit.setPlainText(">patient MRPL39\n" + "\n".join(patient[i:i + 60] for i in range(0, len(patient), 60)))
check("with build, gene and sequence, Check is enabled; Run isn't yet", tab.check_button.isEnabled()
      and not tab.run_button.isEnabled())
tab.check_button.click()
wait("the compatibility check", lambda: tab.compat is not None, 120)
check("check: the report and both differences; the 70 bp insertion applied by default, the SNV not",
      "Identity" in tab.report_label.text() and tab.events_table.rowCount() == 2
      and sorted((tab.events_table.item(i, 1).text(), tab.events_table.item(i, 0).checkState() == Qt.Checked)
                 for i in range(2)) == [("SNV", False), ("insertion (novel sequence)", True)], tab.report_label.text()[:200])
check("check: an estimate of the run time is shown, and Run is enabled", "Estimated run time" in tab.estimate_label.text()
      and tab.run_button.isEnabled(), tab.estimate_label.text())
card = tab.report_label.text()
check("the summary card (spec 17 item 3): identity, the covered region with commas and that the whole transcript is "
      "covered, the orientation as information, the number of differences; no warning for it",
      "Identity" in card and f"chr21:{region[0]:,}-{region[1]:,}; whole transcript covered" in card
      and "MRPL39 is on the − strand; the sequence was given as the reference (+) strand" in card
      and "Differences" in card and "turned around" not in tab.caution_button.toolTip(), card[:600])
check("technical details are collapsed: aligned columns, the method and the sequence positions under 'Details'",
      tab.details_label.isHidden() and "aligned columns" in tab.details_label.text()
      and "of the sequence" in tab.details_label.text() and "aligned columns" not in card, tab.details_label.text()[:200])
tab.details_button.click()
check("...and open with a click", not tab.details_label.isHidden())
tab.details_button.click()
ins_row = next(i for i in range(2) if tab.events_table.item(i, 1).text().startswith("insertion"))
ins = next(e for e in tab.compat.alignment.events if e.kind == "insertion")
head = [tab.events_table.horizontalHeaderItem(i).text() for i in range(tab.events_table.columnCount())]
check("the differences table: Location, shortened Inserted, Origin; the full sequence in the tooltip, in a copy action",
      head == ["Apply", "Type", "Location", "Reference position", "Size", "Inserted", "Origin", "Note"]
      and tab.events_table.item(ins_row, 2).text().startswith("intron 6, c.")
      and tab.events_table.item(ins_row, 5).text() == f"{ins.inserted[:6]}…{ins.inserted[-7:]}, 70 bp"
      and ins.inserted[:60] in tab.events_table.item(ins_row, 5).toolTip()
      and tab.events_table.item(ins_row, 6).text().startswith("novel"),
      [tab.events_table.item(ins_row, c).text() for c in range(1, 7)])
tab.copy_inserted(ins_row)
check("the copy action puts the whole inserted sequence on the clipboard", app.clipboard().text() == ins.inserted)
widest = tab.events_table.fontMetrics().horizontalAdvance("insertion (novel sequence)")
check("the Type column isn't cut off", tab.events_table.columnWidth(1) >= widest,
      f"{tab.events_table.columnWidth(1)} px for {widest} px of text")
snv_row = next(i for i in range(2) if tab.events_table.item(i, 1).text() == "SNV")
tab.events_table.item(snv_row, 0).setCheckState(Qt.Checked)
check("ticking Apply chooses the SNV too", all(e.apply for e in tab.compat.alignment.events))
tab.run_button.click()
check("while running: Cancel on, Run, Check and the differences table off", tab.cancel_button.isEnabled()
      and not tab.run_button.isEnabled() and not tab.events_table.isEnabled())
wait("the run", lambda: tab.result is not None and tab._run_worker is None, 300)
rows = [[tab.results_table.item(r, c).text() for c in range(len(ca.COLUMNS))] for r in range(tab.results_table.rowCount())]
types = [r[0] for r in rows]
check("results: Reference, the wild type's top transcript, the sample's 5 (the wild type's alternatives folded)",
      types == ["Reference", "WT predicted"] + ["Sample"] * 5, str(types))
check("results: the patient's top transcript skips exon 3, briefly with r. and p.",
      rows[2][ca.COLUMNS.index("Change vs Reference")].startswith("Exon 3 skipped")
      and "r.(" in rows[2][ca.COLUMNS.index("Change vs Reference")], rows[2][ca.COLUMNS.index("Change vs Reference")])
ev = ca.COLUMNS.index("Evidence vs WT")
check("results: each patient candidate's evidence against the WT; below the table, how clearly rank 1 leads rank 2",
      rows[2][ev].startswith("max delta") and not rows[1][ev]
      and tab.lead_label.text().startswith("Sample rank 1 leads rank 2 by"), f"{rows[2][ev]} | {tab.lead_label.text()}")
check("results: the local WT concordance around each applied change (spec 14 C2)",
      tab.local_label.text().count("Around the") == 2 and "matches the reference for exons" in tab.local_label.text(),
      tab.local_label.text())
check("results (spec 16): under the table, one line per applied event -- type, size, place in the canonical "
      "transcript, nearest sites, what sample 1 does there; a click opens the viewer there",
      tab.event_lines == [d["line"] for d in ca.event_summary(tab.result)] and len(tab.event_lines) == 2
      and tab.events_label.text().count("href='event:") == 2 and "sample 1:" in tab.events_label.text(),
      tab.events_label.text()[:300])
tab._on_event_link("event:1")
first, last, _ = tab.viewer.data.events[1]
va, vb = tab.viewer.stacked.visible_genomic()
check("...clicking an event line opens the viewer on that event", tab.viewer.tabs.currentIndex() == 1
      and va <= first and last <= vb and tab.viewer.event_combo.currentIndex() == 1)
tab.viewer.close()
tab.viewer = None
check("the weak-evidence setting (0.5) is under Advanced and saved with the run",
      tab.weak_delta.isHidden() and config.load()["cdna_weak_delta"] == 0.5 and tab._settings().weak_delta == 0.5)
tab.show_alternatives.setChecked(True)
check("'WT alternatives' adds them (rank 2 to 5)", tab.results_table.rowCount() == 11)
check("notices (item 5): no warnings banner and no research-use banner; that the reference was computed live is "
      "said in the status line", not hasattr(tab, "warnings_label") and "computed live" in tab.progress_label.text()
      and not any("Research use only" in x.text() for x in tab.findChildren(QLabel)), tab.progress_label.text())
expected = ca.result_notices(tab.result)
icons, wrong = 0, []
for i, r in enumerate(tab.visible_rows()):
    kind = {"Sample": "Patient"}.get(r["Type"], r["Type"])
    want = [t for _, t, where in expected if (kind, r["Rank"] or 0) in where]
    item = tab.results_table.item(i, 0)
    has = not item.icon().isNull()
    icons += has
    if has != bool(want) or (want and not all(t in item.toolTip() for t in want)):
        wrong.append((r["Type"], r["Rank"]))
check("results: a '?' or '!' next to exactly the rows a notice concerns, its text in the tooltip", not wrong,
      f"{icons} icons; {wrong}")
tab._update_status()
check("status: the database's state for the build and the cache in use", "Database hg38: missing"
      in tab.status_label.text() and "cache" in tab.status_label.text() and "Done" in tab.progress_label.text(),
      tab.status_label.text())
check("the sequence-specificity control is off until clicked", not tab.specificity_label.text()
      and tab.specificity_button.isEnabled())
tab.specificity_n.setValue(5)
tab.specificity_button.click()
wait("the specificity control", lambda: tab._run_worker is None, 300)
check("clicked: per change, the observed values next to 5 random stand-ins, said to be no p-value; kept with the "
      "results", "not a p-value" in tab.specificity_label.text() and "random: " in tab.specificity_label.text()
      and "/5 changed" in tab.specificity_label.text() and len(tab.result.specificity) == 2
      and config.load()["cdna_specificity_n"] == 5, tab.specificity_label.text())
tab.specificity_n.setValue(20)

patient_row = types.index("Sample")
tab.show_alternatives.setChecked(False)
rows_now = tab.visible_rows()
third = next(i for i, r in enumerate(rows_now) if r["Type"] == "Sample" and r["Rank"] == 3)
tab.results_table.cellClicked.emit(third, 0)
viewer = tab.viewer
check("clicking a sample row opens the viewer, on the exon schematic", viewer is not None and viewer.isVisible()
      and viewer.data.gene == "MRPL39" and viewer.tabs.currentIndex() == 0)
viewer.tabs.setCurrentIndex(2)
app.processEvents()
check("...and that candidate is chosen in its sequence comparison, built when first shown",
      viewer.sequences is not None and viewer.sequences.candidate.currentIndex() == 2)
check("the viewer gives each candidate's evidence and the lead of rank 1", all(
    "max delta" in viewer.sequences.candidate.itemText(i) or "same as the WT" in viewer.sequences.candidate.itemText(i)
    for i in range(viewer.sequences.candidate.count())), viewer.sequences.candidate.itemText(0))
tab.show_alternatives.setChecked(True)

# --- exports ----------------------------------------------------------------------------------------
csv_path = os.path.join(out_dir, "table.csv")
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (csv_path, ""))
tab.export_csv_button.click()
with open(csv_path, newline="", encoding="utf-8") as fh:
    lines = list(csv.reader(fh))
n_events = len(ca.event_summary(tab.result))
check("Export CSV: the rows shown, with the table's columns (after Excel's sep= line, the research-use note and, "
      "spec 16, one line per applied event with where its inserted bases come from)",
      lines[0] == ["sep=", ""] and lines[1] == [ca.RESEARCH_NOTE]
      and [x[0] for x in lines[2:2 + n_events]] == [f"Event {i}: {d['line']}" for i, d in
                                                     enumerate(ca.event_summary(tab.result), 1)]
      and n_events >= 1 and lines[2 + n_events] == ca.COLUMNS
      and len(lines) == 3 + n_events + tab.results_table.rowCount(), str(lines[:4]))
saved = tab.save_results(os.path.join(out_dir, "patient"))
check("Download results: one bundle file", saved.endswith(bundle.EXTENSION) and os.path.isfile(saved))
before = [[tab.results_table.item(r, c).text() for c in range(len(ca.COLUMNS))] for r in range(tab.results_table.rowCount())]
tab.result = None
tab._fill_table()
mime = QMimeData()
mime.setUrls([QUrl.fromLocalFile(saved)])
QApplication.sendEvent(w, QDragEnterEvent(QPoint(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier))
QApplication.sendEvent(w, QDropEvent(QPointF(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier))
after = [[tab.results_table.item(r, c).text() for c in range(len(ca.COLUMNS))] for r in range(tab.results_table.rowCount())]
check("the bundle dropped onto the window (this tab shown) opens: the same table", tab.result is not None
      and after == before and "Opened" in tab.progress_label.text())
paths = tab.export_bedgraph(os.path.join(out_dir, "tracks"))
with TabixFile(paths[0]) as tf:
    first = next(iter(tf.fetch("chr" + gene.chrom, tab.result.plan["start"], tab.result.plan["end"])), "")
check("Export bedGraph: reference, sample and delta tracks for acceptor and donor, each indexed",
      len(paths) == 6 and all(os.path.exists(p + ".tbi") for p in paths) and first.startswith("chr" + gene.chrom)
      and any("sample_acceptor" in p for p in paths),
      str([os.path.basename(p) for p in paths]))
damaged = os.path.join(out_dir, "not-a-bundle" + bundle.EXTENSION)
with open(damaged, "w") as fh:
    fh.write("nonsense")
shown.clear()
check("a damaged results file is refused with a message", tab.open_results(damaged) is None
      and shown and "isn't a results bundle" in shown[-1][1])

# --- a gene without an NM_ (hg19 BOP1), Cancel, help --------------------------------------------
tab.build_combo.setCurrentText("hg19")
tab.gene_edit.setText("BOP1")
tab.gene_edit.editingFinished.emit()
wait("BOP1", lambda: tab.gene is not None and tab.gene.name == "BOP1", 60)
check("hg19 BOP1: no NM_ there -- explained at the top of the tab, and nothing can be checked or run",
      "no NM_ transcript" in tab.gene_label.text() and not tab.gene_note.isHidden()
      and "Nothing can be predicted" in tab.gene_note.text() and not tab.check_button.isEnabled()
      and not tab.run_button.isEnabled(), tab.gene_note.text()[:200])
tab.gene_edit.setText("RHBG")
tab.gene_edit.editingFinished.emit()
wait("RHBG", lambda: tab.gene is not None and tab.gene.name == "RHBG", 60)
check("hg19 RHBG: the note at the top names the fallback NM_; Check is possible",
      tab.gene_note.text() == "No MANE Select transcript is available for this gene on hg19. The reference cDNA and "
      "the anchors (transcript start, end and ATG) are taken from NM_001412175.1." and tab.check_button.isEnabled(),
      tab.gene_note.text())
# --- Coordinates: no gene entered -> the gene(s) there suggested; a single position searched around (14 C5) ---
tab.build_combo.setCurrentText("hg38")
tab.gene_edit.setText("")
tab.window_edit.setText("chr21:14,533,713")
tab.window_edit.editingFinished.emit()
wait("the genes there", lambda: tab.gene is not None, 60)
check("coordinates with no gene: the genes there are suggested, the MANE Select one chosen and filled in",
      tab.gene.name == "SAMSN1" and tab.gene_edit.text() == "SAMSN1" and tab.gene_label.text().startswith(
          "Genes at chr21:14,533,713: SAMSN1, AF165138.7; SAMSN1 is chosen (change it in the list).")
      and not tab.gene_choice.isHidden(), tab.gene_label.text()[:160])
s_gene = tab.gene
fa = Fasta(str(helpers.HG38))["chr" + s_gene.chrom]
tab.sequence_edit.setPlainText(fa[14533713 - 1500 - 1:14533713 + 1500].seq)
tab.check_button.click()
wait("the check around the position", lambda: tab.compat is not None, 120)
check("a single position: the check searched around it and filled in the range it found",
      tab.window_edit.text() == "chr21:14,532,213-14,535,213" and "Found within" in tab.details_label.text(),
      f"{tab.window_edit.text()} | {tab.details_label.text()[:120]}")

# --- From coordinates (spec 14 C6) -------------------------------------------------------------------
tab.window_edit.setText("")
e5 = list(reversed(s_gene.canonical.exons))[4]
tab.input_tabs.setCurrentIndex(1)
tab.coords_edit.setPlainText(f"# exon 5's acceptor\ndeletion chr21:{e5[1] - 19:,}-{e5[1] + 60:,}")
tab.generate_button.click()
wait("the generated sequence", lambda: tab.compat is not None or "✗" in tab.coords_status.text()
     or "isn't" in tab.coords_status.text(), 120)
check("From coordinates: 'Generate sequence' puts the sample sequence (FASTA) in the paste box and the event into "
      "the differences table directly, applied and noted as given as coordinates",
      tab.compat is not None and tab.sequence_edit.toPlainText().startswith(f">SAMSN1_hg38_sample chr21:")
      and tab.events_table.rowCount() == 1 and "given as coordinates" in tab.events_table.item(0, 7).text()
      and tab.input_tabs.currentIndex() == 1 and "1-based and inclusive" in tab.input_tabs.widget(1).findChild(
          type(tab.coords_status)).text(), f"{tab.coords_status.text()} | {tab.sequence_edit.toPlainText()[:60]}")
fasta_path = os.path.join(str(home), "sample.fa")
QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (fasta_path, ""))
tab.save_fasta_button.click()
saved_fasta = open(fasta_path, encoding="ascii").read().split("\n")
check("'Save as FASTA...' writes that sequence, 60 bases a line", saved_fasta[0].startswith(">SAMSN1_hg38_sample")
      and len(saved_fasta[1]) == 60 and "".join(saved_fasta[1:]) == tab.compat.sequence)
tab.coords_edit.setPlainText(f"deletion chr21:{e5[1] - 19}-{e5[1] - 16} AAAA" if fa[e5[1] - 20:e5[1] - 16].seq
                             .upper() != "AAAA" else f"deletion chr21:{e5[1] - 19}-{e5[1] - 16} CCCC")
tab.generate_button.click()
wait("the refused event", lambda: "not" in tab.coords_status.text(), 60)
check("an event whose given reference bases differ is refused in red, saying what the reference has",
      "the reference at chr21:" in tab.coords_status.text() and tab.coords_status.styleSheet() == theme.STATUS_CRIT,
      tab.coords_status.text())
atg = s_gene.canonical.cds[1]                   # minus strand: the ATG's A is the CDS's highest base
tab.coords_edit.setPlainText(f"deletion chr21:{atg - 2}-{atg + 3}")
tab.generate_button.click()
wait("the generated sequence", lambda: tab.compat is not None, 120)
check("a change touching the ATG is critical (item 5): a red '!' with the text, Run blocked until 'Run anyway'",
      not tab.critical_button.isHidden() and "ATG" in tab.critical_button.toolTip()
      and not tab.run_button.isEnabled() and not tab.run_anyway.isHidden(), tab.check_status.text())
tab.run_anyway.setChecked(True)
check("...ticking 'Run anyway' allows Run", tab.run_button.isEnabled())
tab.coords_edit.setPlainText(f"deletion chr21:{e5[1] - 19}-{e5[1] + 60}")
tab.generate_button.click()
wait("the generated sequence", lambda: tab.compat is not None and tab.run_anyway.isHidden(), 120)
check("a new check clears 'Run anyway'", not tab.run_anyway.isChecked())
tab.run_button.click()
wait("the run", lambda: tab.result is not None and tab._run_worker is None and tab.result.gene.name == "SAMSN1", 600)
tab.breakpoints_button.click()
wait("the breakpoint sensitivity", lambda: tab._run_worker is None and tab.breakpoints_label.text(), 600)
check("after Run, 'Breakpoint sensitivity' moves both breakpoints by 10 and 50 bases either way and says whether the "
      "top transcript changes", "of 8 shifted breakpoints" in tab.breakpoints_label.text()
      and len(tab.result.breakpoints) == 1, tab.breakpoints_label.text())
tab.input_tabs.setCurrentIndex(0)
tab.window_edit.setText("")
tab.gene_edit.setText("USP25")                  # not computed yet in this test: its live computation can be cancelled
tab.gene_edit.editingFinished.emit()
wait("USP25", lambda: tab.gene is not None and tab.gene.name == "USP25", 60)
u_gene = tab.gene
tab.sequence_edit.setPlainText(Fasta(str(helpers.HG38))["chr" + u_gene.chrom][u_gene.mask[0] + 2000 - 1:
                                                                             u_gene.mask[0] + 30000].seq)
tab.check_button.click()
wait("USP25's check", lambda: tab.compat is not None, 120)
check("check notices (item 5): a paste of part of the gene is a caution -- a '?' on the right with the text, no '!', "
      "Run allowed", not tab.caution_button.isHidden() and "doesn't cover the whole gene" in tab.caution_button.toolTip()
      and tab.critical_button.isHidden() and tab.run_button.isEnabled() and "caution" in tab.check_status.text()
      and "doesn't cover" not in tab.report_label.text(), tab.check_status.text())
tab.run_button.click()
tab.cancel_button.click()
wait("the cancel", lambda: tab._run_worker is None, 120)
check("Cancel stops a run (here while the reference is computed live)", tab.progress_label.text() == "Cancelled.",
      tab.progress_label.text())
check("help texts: what the tab does, RNA confirmation, relative scores as rankings, typical use (from the same "
      "texts as the introduction)", all(k in cdna_help.TAB_HELP_TEXT + results_help
                                        for k in ("RNA (cDNA)", "a ranking", "typical use", "ATG")))

# --- the RefSeq tables missing: a clear error, never an empty gene list (spec 14 A3) ---------------
import spliceai_gui.cdna_tab as cdna_tab_module  # noqa: E402

real_refseq_dir, empty = cdna_tab_module.refseq_dir, tempfile.mkdtemp()
cdna_tab_module.refseq_dir = lambda: empty
tab.tables.clear()
tab.build_combo.setCurrentText("hg19")
on_build = tab.gene_label.text()
tab.gene_edit.setText("BRCA2")
tab.gene_edit.editingFinished.emit()
app.processEvents()
check("RefSeq table missing: choosing the build and entering a gene say which table is missing and what to do, "
      "in red; no gene, nothing to check", "is missing from" in on_build and "is missing from" in tab.gene_label.text()
      and "fetch-thirdparty.ps1" in tab.gene_label.text() and tab.gene is None and tab.gene_choice.count() == 0
      and not tab.check_button.isEnabled() and tab.gene_label.styleSheet() == theme.STATUS_CRIT, tab.gene_label.text())
cdna_tab_module.refseq_dir = real_refseq_dir
w.close()
helpers.finish()
