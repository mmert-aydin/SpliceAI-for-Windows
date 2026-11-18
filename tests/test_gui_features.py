"""Run controls and the other window features, in the real main window (Qt
offscreen, throw-away home folder) with real scoring runs: Pause / Resume /
End, End during SnpEff, closing during a run, the CPU display, the MANE
folder, the SpliceAI check, and scrolling.

    .venv\\Scripts\\python tests\\test_gui_features.py
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import time

import helpers
from helpers import check, skip

helpers.watchdog(900)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg19",
               "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True}, fh)

from pyfaidx import Fasta  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox, QScrollArea  # noqa: E402

from spliceai_gui import config, spliceai_setup  # noqa: E402
from spliceai_gui import theme  # noqa: E402
from spliceai_gui.theme import STYLESHEET as theme_sheet  # noqa: E402
from spliceai_gui.cpu_meter import CpuMeter  # noqa: E402
from spliceai_pipeline.control import RunCancelled, RunControl  # noqa: E402
from spliceai_pipeline.mane import is_summary_filename  # noqa: E402

# --- units ---------------------------------------------------------------------------
c = RunControl()
c.pause()
check("RunControl: pause / resume", c.paused and (c.resume() or not c.paused))
c.cancel()
try:
    c.check()
    check("RunControl: an ended run raises RunCancelled", False)
except RunCancelled:
    check("RunControl: an ended run raises RunCancelled", True)
meter = CpuMeter()
t_end = time.time() + 0.6
while time.time() < t_end:
    sum(i * i for i in range(10000))
app_pct, pc_pct = meter.sample()
check("CpuMeter reports this program's and the whole PC's CPU",
      app_pct is not None and app_pct > 0 and (pc_pct is not None or sys.platform != "win32"), f"{app_pct} / {pc_pct}")
check("is_summary_filename", is_summary_filename(r"x\MANE.GRCh38.v1.5.summary.txt.gz")
      and not is_summary_filename("MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"))

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402


def settle(seconds=0.5):
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents()
        time.sleep(0.02)


def wait_for(cond, timeout, step=0.1):
    end = time.time() + timeout
    while time.time() < end:
        app.processEvents()
        if cond():
            return True
        time.sleep(step)
    return cond()


def scored(w):
    m = re.search(r"live-scored (\d+)/", w.progress_label.text())
    return int(m.group(1)) if m else 0


def snpeff_java_running():
    out = subprocess.run(
        ["powershell", "-NoProfile", "-Command",
         "@(Get-CimInstance Win32_Process -Filter \"Name='java.exe'\" | "
         "Where-Object { $_.CommandLine -match 'spliceai_snpeff_' }).Count"],
        capture_output=True, text=True).stdout.strip()
    return int(out or 0)


w = MainWindow()
w.show()
w.resize(1000, 620)  # a small laptop screen
settle(1.0)
scroll = w.centralWidget()
check("the whole window scrolls", isinstance(scroll, QScrollArea))
check("a small window gets a vertical scrollbar", scroll.verticalScrollBar().maximum() > 0)

check("the SpliceAI pill shows it's found", w.spliceai_status_button.text() == "SpliceAI ready",
      w.spliceai_status_button.text())
check("...and says what clicking it does", "Click to see its version" in w.spliceai_status_button.toolTip())

# Setup puts snpeff, java and mane side by side. Java and MANE are found beside
# the chosen SnpEff folder, not only beside the app -- without that, running
# from anywhere but the install reports "no Java" with a good one right there,
# and SnpEff silently switches itself off.
from spliceai_pipeline.mane import sibling_mane_dir  # noqa: E402
from spliceai_pipeline.snpeff import find_java, sibling_java_dir  # noqa: E402

if helpers.SNPEFF_DIR is None:
    helpers.skip("Java and MANE beside SnpEff", "no installer\thirdparty")
else:
    snpeff_dir = str(helpers.SNPEFF_DIR)
    check("the java folder beside SnpEff is found", sibling_java_dir(snpeff_dir) is not None,
          str(sibling_java_dir(snpeff_dir)))
    try:
        check("...and a Java 21+ in it", bool(find_java(snpeff_dir)), find_java(snpeff_dir))
    except Exception as exc:
        check("...and a Java 21+ in it", False, str(exc))
    check("the mane folder beside SnpEff is found", sibling_mane_dir(snpeff_dir) is not None,
          str(sibling_mane_dir(snpeff_dir)))
    check("nothing is claimed for a folder with no siblings",
          sibling_java_dir(None) is None and sibling_mane_dir("") is None)
info = spliceai_setup.check_installation()
check("the SpliceAI check: version 1.3.1, all model and annotation files",
      info and info["version"] == "1.3.1" and not info["missing_models"] and not info["missing_annotations"])

# --- 1.2.0: scoring only chosen genes, and the Zygosity column -------------------
from PySide6.QtCore import Qt  # noqa: E402

from spliceai_gui.table_model import COLUMNS  # noqa: E402

zyg_col = [a for _, a, _ in COLUMNS].index("zygosity")
check("Zygosity is the column right after Gene", [a for _, a, _ in COLUMNS][:2] == ["gene", "zygosity"])
w._apply_column_order([a for _, a, _ in COLUMNS if a != "zygosity"])  # saved by 1.1.1
check("...also for a column order saved before it existed",
      w.table_view.horizontalHeader().visualIndex(zyg_col) == 1)
# Genes to score: in the Input VCF section, and a choice Run waits for.
from PySide6.QtCore import QPoint, QTimer  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QGroupBox  # noqa: E402


def group_title(widget):
    p = widget.parentWidget()
    while p is not None and not isinstance(p, QGroupBox):
        p = p.parentWidget()
    return p.title() if p is not None else None


def click(button):
    QTest.mouseClick(button, Qt.LeftButton, Qt.NoModifier, QPoint(button.width() // 2, button.height() // 2))
    app.processEvents()


def in_gene_dialog(act, seen):
    """Runs act(dialog) on the modal Genes-to-score dialog once it is open (Qt's exec() blocks the test)."""
    def step():
        dialog = QApplication.activeModalWidget()
        if dialog is None or not hasattr(dialog, "paste_box"):
            QTimer.singleShot(100, step)
            return
        act(dialog, seen)
    QTimer.singleShot(100, step)


check("Genes to score sits in the Input VCF section, both choices next to each other, Run right of them",
      group_title(w.all_genes_button) == group_title(w.choose_score_genes_button) == group_title(w.run_button)
      == "Input VCF" and w.all_genes_button.x() < w.choose_score_genes_button.x() < w.run_button.x())
other = helpers.snv_lines(Fasta(str(helpers.HG19)), "chr2", [100_000_000])
w.vcf_text.setPlainText(
    "##fileformat=VCFv4.2\n##reference=hg19\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "chr10\t27326999\t.\tT\tC\t.\t.\t.\tGT\t0/1\n" + other[0].rstrip("\n") + "\tGT\t1/1\n")
settle(0.5)
check("...nothing chosen yet: Run is off and says why, neither button is marked",
      not w._score_genes_chosen and not w.run_button.isEnabled()
      and "choose All genes or a gene list" in w.run_button.toolTip()
      and not w.all_genes_button.isChecked() and not w.choose_score_genes_button.isChecked()
      and w.score_genes_label.text() == "Choose one to enable Run.", w.run_button.toolTip())
check("...All genes in the deeper orange, the gene list in the lighter one",
      w.all_genes_button.objectName() == "allGenesButton" and w.choose_score_genes_button.objectName() == "geneListButton"
      and theme.ORANGE in theme_sheet and theme.ORANGE_SOFT in theme_sheet)
click(w.all_genes_button)
check("All genes clicked: marked with a tick, Run on",
      w._score_genes_chosen and w._score_genes == [] and w.all_genes_button.isChecked()
      and w.all_genes_button.text() == "✓ All genes" and not w.choose_score_genes_button.isChecked()
      and w.run_button.isEnabled(), w.run_button.toolTip())


def pick_list(dialog, seen):
    dialog.paste_box.setPlainText("ANKRD26, NOTAGENE1")
    click(dialog.check_button)
    seen["summary"] = dialog.summary_label.text()
    seen["apply"] = dialog.apply_button.isEnabled()
    click(dialog.apply_button)


seen = {}
in_gene_dialog(pick_list, seen)
click(w.choose_score_genes_button)
check("Choose gene list... (the real dialog): a pasted list is checked, the unknown symbol named, Apply enabled",
      seen.get("apply") and "1 recognized, 1 not recognized" in seen.get("summary", ""), str(seen))
check("...the list is the choice: marked with a tick, All genes no longer, shown, Run on",
      w._score_genes_chosen and w._score_genes == ["ANKRD26"] and w.choose_score_genes_button.isChecked()
      and w.choose_score_genes_button.text() == "✓ Choose gene list..." and not w.all_genes_button.isChecked()
      and "1 gene(s): ANKRD26" in w.score_genes_label.text() and w.run_button.isEnabled())


def cancel(dialog, seen):
    seen["offered"] = dialog.paste_box.toPlainText()
    click(dialog.cancel_button)


seen = {}
in_gene_dialog(cancel, seen)
click(w.choose_score_genes_button)
check("...opened again it offers the same list; Cancel leaves the choice as it was",
      seen.get("offered") == "ANKRD26" and w._score_genes == ["ANKRD26"] and w.choose_score_genes_button.isChecked()
      and w.run_button.isEnabled(), str(seen))
w._save_settings()
check("...and the list is remembered", config.load()["score_genes"] == ["ANKRD26"])
vcf_now = w.vcf_text.toPlainText()
click(w.clear_vcf_button)
check("Clear: the choice is asked again (the list kept, to be offered)",
      not w._score_genes_chosen and w._score_genes == ["ANKRD26"] and not w.choose_score_genes_button.isChecked())
w.vcf_text.setPlainText(vcf_now)
settle(0.5)
check("...a VCF alone doesn't turn Run on", not w.run_button.isEnabled())
w._set_score_genes(["ANKRD26"], chosen=True)
w._on_run_clicked()
check("a run with one chosen gene finishes", wait_for(lambda: w._worker is None, 300), w.progress_label.text())
check("...the worker got the genes", "ANKRD26" in str(w._last_run_settings.get("score_genes")))
model = w.table_model
genes = [model.row_at(i).gene for i in range(model.rowCount())]
check("...only that gene's variant is in the results", genes == ["ANKRD26"], str(genes))
if genes:
    cell = model.index(0, zyg_col)
    check("...its Zygosity cell: 'Het' with the drawn icon and a tooltip",
          cell.data(Qt.DisplayRole) == "Het" and cell.data(Qt.DecorationRole) is not None
          and "Heterozygous" in (cell.data(Qt.ToolTipRole) or ""))
w._set_score_genes([])
w.vcf_text.setPlainText("")
settle(0.3)

w.use_snpeff_checkbox.setChecked(True)
if helpers.MANE_DIR:
    w.mane_dir_edit.setText(str(helpers.MANE_DIR))
    settle()
    check("MANE: a folder with the summary file -> found", w.mane_status_label.text().startswith("✓ Found (MANE v"),
          w.mane_status_label.text())
else:
    skip("MANE: found", "no installer/thirdparty/mane")
w.mane_dir_edit.setText(tempfile.mkdtemp())
settle()
check("MANE: a folder without it -> not found", w.mane_status_label.text() == "Not found")
w.mane_dir_edit.setText(str(helpers.MANE_DIR or ""))
w._save_settings()
check("MANE: the folder is remembered", json.load(open(config.CONFIG_PATH)).get("mane_dir") == str(helpers.MANE_DIR or ""))
w.use_snpeff_checkbox.setChecked(False)

# A real run: 40 SNVs in ANKRD26 (hg19), generated from the reference.
w.vcf_text.setPlainText(helpers.vcf_text(helpers.ankrd26_snvs(Fasta(str(helpers.HG19), rebuild=False), 40),
                                         ["##reference=hg19"]))
settle(1.0)
check("run setup: hg19 selected from the VCF, Run enabled",
      w.build_combo.currentText() == "hg19" and w.run_button.isEnabled(), w.run_button.toolTip())

w._on_run_clicked()
check("Pause and End are enabled while running", w.pause_button.isEnabled() and w.end_button.isEnabled())
check("live scoring starts", wait_for(lambda: scored(w) >= 10, 180), w.progress_label.text())
w._on_pause_clicked()
check("Pause turns into Resume", w.pause_button.text() == "Resume")
settle(4)  # let the variant in progress finish
# A meter of our own: the window's is also read by its 1-second display timer,
# and a reading taken right after that one would come back empty.
paused_meter = CpuMeter()
settle(6)
cpu_while_paused, _ = paused_meter.sample()
check("while paused, the program uses (almost) no CPU -- nothing is scored",
      cpu_while_paused is not None and cpu_while_paused < 5, f"{cpu_while_paused}%")
check("the CPU display is shown during a run", w.cpu_label.isVisible() and w.cpu_label.text().startswith("CPU:"))
n_paused = scored(w)
w._on_pause_clicked()
check("Resume continues scoring", wait_for(lambda: scored(w) > n_paused, 90), f"{n_paused} -> {scored(w)}")

QMessageBox.question = lambda *a, **k: QMessageBox.Yes
t0 = time.time()
w._on_end_clicked()
ended = wait_for(lambda: w._worker is None, 60)
check("End stops the run, with no results", ended and w.progress_label.text() == "Run ended -- no results."
      and w.table_model.rowCount() == 0 and not w.pause_button.isEnabled() and w.run_button.isEnabled(),
      f"{time.time() - t0:.1f}s")
check("the CPU display hides after the run", not w.cpu_label.isVisible())

if helpers.SNPEFF_DIR is None:
    skip("End during SnpEff", "no installer/thirdparty/snpeff")
else:
    w.use_snpeff_checkbox.setChecked(True)
    w.snpeff_dir_edit.setText(str(helpers.SNPEFF_DIR))
    settle()
    w._on_run_clicked()
    in_snpeff = wait_for(lambda: "Running SnpEff (" in w.progress_label.text() and snpeff_java_running() > 0,
                         120, 0.5)
    check("the run reaches SnpEff (Java running)", in_snpeff, w.progress_label.text())
    t0 = time.time()
    w._on_end_clicked()
    ended = wait_for(lambda: w._worker is None, 30)
    check("End during SnpEff stops within seconds and kills Java",
          ended and time.time() - t0 < 10 and snpeff_java_running() == 0, f"{time.time() - t0:.1f}s")
    w.use_snpeff_checkbox.setChecked(False)

w._on_run_clicked()
check("(run started again)", wait_for(lambda: scored(w) >= 1, 180), w.progress_label.text())
QMessageBox.warning = lambda *a, **k: QMessageBox.No
w.close()
settle(1)
check("closing during a run -> 'No' keeps the window and the run",
      w.isVisible() and w._worker is not None and w._worker.isRunning())
worker = w._worker
QMessageBox.warning = lambda *a, **k: QMessageBox.Yes
t0 = time.time()
w.close()
check("closing during a run -> 'Yes' stops the run cleanly, then closes",
      not w.isVisible() and not worker.isRunning(), f"{time.time() - t0:.1f}s")

# Scrolling a page down past a table: the wheel turn over the table's sideways bar doesn't slide its columns.
# (Passing the turn on to the page is Qt's for real wheel events; a sent one stops where it is ignored.)
from PySide6.QtCore import QPoint, QPointF  # noqa: E402
from PySide6.QtGui import QWheelEvent  # noqa: E402
from PySide6.QtWidgets import QTableWidget  # noqa: E402

from spliceai_gui import theme  # noqa: E402

theme.apply_theme(app)
table = QTableWidget(2, 40)            # forty columns: a sideways bar
table.resize(500, 160)
table.show()
app.processEvents()
bar = table.horizontalScrollBar()


def turn(widget, dx, dy, modifiers=Qt.NoModifier):
    event = QWheelEvent(QPointF(5, 5), QPointF(widget.mapToGlobal(QPoint(5, 5))), QPoint(0, 0), QPoint(dx, dy),
                        Qt.NoButton, modifiers, Qt.NoScrollPhase, False)
    app.sendEvent(widget, event)
    app.processEvents()
    return event


plain = turn(bar, 0, -120)
check("a plain wheel turn over a table's sideways bar leaves the columns where they are (the page scrolls)",
      bar.value() == 0 and not plain.isAccepted(), str(bar.value()))
turn(bar, -120, 0)
check("a sideways swipe (touchpad, tilt wheel) still moves the bar", bar.value() > 0, str(bar.value()))
table.close()

helpers.finish()
