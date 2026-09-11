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

check("the SpliceAI button shows it's found", w.spliceai_status_button.text().startswith("✓ SpliceAI found"))
info = spliceai_setup.check_installation()
check("the SpliceAI check: version 1.3.1, all model and annotation files",
      info and info["version"] == "1.3.1" and not info["missing_models"] and not info["missing_annotations"])

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

helpers.finish()
