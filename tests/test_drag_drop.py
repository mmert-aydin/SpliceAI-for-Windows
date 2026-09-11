"""Drag and drop from Explorer, in the real main window (Qt offscreen,
throw-away home folder): real QDragEnterEvent / QDropEvent events with file
URLs, as Explorer sends them.

    .venv\\Scripts\\python tests\\test_drag_drop.py
"""
import json
import os
import sys
import tempfile
import time

import helpers
from helpers import check

helpers.watchdog(120)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup(need_spliceai=False)

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"build": "hg38", "precomputed_mode": "none", "hide_reference_note": True}, fh)

from PySide6.QtCore import QMimeData, QPoint, QPointF, Qt, QUrl  # noqa: E402
from PySide6.QtGui import QDragEnterEvent, QDropEvent  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402

w = MainWindow()
w.show()


def settle(seconds=0.6):
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents()
        time.sleep(0.02)


def mime_for(path=None, text=None):
    mime = QMimeData()
    if path:
        mime.setUrls([QUrl.fromLocalFile(str(path))])
    if text:
        mime.setText(text)
    return mime


def drag(widget, mime):
    """Drag enter + drop on `widget`; returns whether it accepted the drag."""
    enter = QDragEnterEvent(QPoint(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier)
    QApplication.sendEvent(widget, enter)
    accepted = enter.isAccepted()
    if accepted:
        QApplication.sendEvent(widget, QDropEvent(QPointF(5, 5), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier))
    settle()
    return accepted


def file_with(name, text):
    path = os.path.join(tempfile.mkdtemp(), name)
    with open(path, "w") as fh:
        fh.write(text)
    return path


hg19_vcf = file_with("sample.vcf", helpers.vcf_text(["chr10\t27326999\t.\tT\tC\t.\t.\t.\n"], ["##reference=hg19"]))
hg38_bgz = file_with("other.vcf.bgz", helpers.vcf_text(["chr1\t100\t.\tA\tG\t.\t.\t.\n"], ["##reference=hg38"]))
settle()

check("a VCF dropped onto the input box loads", drag(w.vcf_text, mime_for(hg19_vcf))
      and w.vcf_text.toPlainText().startswith("##fileformat") and "sample.vcf" in w.load_status_label.text())
check("...and its build is selected (hg19)", w.build_combo.currentText() == "hg19")

check("a VCF dropped anywhere else on the window loads", drag(w, mime_for(hg38_bgz))
      and "##reference=hg38" in w.vcf_text.toPlainText() and "other.vcf.bgz" in w.load_status_label.text())
check("...a .vcf.bgz name counts as a VCF, and its build is selected (hg38)", w.build_combo.currentText() == "hg38")

# The box refuses other files; Qt then offers them to the window, which explains.
shown = []
QMessageBox.information = lambda *a, **k: shown.append(a[1])
before = w.vcf_text.toPlainText()
drag(w.vcf_text, mime_for(file_with("notes.txt", "x")))
check("a non-VCF file dropped onto the input box is explained, not pasted in as file:/// text",
      shown == ["Can't use this file"] and w.vcf_text.toPlainText() == before and "file:///" not in before)

check("a FASTA dropped onto the window goes to its build's row (by name)",
      drag(w, mime_for(helpers.HG19)) and w.fasta_edits["hg19"].text() == str(helpers.HG19))
check("a FASTA dropped onto a row goes into that row, as a path",
      drag(w.fasta_edits["hg38"], mime_for(helpers.HG38)) and w.fasta_edits["hg38"].text() == str(helpers.HG38))

asked = []
QMessageBox.question = lambda *a, **k: (asked.append(a[1]), QMessageBox.No)[1]
drag(w.fasta_edits["hg19"], mime_for(helpers.HG38))
check("a FASTA whose name says hg38, dropped onto the hg19 row, asks first (No keeps the old one)",
      asked == ["Different build?"] and w.fasta_edits["hg19"].text() == str(helpers.HG19))

if helpers.MANE_DIR:
    mane_file = next(p for p in helpers.MANE_DIR.iterdir() if p.name.startswith("MANE.GRCh38"))
    check("a MANE summary file dropped onto the window sets the MANE folder",
          drag(w, mime_for(mane_file)) and w.mane_dir_edit.text() == str(helpers.MANE_DIR))
folder = tempfile.mkdtemp()
check("a folder dropped onto the SnpEff folder field sets it",
      drag(w.snpeff_dir_edit, mime_for(folder)) and os.path.normcase(w.snpeff_dir_edit.text()) == os.path.normcase(folder))

shown.clear()
before = w.vcf_text.toPlainText()
drag(w, mime_for(file_with("report.pdf", "x")))
check("any other file dropped onto the window gets a short explanation, and changes nothing",
      shown == ["Can't use this file"] and w.vcf_text.toPlainText() == before)
w.close()

helpers.finish()
