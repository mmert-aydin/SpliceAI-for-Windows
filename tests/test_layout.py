"""Spec 14 C4: the Sequence-to-cDNA tab fits the window at 1366 and 1920 px
screen widths with Windows display scaling of 100%, 125% and 150% (a
logical width of the screen width divided by the scaling): no sideways
scrolling, and every control and help icon of the tab inside the window.

Windows' own font metrics decide these widths (the offscreen platform draws
text as boxes about twice as wide), so the window is laid out by the real
Windows platform but never shown on screen (WA_DontShowOnScreen).

    .venv\\Scripts\\python tests\\test_layout.py
"""
import json
import os
import sys

import helpers
from helpers import check

helpers.watchdog(600)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"
helpers.setup()

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg19",
               "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True,
               "mane_dir": str(helpers.MANE_DIR)}, fh)

from PySide6.QtCore import QPoint, Qt  # noqa: E402
from PySide6.QtWidgets import QAbstractButton, QApplication, QScrollArea, QWidget  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402

if app.platformName() != "windows":
    print("SKIP the layout check needs the Windows platform's fonts")
    helpers.finish()
    sys.exit(0)

w = MainWindow()
w.setAttribute(Qt.WA_DontShowOnScreen)
w.show()
tab = w.cdna_tab
w.tabs.setCurrentWidget(tab)
tab._on_continue_without_database()
tab.advanced.setChecked(True)                   # the widest state of the parameters
scroll = next(s for s in w.findChildren(QScrollArea) if s.isAncestorOf(tab))
problems = []
for screen in (1366, 1920):
    for scale in (1.0, 1.25, 1.5):
        width = round(screen / scale)
        w.resize(width, 760)
        for _ in range(3):
            app.processEvents()
        viewport = scroll.viewport()
        cut = [b.text() or b.toolTip()[:30] for b in tab.findChildren(QAbstractButton)
               if b.window() is w and b.isVisibleTo(tab)
               and b.mapTo(viewport, QPoint(b.width(), 0)).x() > viewport.width()]
        sideways = scroll.horizontalScrollBar().maximum()
        if sideways or cut:
            problems.append(f"{screen} px at {scale:.0%} ({width} logical): scrolls {sideways} px sideways, "
                            f"outside: {cut}")
check("the tab fits 1366 and 1920 px screens at 100%, 125% and 150% scaling: no sideways scrolling, every button "
      "and help icon inside (Advanced open)", not problems, "; ".join(problems))
inputs = tab.build_combo.parentWidget()
build_row = tab.build_combo.mapTo(inputs, QPoint(0, 0)).y()
check("the genome build has its own row, above Gene and Coordinates",
      build_row < tab.gene_edit.mapTo(inputs, QPoint(0, 0)).y()
      and abs(tab.gene_edit.mapTo(inputs, QPoint(0, 0)).y() - tab.window_edit.mapTo(inputs, QPoint(0, 0)).y()) < 4)
w.resize(round(1366 / 1.5), 760)
for _ in range(3):
    app.processEvents()
text_px = tab.gene_edit.fontMetrics().horizontalAdvance(tab.gene_edit.placeholderText())
check("the Gene field's placeholder is not cut off, even at the narrowest width (spec 15 item 4)",
      text_px + 16 <= tab.gene_edit.width(), f"{text_px} px of text in a {tab.gene_edit.width()} px field")
print(f"(the tab's minimum width: {tab.minimumSizeHint().width()} px)")

# VCF tab: closing Advanced settings gives back the height opening them took.
w.tabs.setCurrentIndex(0)
w.resize(1400, 900)
for _ in range(3):
    app.processEvents()
before = w.panes.sizes()
for checked in (True, False):
    w.advanced_button.setChecked(checked)
    for _ in range(3):
        app.processEvents()
    if checked:
        opened = w.panes.sizes()
after = w.panes.sizes()
check("Advanced settings: open makes Settings taller, close gives the height back to the results",
      opened[2] > before[2] and after[2] <= before[2] + 2 and after[0] >= before[0] - 2, f"{before} {opened} {after}")
w.close()

# The welcome notice: switched to Türkçe while open (longer text), nothing is cut off.
from spliceai_gui import i18n  # noqa: E402
from spliceai_gui.first_launch_dialog import FirstLaunchDialog  # noqa: E402

i18n._language = "en"
d = FirstLaunchDialog()
d.show()
d.language_combo.setCurrentIndex(list(i18n.LANGUAGES).index("tr"))
for _ in range(5):
    app.processEvents()
need = d.body_label.heightForWidth(d.body_label.width())
check("welcome notice in Türkçe (switched while open): the whole text fits", d.body_label.height() >= need,
      f"{d.body_label.height()} px for {need} px of text")
d.close()
i18n._language = "en"
helpers.finish()
