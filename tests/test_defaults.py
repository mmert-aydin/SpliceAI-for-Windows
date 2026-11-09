"""Spec 17 item 2: on a fresh install SnpEff and MANE Select are on in the VCF
tab; settings already saved stay as they are; a box turned off only because
SnpEff isn't there isn't saved as the user's choice.

    .venv\\Scripts\\python tests\\test_defaults.py
"""
import json
import os
import sys

import helpers
from helpers import check, skip

helpers.watchdog(300)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

from PySide6.QtWidgets import QApplication  # noqa: E402

from spliceai_gui import config  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402

check("config.DEFAULTS: SnpEff and MANE Select on", config.DEFAULTS["use_snpeff"] and config.DEFAULTS["use_mane"])
from spliceai_gui.cdna_help import PARAMETERS  # noqa: E402
from spliceai_pipeline import cdna_analysis, splice_paths  # noqa: E402

help_exon = next(x for x in PARAMETERS if x[0] == "Minimum internal exon")
check("the default minimum internal exon is 3 nt everywhere, and its help names micro-exons and incorrect extra exons",
      config.DEFAULTS["cdna_min_exon"] == splice_paths.MIN_EXON == cdna_analysis.Settings().min_exon == 3
      and help_exon[1] == "3 nt" and "icro-exons" in help_exon[2] and "incorrect ones" in help_exon[2])


def window(saved=None):
    if config.CONFIG_PATH.exists():
        config.CONFIG_PATH.unlink()
    if saved is not None:
        os.makedirs(config.CONFIG_PATH.parent, exist_ok=True)
        config.CONFIG_PATH.write_text(json.dumps(saved), encoding="utf-8")
    return MainWindow()


w = window()
if helpers.SNPEFF_DIR is None:
    skip("a fresh install annotates", "no installer/thirdparty")
else:
    check("a fresh install: 'Annotate with SnpEff' and 'MANE Select' ticked, from the default folders",
          w.use_snpeff_checkbox.isChecked() and w.use_mane_checkbox.isChecked() and w.use_mane_checkbox.isEnabled())
w._save_settings()
saved = json.loads(config.CONFIG_PATH.read_text(encoding="utf-8"))
check("...and saved as on", saved["use_snpeff"] and saved["use_mane"])

w = window({"use_snpeff": False, "use_mane": False, "snpeff_default_applied": True})
check("SnpEff is always on: no box to untick, and an 'off' saved by an older version is ignored",
      not w.use_snpeff_checkbox.isVisibleTo(w)
      and (w.use_snpeff_checkbox.isChecked() or helpers.SNPEFF_DIR is None))
w._save_settings()
check("...and it is saved as on", json.loads(config.CONFIG_PATH.read_text(encoding="utf-8"))["use_snpeff"])

w = window({"score_genes": ["ANKRD26"]})
check("Genes to score at a new start: nothing chosen (Run waits), the saved list kept to be offered",
      not w._score_genes_chosen and w._score_genes == ["ANKRD26"] and not w.run_button.isEnabled()
      and "choose All genes or a gene list" in w.run_button.toolTip())

w = window({"cdna_min_exon": 37})
check("a minimum exon saved earlier stays as it is", w.cdna_tab.min_exon.value() == 37)

w = window({"use_snpeff": True, "snpeff_dir": str(home / "no_snpeff_here"), "snpeff_default_applied": True})
check("SnpEff not found: annotation is off (and said so), so a run isn't blocked",
      not w.use_snpeff_checkbox.isChecked())
w._save_settings()
check("...but that isn't saved as a choice: on again once SnpEff is there",
      json.loads(config.CONFIG_PATH.read_text(encoding="utf-8"))["use_snpeff"])
helpers.finish()
