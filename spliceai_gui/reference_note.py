"""Non-blocking reminder about reference data, shown after the main window
opens while no usable reference FASTA is set (for either build).

The installer can't include the reference genome (several GB) or Illumina's
precomputed score files, so a fresh install starts with an empty "Reference
FASTA" field. Without this note that is easy to mistake for a broken or
incomplete install. It's a modeless message box (the main window stays
usable), it stops appearing once a FASTA is set, and "Don't show this again"
turns it off for good (config key "hide_reference_note").
"""
import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox, QMessageBox

from . import config

REFERENCE_NOTE_HTML = (
    "<p><b>SpliceAI Variant Scoring is installed and ready.</b></p>"
    "<p>To score variants it needs a <b>reference genome FASTA</b> (hg19.fa or hg38.fa, about "
    "3 GB). It isn't part of the installer because of its size &mdash; this is not a missing "
    "or broken installation.</p>"
    "<ul>"
    "<li>If you were given reference files (for example on the same USB drive as Setup), copy "
    "them to this computer &mdash; the <b>.fa</b> file together with its <b>.fa.fai</b> file "
    "&mdash; and select each .fa in its own row, <b>hg19 reference FASTA</b> or <b>hg38 "
    "reference FASTA</b> (<b>Browse...</b>). The program picks the right build from each VCF.</li>"
    "<li>Otherwise, <b>Download...</b> in a row fetches that build from UCSC (about 1 GB download).</li>"
    "<li>Illumina's precomputed score files are optional. If you don't have them, choose "
    "<b>I don't have precomputed data</b>; every variant is then scored directly by SpliceAI.</li>"
    "</ul>"
)


def should_show(settings):
    if settings.get("hide_reference_note"):
        return False
    paths = [settings.get(config.fasta_key(build)) or "" for build in ("hg19", "hg38")]
    return not any(path and os.path.isfile(path) for path in paths)


def _remember_hidden():
    settings = config.load()
    settings["hide_reference_note"] = True
    config.save(settings)


def show_reference_note_if_needed(parent):
    """Shows the note modelessly over `parent` if needed; returns the box (or
    None) -- parented to `parent`, so it lives as long as the main window."""
    if not should_show(config.load()):
        return None
    box = QMessageBox(parent)
    box.setWindowTitle("Reference data")
    box.setIcon(QMessageBox.Information)
    box.setTextFormat(Qt.RichText)
    box.setText(REFERENCE_NOTE_HTML)
    box.setStandardButtons(QMessageBox.Ok)
    dont_show = QCheckBox("Don't show this again")
    box.setCheckBox(dont_show)
    box.setWindowModality(Qt.NonModal)
    box.finished.connect(lambda _result: dont_show.isChecked() and _remember_hidden())
    box.show()
    return box
