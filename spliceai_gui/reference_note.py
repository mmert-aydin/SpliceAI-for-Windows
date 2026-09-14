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

from . import config, reference_locator

REFERENCE_NOTE_HTML = (
    "<p><b>SpliceAI Variant Scoring is installed and ready.</b></p>"
    "<p>To score variants it needs a <b>reference genome FASTA</b> (hg19.fa or hg38.fa, about "
    "3 GB each). They aren't inside Setup because of their size &mdash; this is not a missing or "
    "broken installation. The program looks for them every time it opens, so this only appears "
    "when it couldn't find them.</p>"
    "<ul>"
    "<li><b>If you have the USB drive Setup came from</b>, the quickest fix is to run Setup again "
    "from it and leave <b>Copy the reference genomes to this PC</b> ticked. After that the "
    "program finds them by itself and never asks again.</li>"
    "<li><b>If you already have the files somewhere</b>, put the <b>.fa</b> and its <b>.fa.fai</b> "
    "together in a folder named <b>SpliceAI_reference_data</b> in your user folder, and they will "
    "be picked up next time you start. Or choose them by hand: <b>Advanced settings...</b> at the "
    "bottom of <b>Settings</b>, then <b>Browse...</b> on the hg19 or hg38 row. You can also drag "
    "a .fa file straight onto this window.</li>"
    "<li><b>If you have neither</b>, <b>Download...</b> on a row fetches that build from UCSC "
    "(about 1 GB). On a hospital or company network this can be refused with a certificate "
    "message; copying the files from the drive avoids that.</li>"
    "<li>Illumina's precomputed score files are optional. If you don't have them, choose "
    "<b>I don't have precomputed data</b>; every variant is then scored directly by SpliceAI.</li>"
    "</ul>"
)


def should_show(settings, paths=None):
    """True when the note is worth showing. `paths` is the reference FASTA of
    each build as the window actually has it -- which includes the ones found
    automatically (reference_locator), not only the ones saved in settings, so
    a machine where Setup copied the files never sees this note at all."""
    if settings.get("hide_reference_note"):
        return False
    if paths is None:
        paths, _detected = reference_locator.resolve(settings)
    return not any(path and os.path.isfile(path) for path in paths.values())


def _remember_hidden():
    settings = config.load()
    settings["hide_reference_note"] = True
    config.save(settings)


def show_reference_note_if_needed(parent, paths=None):
    """Shows the note modelessly over `parent` if needed; returns the box (or
    None) -- parented to `parent`, so it lives as long as the main window."""
    if not should_show(config.load(), paths):
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
