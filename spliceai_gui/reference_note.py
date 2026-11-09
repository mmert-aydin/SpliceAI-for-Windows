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
from .i18n import t

REFERENCE_NOTE_HTML = t(
    "<p><b>SpliceAI Variant Scoring is installed and ready.</b></p><p>Scoring variants requires a <b>reference "
    "genome FASTA</b> (hg19.fa or hg38.fa, about 3 GB each). Because of their size, these files are not included "
    "in Setup; this does not indicate a missing or faulty installation. The program looks for them each time it "
    "starts, and this note appears only when they cannot be found.</p><ul><li><b>If the USB drive used for Setup "
    "is available</b>, the quickest solution is to run Setup again from it and leave <b>Copy the reference "
    "genomes to this PC</b> selected. The program then finds the files by itself.</li><li><b>If the files are "
    "already on this computer</b>, place the <b>.fa</b> file and its <b>.fa.fai</b> file together in a folder "
    "named <b>SpliceAI_reference_data</b> in the user folder; they will be found at the next start. They can "
    "also be selected manually: <b>Advanced settings...</b> at the bottom of <b>Settings</b>, then "
    "<b>Browse...</b> on the hg19 or hg38 row. A .fa file can also be dragged onto this "
    "window.</li><li><b>Otherwise</b>, <b>Download...</b> on a row downloads that build from UCSC (about 1 GB). "
    "On hospital or institutional networks, this download may be refused with a certificate message; copying the "
    "files from the USB drive avoids this.</li><li>Illumina's precomputed score files are optional. If they are "
    "not available, select <b>I don't have precomputed data</b>; every variant is then scored directly by "
    "SpliceAI.</li></ul>"
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
    box.setWindowTitle(t("Reference data"))
    box.setIcon(QMessageBox.Information)
    box.setTextFormat(Qt.RichText)
    box.setText(REFERENCE_NOTE_HTML)
    box.setStandardButtons(QMessageBox.Ok)
    dont_show = QCheckBox(t("Don't show this again"))
    box.setCheckBox(dont_show)
    box.setWindowModality(Qt.NonModal)
    box.finished.connect(lambda _result: dont_show.isChecked() and _remember_hidden())
    box.show()
    return box
