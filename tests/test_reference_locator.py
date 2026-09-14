"""Finding the reference genomes without asking the user, and the folded-away
reference section that replaced the two always-visible FASTA rows.

Needs no reference data of its own: the locator only looks at file names and at
whether a .fa.fai sits next to the file, so stand-in files are enough.

    .venv\\Scripts\\python tests\\test_reference_locator.py
"""
import os
import sys
import tempfile

import helpers
from helpers import check

helpers.watchdog(120)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
# The override is the first folder searched; unset here so the order below is
# the one a normal machine has.
os.environ.pop("SPLICEAI_REF_DIR", None)

from PySide6.QtWidgets import QApplication  # noqa: E402

from spliceai_gui import config, reference_locator, reference_note  # noqa: E402

tmp = tempfile.mkdtemp(prefix="spliceai_test_ref_")
DRIVE = os.path.join(tmp, "usb", "SpliceAI-USB-1.0.0", "reference-data")
LOCAL = os.path.join(tmp, "SpliceAI_reference_data")
EMPTY = os.path.join(tmp, "empty")
os.makedirs(EMPTY)


def fake_reference(directory, build, name=None, with_index=True):
    """A stand-in FASTA, optionally with its .fai index."""
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, name or f"{build}.fa")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">chr1 {build}\nACGT\n")
    if with_index:
        with open(path + ".fai", "w", encoding="utf-8") as fh:
            fh.write("chr1\t4\t12\t4\t5\n")
    return path


# --- what counts as usable -------------------------------------------------
hg19 = fake_reference(DRIVE, "hg19")
check("a FASTA with its .fai next to it is usable", reference_locator.is_usable(hg19))
no_index = fake_reference(DRIVE, "hg38", with_index=False)
check("without the .fa.fai it isn't -- building one takes minutes, and fails "
      "outright on a read-only drive", not reference_locator.is_usable(no_index))
check("so only hg19 is offered", reference_locator.find_references([DRIVE]) == {"hg19": hg19})

hg38 = fake_reference(DRIVE, "hg38")  # now with its index
check("both builds are found, each for its own row",
      reference_locator.find_references([DRIVE]) == {"hg19": hg19, "hg38": hg38})

# --- names other than hg19.fa ----------------------------------------------
other = os.path.join(tmp, "other")
named = fake_reference(other, "hg38", name="GRCh38_no_alt_analysis_set.fna")
check("a GRCh38-named .fna is recognised as hg38",
      reference_locator.find_references([other]) == {"hg38": named})
fake_reference(other, "hg38")
check("...but an exact hg38.fa in the same folder wins",
      reference_locator.find_references([other])["hg38"] == os.path.join(other, "hg38.fa"))
unclear = os.path.join(tmp, "unclear")
fake_reference(unclear, "hg19", name="genome.fa")
check("a name that doesn't say which build is left alone -- the wrong build "
      "would give wrong scores", reference_locator.find_references([unclear]) == {})

# --- order ------------------------------------------------------------------
local_hg19 = fake_reference(LOCAL, "hg19")
check("a copy on the PC is preferred over the same file on a drive",
      reference_locator.find_references([LOCAL, DRIVE])["hg19"] == local_hg19)
check("the search stops at the first usable file per build",
      reference_locator.find_references([EMPTY, DRIVE]) == {"hg19": hg19, "hg38": hg38})
real_dirs = reference_locator.candidate_dirs()
check("~\\SpliceAI_reference_data -- where Setup copies the drive's files, and "
      "where \"Download reference...\" puts them -- is looked in first",
      real_dirs[0] == str(reference_locator.USER_REFERENCE_DIR), real_dirs[0])
check("no network drive is ever searched (listing a dead share can hang for a "
      "long time)", all(":" not in d[2:] for d in real_dirs) and
      not any(d.startswith("\\\\") for d in real_dirs))

# --- resolve(): saved settings win, gaps are filled -------------------------
settings = dict(config.DEFAULTS)
paths, detected = reference_locator.resolve(settings, [DRIVE])
check("nothing saved: both builds are filled in automatically",
      paths == {"hg19": hg19, "hg38": hg38})
check("...and both are reported as found", sorted(detected) == ["hg19", "hg38"])

settings[config.fasta_key("hg19")] = local_hg19
paths, detected = reference_locator.resolve(settings, [DRIVE])
check("a saved path that still works is kept", paths["hg19"] == local_hg19)
check("...and is not reported as found automatically", "hg19" not in detected)
check("the build that has no setting is still filled in", paths["hg38"] == hg38)

settings[config.fasta_key("hg19")] = os.path.join(tmp, "gone", "hg19.fa")
paths, _ = reference_locator.resolve(settings, [DRIVE])
check("a saved path whose file has gone (drive unplugged) is replaced", paths["hg19"] == hg19)

# --- the first-launch reference note ----------------------------------------
check("the note is not shown when the files were found",
      reference_note.should_show(dict(config.DEFAULTS), {"hg19": hg19, "hg38": ""}) is False)
check("the note is shown when there is nothing to find",
      reference_note.should_show(dict(config.DEFAULTS), {"hg19": "", "hg38": ""}) is True)

# --- the window --------------------------------------------------------------
app = QApplication.instance() or QApplication([])
from spliceai_gui.main_window import MainWindow  # noqa: E402

# Fixed search folders, so this checks the window and not whatever happens to
# be on the machine running it.
reference_locator.candidate_dirs = lambda: [DRIVE]

w = MainWindow()
check("the window opens with the hg19 it found by itself", w.fasta_edits["hg19"].text() == hg19)
check("...and the hg38", w.fasta_edits["hg38"].text() == hg38)
check("nothing is wrong, so there is nothing to sort out", w.setup_problems() == [],
      "; ".join(w.setup_problems()))
check("advanced settings start folded away", w.advanced_widget.isHidden())
summary = w.setup_status_label.text()
check("the status line says the genomes are ready", "reference genome hg19 and hg38" in summary, summary)
check("...and which folder they came from", DRIVE in summary, summary)
check("the button offers to open advanced settings",
      w.advanced_button.text() == "Advanced settings...", w.advanced_button.text())

w.advanced_button.setChecked(True)
app.processEvents()
check("pressing it shows them", not w.advanced_widget.isHidden())
check("...and it then offers to fold them away again",
      w.advanced_button.text() == "Hide advanced settings")

local_hg38 = fake_reference(LOCAL, "hg38")
w.fasta_edits["hg19"].setText(local_hg19)
w.fasta_edits["hg38"].setText(local_hg38)
app.processEvents()
check("after choosing files by hand, the status stops saying it found them",
      "found in" not in w.setup_status_label.text(), w.setup_status_label.text())

w._save_settings()
saved = config.load()
check("what is in the rows is saved, so it stays bound to the program",
      saved[config.fasta_key("hg38")] == local_hg38, saved[config.fasta_key("hg38")])
w.close()

# Nothing to find: the rows start open, so a new user isn't left facing a
# button they would have to know to press.
reference_locator.candidate_dirs = lambda: [EMPTY]
for build in ("hg19", "hg38"):
    settings = config.load()
    settings[config.fasta_key(build)] = ""
    config.save(settings)
w = MainWindow()
check("with no reference anywhere, advanced settings start open",
      not w.advanced_widget.isHidden())
check("...and the FASTA rows are inside them", not w.fasta_rows_widget.isHidden())
check("...and the status line says what is missing",
      "reference genome hasn't been found" in w.setup_status_label.text(),
      w.setup_status_label.text())
check("...and the first-launch note would be shown",
      reference_note.should_show(config.load(), w.reference_paths()) is True)
check("...and the panes are the three sections", w.panes.count() == 3)

# --- draggable panes ---------------------------------------------------------
from spliceai_gui.main_window import DEFAULT_PANE_SIZES  # noqa: E402

check("the VCF box starts small -- it only has to show which file was loaded",
      DEFAULT_PANE_SIZES[1] < DEFAULT_PANE_SIZES[0], str(DEFAULT_PANE_SIZES))
check("no pane can be collapsed to nothing by dragging", not w.panes.childrenCollapsible())

# A splitter only has real heights to divide once the window has a size, so
# this part needs a shown window (offscreen, at a fixed size, so the same
# numbers come back both times).
WINDOW = (1100, 950)
w.resize(*WINDOW)
w.show()
app.processEvents()
w.panes.setSizes([400, 300, 200])
app.processEvents()
dragged = w.panes.sizes()
check("dragging a divider gives one section more room than the default",
      dragged != DEFAULT_PANE_SIZES and all(n > 0 for n in dragged), str(dragged))
w._save_settings()
check("where the dividers were left is saved", config.load()["pane_sizes"] == dragged,
      str(config.load()["pane_sizes"]))

w.close()
w = MainWindow()
w.resize(*WINDOW)
w.show()
app.processEvents()
check("...and restored next time the window opens", w.panes.sizes() == dragged,
      str(w.panes.sizes()))

settings = config.load()
settings["pane_sizes"] = [0, 0]          # an older or hand-edited settings file
config.save(settings)
w.close()
w = MainWindow()
w.resize(*WINDOW)
w.show()
app.processEvents()
check("a pane_sizes that doesn't fit is ignored rather than hiding a section",
      all(n > 0 for n in w.panes.sizes()), str(w.panes.sizes()))
w.close()

helpers.finish()
