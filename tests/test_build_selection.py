"""Automatic genome build selection and the per-build FASTA rows, in the real
main window (Qt offscreen, throw-away home folder).

    .venv\\Scripts\\python tests\\test_build_selection.py
"""
import json
import os
import sys
import time

import helpers
from helpers import check

helpers.watchdog(120)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup(need_spliceai=False)
HG19, HG38 = str(helpers.HG19), str(helpers.HG38)

# Settings as an older version left them: one FASTA field, pointing at hg38.fa.
os.makedirs(home / ".spliceai_gui")
# Written with a byte-order mark, as some Windows editors do: it must still be read.
with open(home / ".spliceai_gui" / "config.json", "w", encoding="utf-8-sig") as fh:
    json.dump({"fasta_path": HG38, "build": "hg38", "precomputed_mode": "none"}, fh)

from pyfaidx import Fasta  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from spliceai_gui import config, reference_note  # noqa: E402
from spliceai_gui.vcf_info import detect_vcf_build  # noqa: E402

assert str(config.CONFIG_PATH).startswith(str(home)), config.CONFIG_PATH

HEAD, COLS = "##fileformat=VCFv4.2\n", helpers.VCF_COLUMNS + "chr1\t100\t.\tA\tG\t.\t.\t.\n"
cases = {
    "##reference=GRCh37": "hg19",
    "##reference=file:///ref/human_g1k_v37.fasta": "hg19",
    "##reference=/data/b37/genome.fa": "hg19",
    "##reference=hs37d5": "hg19",
    "##assembly=GRCh38": "hg38",
    "##reference=Homo_sapiens_assembly38.fasta": None,   # doesn't name the build...
    "##contig=<ID=chr1,length=249250621>": "hg19",        # ...but chr1's length does
    "##contig=<ID=1,length=248956422,assembly=x>": "hg38",
    "##contig=<ID=chr2,length=242193529>": None,
    "": None,
}
got = {h: detect_vcf_build(HEAD + (h + "\n" if h else "") + COLS)[0] for h in cases}
check("header rules", got == cases, str({h: g for h, g in got.items() if g != cases[h]}))
both = HEAD + "##reference=Homo_sapiens_assembly38.fasta\n##contig=<ID=chr1,length=248956422>\n" + COLS
check("an unhelpful ##reference falls through to ##contig", detect_vcf_build(both)[0] == "hg38")
check("guess_build_from_name",
      [config.guess_build_from_name(n) for n in ("hg38.fa", "GCF_000001405.13_GRCh37_genomic.fna", "genome.fa")]
      == ["hg38", "hg19", None])

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402

w = MainWindow()
w.show()


def settle(seconds=0.6):
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents()
        time.sleep(0.02)


settle()
check("a settings file saved with a byte-order mark is read", w.fasta_edits["hg38"].text() == HG38)
check("an old single FASTA setting moves to the hg38 row",
      w.fasta_edits["hg38"].text() == HG38 and w.fasta_edits["hg19"].text() == "")
check("the hg38 row is marked as used", "(used for this run)" in w.fasta_row_labels["hg38"].text()
      and "(used" not in w.fasta_row_labels["hg19"].text())

w.vcf_text.setPlainText(HEAD + "##reference=GRCh37\n" + COLS)
settle()  # through the debounce timer, like a real paste
check("pasting an hg19 VCF selects hg19", w.build_combo.currentText() == "hg19")
check("the label says where it came from", w.build_source_label.text().startswith("✓ from the VCF"))
check("the hg19 row is now marked as used", "(used for this run)" in w.fasta_row_labels["hg19"].text())
check("Run asks for the hg19 FASTA", "set the hg19 reference FASTA" in w.run_button.toolTip(), w.run_button.toolTip())
w.fasta_edits["hg19"].setText(HG19)
check("the selected build's FASTA is the one used", w._selected_fasta() == HG19)

# An hg38 VCF known only by chr1's length (no ##reference line), loaded from a file.
hg38_vcf = helpers.write_vcf(helpers.snv_lines(Fasta(HG38, rebuild=False), "chr1", [1_000_000, 2_000_000]),
                             ["##contig=<ID=chr1,length=248956422>"])
with open(hg38_vcf, "rb") as fh:  # saved with a byte-order mark, as some editors do
    raw = fh.read()
with open(hg38_vcf, "wb") as fh:
    fh.write(b"\xef\xbb\xbf" + raw)
w._load_vcf_from_path(hg38_vcf)
settle()
check("a VCF file saved with a byte-order mark loads without it", w.vcf_text.toPlainText().startswith("##fileformat"))
check("loading an hg38 VCF (known by chr1's length) selects hg38", w.build_combo.currentText() == "hg38",
      w.build_source_label.text())
check("...and uses the hg38 FASTA", w._selected_fasta() == HG38)

w.build_combo.setCurrentText("hg19")  # changed by hand
settle()
check("a manual choice is kept and flagged", w.build_combo.currentText() == "hg19"
      and w.build_source_label.text().startswith("⚠ the VCF says hg38"), w.build_source_label.text())
w.vcf_text.insertPlainText("\n")  # editing the same VCF doesn't undo it
settle()
check("editing the same VCF doesn't undo a manual choice", w.build_combo.currentText() == "hg19")

w.vcf_text.setPlainText(HEAD + COLS)  # no build information at all
settle()
check("a VCF without build information leaves the selection, and says so",
      w.build_combo.currentText() == "hg19" and "not stated" in w.build_source_label.text())
w.vcf_text.clear()
settle()
check("empty input clears the label", w.build_source_label.text() == "")

w._save_settings()
saved = json.load(open(config.CONFIG_PATH))
check("settings are saved per build; the old key is emptied",
      saved.get("fasta_path_hg19") == HG19 and saved.get("fasta_path_hg38") == HG38 and saved.get("fasta_path") == "")
check("the reference note stays hidden once a FASTA is set", reference_note.should_show(config.load()) is False)
check("the reference note shows when neither FASTA is set",
      reference_note.should_show({"fasta_path_hg19": "", "fasta_path_hg38": ""}) is True)
w.close()

helpers.finish()
