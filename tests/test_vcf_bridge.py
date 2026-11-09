"""Spec 16 item 4: from the VCF tab to the Sequence-to-cDNA tab. The variant
as a "From coordinates" event, its genotype read from the VCF, the engine's
DS / DP against SpliceAI's own (get_delta_scores), and in the real window
(Qt offscreen, throw-away home, no database: genes computed live): variants
scored by the VCF pipeline sent with "Analyse transcript effect" -- the
confirmation's contents, the prefilled tab, the run, the header with both
scores and the genotype, a variant in two genes, other variants of the gene
applied too, and a queue saved as one results file each.

The three regression variants (ClinVar demo variants: build, position, REF and
ALT only) and ClinVar demo variants are used; what their transcripts show
is reported, not checked as truth.

    .venv\\Scripts\\python tests\\test_vcf_bridge.py
"""
import json
import os
import sys
import time
import zipfile

import helpers
from helpers import check

helpers.watchdog(3000)
home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

os.makedirs(home / ".spliceai_gui")
with open(home / ".spliceai_gui" / "config.json", "w") as fh:
    json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg38",
               "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True,
               "mane_dir": str(helpers.MANE_DIR), "cdna_show_intro": False}, fh)

from types import SimpleNamespace  # noqa: E402

from spliceai_pipeline import bundle, vcf_bridge as vb  # noqa: E402
from spliceai_pipeline.cli import run_pipeline_core  # noqa: E402
from spliceai_pipeline.refdb_validate import delta_scores  # noqa: E402
from spliceai_pipeline.score import LiveScorer  # noqa: E402

# --- without Qt ----------------------------------------------------------------------------------------
check("event lines: the VCF tab's variant as HGVS g. on the chromosome's accession -- SNV, deletion, insertion, "
      "delins, a longer REF/ALT trimmed to the change",
      vb.event_line("hg38", "chr19", 45770000, "A", "G") == "NC_000019.10:g.45770000A>G"
      and vb.event_line("hg19", "chr22", 1000000, "CTGGGG", "C") == "NC_000022.10:g.1000001_1000005del"
      and vb.event_line("hg19", "chr7", 2000000, "G", "GC") == "NC_000007.13:g.2000000_2000001insC"
      and vb.event_line("hg38", "1", 100, "AT", "G") == "NC_000001.11:g.100_101delinsG"
      and vb.event_line("hg38", "chrX", 100, "ACG", "ATG") == "NC_000023.11:g.101C>T"
      and vb.event_line("hg38", "chr17", 3000000, "AC", "A") == "NC_000017.11:g.3000001del")
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
       "19\t100\t.\tA\tG,T\t.\t.\t.\tGT:DP\t0/1:20\t1/1:3\n19\t200\t.\tA\tG\t.\t.\t.\tGT\t1|1\t0/0\n"
       "19\t300\t.\tA\tG,T\t.\t.\t.\tGT\t1/2\t0/0\nX\t400\t.\tC\tT\t.\t.\t.\tGT\t1\t0\n"
       "19\t500\t.\tCAA\tCA\t.\t.\t.\tGT\t0/1\t0/0\n")
gt = lambda *a, **k: vb.genotype(VCF, *a, **k)
check("genotype from the VCF (first sample): heterozygous -- the other allele expected to give the WT transcript --, "
      "homozygous, another ALT on the other allele, hemizygous, not found",
      gt("chr19", 100, "A", "G").heterozygous and "expected to give the WT transcript" in gt("chr19", 100, "A", "G").text()
      and "the first of 2 samples" in gt("chr19", 100, "A", "G").text()
      and gt("chr19", 200, "A", "G").meaning == "homozygous" and gt("19", 300, "A", "T").meaning.startswith(
          "heterozygous, another ALT") and not gt("19", 300, "A", "T").heterozygous
      and gt("X", 400, "C", "T").meaning == "hemizygous" and gt("19", 999, "A", "G").gt == "")
check("genotype: a record the VCF tab normalized (left-aligned, trimmed) is found through the same normalization",
      gt("19", 500, "CA", "C", normalize=lambda c, p, r, a: (p, "CA", "C") if (p, r, a) == (500, "CAA", "CA") else
         (p, r, a)).heterozygous)

# --- the engine against get_delta_scores ---------------------------------------------------------------
# ClinVar demo variants: CHEK2 c.444+1G>A and MYBPC3 c.1227-13G>A in hg19, BRCA2 c.8488-1G>A in hg38.
REGRESSION = [("hg19", "chr22", 29121230, "C", "T"), ("hg19", "chr11", 47364709, "C", "T"),
              ("hg38", "chr13", 32370955, "G", "A")]
DEMO = [("hg38", "chr11", 47343158, "C", "T"), ("hg38", "chr22", 28725242, "C", "T"),
        ("hg38", "chr5", 112827247, "G", "A")]       # the last one lies in two of SpliceAI's genes

from PySide6.QtCore import QPoint, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication, QDialog, QMessageBox  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui.main_window import MainWindow  # noqa: E402
from spliceai_gui.transcript_effect_dialog import TranscriptEffectDialog  # noqa: E402

shown = []
for name in ("warning", "information", "critical"):
    setattr(QMessageBox, name, staticmethod(lambda _p, title, text, *a, **k: shown.append((title, text)) or QMessageBox.Ok))
QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Yes)      # closing the window at the end
w = MainWindow()
w.show()
tab = w.cdna_tab
app.processEvents()


def wait(until, seconds=900):
    deadline = time.time() + seconds
    while time.time() < deadline:
        app.processEvents()
        if until():
            return True
        time.sleep(0.05)
    return False


def vcf_text(variants, extra=()):
    lines = ["##fileformat=VCFv4.2", f"##reference={variants[0][0]}",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"]
    for build, chrom, pos, ref, alt in list(variants) + list(extra):
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0/1")
    return "\n".join(lines) + "\n"


def score(variants, build, mode="masked", extra=()):
    """The VCF tab's run, as its worker does it: the rows go into the results table."""
    text = vcf_text(variants, extra)
    path = helpers.write_vcf([])
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    rows = run_pipeline_core(path, build, mode, str(helpers.HG19 if build == "hg19" else helpers.HG38),
                             on_progress=lambda *a, **k: None)
    w.vcf_text.setPlainText(text)
    w.build_combo.setCurrentText(build)
    w._last_run_settings = {"build": build, "mode": mode, "used_precomputed": False, "used_snpeff": False}
    w._on_worker_finished(rows)
    app.processEvents()
    return rows


def select(variant_ids):
    from PySide6.QtCore import QItemSelectionModel

    view = w.table_view
    view.clearSelection()
    for r in range(w.page_table_model.rowCount()):
        if w.page_table_model.row_at(r).variant_id in variant_ids:
            view.selectionModel().select(w.page_table_model.index(r, 0),
                                         QItemSelectionModel.Select | QItemSelectionModel.Rows)
    return w.selected_result_rows()


def send(rows, hook):
    seen = {}

    def answer(dialog):
        seen["dialog"] = dialog
        return hook(dialog)
    before = tab.result
    w.analyse_transcript_effect(rows, dialog_hook=answer)
    ran = wait(lambda: tab.result is not before and tab._run_worker is None and not tab._vcf_queue
               and tab._pending_vcf is None and tab._vcf_folder is None)     # a queue's last message said
    return seen.get("dialog"), ran


reports = []


def report(label, result):
    req, comparison = result.from_vcf
    sample = next(t for t in result.transcripts if t.kind == "Patient" and t.rank == 1)
    wt = next(t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1)
    reports.append(f"{label} ({req.build} {req.name}, {req.mode}): {comparison.lines()[0]} | {comparison.lines()[1]} | "
                   f"{'within tolerance' if comparison.agree else 'DIFFERS: ' + '; '.join(comparison.reasons)} | "
                   f"sample 1 vs reference: {sample.change_vs_ref} | vs WT: {sample.change_vs_wt} | "
                   f"WT 1 vs reference: {wt.change_vs_ref}")


# hg38: regression variant 3 and two ClinVar demo variants, masked as the VCF tab scores by default.
hg38 = [v for v in REGRESSION + DEMO if v[0] == "hg38"]
rows38 = score(hg38, "hg38")
annotator = LiveScorer(str(helpers.HG38), "hg38", "masked")._annotator
first = next(r for r in rows38 if (r.chrom, r.pos) == (REGRESSION[2][1], REGRESSION[2][2]))
picked = select({first.variant_id})


def confirm(dialog):
    confirm.text = " | ".join(lbl.text() for lbl in dialog.findChildren(type(dialog.problem)))
    return dialog.compute_button.isEnabled()


dialog, ran = send(picked, confirm)
line = vb.event_line("hg38", first.chrom, first.pos, first.ref, first.alt)
check("the confirmation: build, variant (and the event it is entered as), gene, canonical transcript, genotype from "
      "the VCF, estimated time, and that only this variant is applied; Compute enabled",
      isinstance(dialog, TranscriptEffectDialog) and "hg38" in confirm.text and line in confirm.text
      and first.gene in confirm.text and "NM_" in dialog.transcript_label.text() and "0/1 (heterozygous" in confirm.text
      and "about " in dialog.estimate_label.text() and "Only this variant is applied" in confirm.text
      and not dialog.also.isChecked() and dialog.others_list.isHidden(), confirm.text[:400])
check("on Compute: the Sequence to cDNA tab, build and gene filled in, the variant entered in 'From coordinates', "
      "the sequence generated, checked and run", ran and w.tabs.currentWidget() is tab
      and tab.build_combo.currentText() == "hg38" and tab.gene is not None and tab.gene.name == first.gene
      and tab.input_tabs.currentIndex() == 1 and tab.coords_edit.toPlainText() == line
      and getattr(tab.result, "from_vcf", None) is not None, f"{ran} {tab.coords_edit.toPlainText()!r} {shown[-2:]}")
req, comparison = tab.result.from_vcf
_, live = delta_scores(annotator, first.chrom, first.pos, first.ref, first.alt, masked=True)
unrounded = live[first.gene]["ds"]
check("the engine's DS / DP for the variant match the VCF tab's within the spec 5.6 tolerance, and SpliceAI's "
      "unrounded values within float16 rounding", comparison.agree
      and max(abs(a - b) for a, b in zip(comparison.engine[0], unrounded)) <= 2 ** -12 + 1e-5 + 0.001,
      f"{comparison.lines()} vs unrounded {[round(x, 6) for x in unrounded]}")
header = tab.vcf_label.text()
check("the result header: the variant, its genotype -- the other allele expected to give the WT transcript --, the "
      "VCF tab's DS_AG/DS_AL/DS_DG/DS_DL with DP next to the engine's, within tolerance (no ?)",
      tab.vcf_header.isVisibleTo(tab) and "From the VCF tab" in header and "expected to give the WT transcript" in header
      and "VCF tab (masked): DS_AG" in header and "cDNA engine: DS_AG" in header and "within the rounding tolerance" in header
      and not tab.vcf_why.isVisibleTo(tab), header[:500])
report("regression variant 3", tab.result)
saved = tab.save_results(os.path.join(str(home), "one"))
with zipfile.ZipFile(saved) as z:
    manifest = json.loads(z.read("manifest.json"))
check("the results file keeps where it came from: variant, genotype, both scores, the verdict",
      manifest.get("from_vcf", {}).get("variant") == req.name and manifest["from_vcf"]["within_tolerance"]
      and manifest["from_vcf"]["vcf_tab"].startswith("VCF tab: DS_AG"))

# The ClinVar demo variants, one after the other: a queue, each saved as its own results file.
demo_rows = select({r.variant_id for r in rows38 if (r.chrom, r.pos) in {(v[1], v[2]) for v in DEMO[:2]}})
folder = str(home / "queue")
os.makedirs(folder)


def queue(dialog):
    dialog.folder_edit.setText(folder)
    return dialog.compute_button.isEnabled()


dialog, ran = send(demo_rows, queue)
files = sorted(f for f in os.listdir(folder) if f.endswith(bundle.EXTENSION))
check("a queue of 2: each variant run on its own and saved as its own results file, the list said at the end",
      ran and len(demo_rows) == 2 and len(files) == 2 and any("2 variant(s), results in" in t for _, t in shown[-2:]),
      f"{files} {shown[-1:]}")
for f in files:
    loaded, man, _ = bundle.load(os.path.join(folder, f))
    entry = man.get("from_vcf", {})
    reports.append(f"ClinVar demo {entry.get('variant')} ({entry.get('gene')}, masked): {entry.get('vcf_tab')} | "
                   f"{entry.get('cdna_engine')} | {'within tolerance' if entry.get('within_tolerance') else entry.get('reasons')}"
                   f" | sample 1 vs reference: {next(t for t in loaded.transcripts if t.kind == 'Patient').change_vs_ref}"
                   f" | vs WT: {next(t for t in loaded.transcripts if t.kind == 'Patient').change_vs_wt}"
                   f" | WT 1 vs reference: {next(t for t in loaded.transcripts if t.kind == 'WT predicted').change_vs_ref}")
check("...and in each, the engine's scores agree with the VCF tab's", all(
    bundle.load(os.path.join(folder, f))[1]["from_vcf"]["within_tolerance"] for f in files))

# A variant in two of SpliceAI's genes: the user chooses.
two = [r for r in rows38 if (r.chrom, r.pos) == (DEMO[2][1], DEMO[2][2])]
dialog = w.analyse_transcript_effect([two[0]], dialog_hook=lambda d: False)
genes, missing = dialog.found[0]
check("a variant scored in two genes: both offered (or the one the cDNA table lacks named), Cancel runs nothing",
      len(two) == 2 and len(genes) + len(missing) == 2 and (dialog.gene_combo.isVisibleTo(dialog) == (len(genes) > 1))
      and tab._pending_vcf is None, f"{[g.name for g in genes]} missing {missing}")

# Other variants of the gene applied too (unchecked by default).
near = (first.chrom, first.pos + 400)
from pyfaidx import Fasta  # noqa: E402

base = Fasta(str(helpers.HG38))[near[0]][near[1] - 1].seq.upper()
other_variant = ("hg38", near[0], near[1], base, "A" if base != "A" else "G")
rows38 = score(hg38, "hg38", extra=[other_variant])
first = next(r for r in rows38 if (r.chrom, r.pos) == (REGRESSION[2][1], REGRESSION[2][2]))


def click_box(box):
    """A real left click on a checkbox's own square."""
    QTest.mouseClick(box, Qt.LeftButton, Qt.NoModifier, QPoint(8, box.height() // 2))
    app.processEvents()


def with_other(dialog):
    with_other.before = (dialog.also.isChecked(), dialog.others_list.count(), dialog.others_list.isHidden())
    dialog.show()
    app.processEvents()
    click_box(dialog.also)
    with_other.shown = dialog.also.isChecked() and not dialog.others_list.isHidden()
    for k in range(dialog.others_list.count()):     # a click on the row's text, not only on its small box
        rect = dialog.others_list.visualItemRect(dialog.others_list.item(k))
        QTest.mouseClick(dialog.others_list.viewport(), Qt.LeftButton, Qt.NoModifier, rect.center())
    with_other.note = dialog.only_label.text()
    return True


dialog, ran = send(select({first.variant_id}), with_other)
req, comparison = tab.result.from_vcf
check("'Also apply other variants in this gene': unchecked by default, the VCF's other variants of the gene listed; "
      "ticked, they are entered together, and the header says so (phase unknown, one haplotype)",
      with_other.before == (False, 1, True) and with_other.shown and ran and len(req.others) == 1
      and len(tab.coords_edit.toPlainText().splitlines()) == 2
      and vb.event_line("hg38", req.others[0].chrom, req.others[0].pos, req.others[0].ref, req.others[0].alt)
      in tab.coords_edit.toPlainText() and "Applied together with 1 other variant" in tab.vcf_label.text()
      and "phase is unknown" in tab.vcf_label.text() and "ticked below are applied together" in with_other.note,
      tab.vcf_label.text()[:300])
check("...and they are in the sample: both variants are applied events",
      len(tab.result.edits) == 2 and tab.result.edits[1].start - tab.result.edits[0].start >= 300,
      str([(e.start, e.kind) for e in tab.result.edits]))


def alone(dialog):
    """A gene with no other scored variant: the box still answers, and says so."""
    dialog.show()
    app.processEvents()
    click_box(dialog.also)
    alone.state = (dialog.also.isChecked(), dialog.others_list.isHidden(), dialog.phase_label.isHidden(),
                   dialog.phase_label.text(), dialog.others())
    return False


w.analyse_transcript_effect([next(r for r in rows38 if (r.chrom, r.pos) == (DEMO[1][1], DEMO[1][2]))], dialog_hook=alone)
check("a gene with no other variant: the box still ticks, and says that there is none",
      alone.state[:3] == (True, True, False) and "no other scored variant in" in alone.state[3] and alone.state[4] == [],
      str(alone.state))
renamed = rows38[:]
for r in renamed:
    if r.chrom == near[0] and r.pos == near[1]:
        r.gene = "AN_OLDER_SYMBOL"
w._all_rows = renamed
dialog = w.analyse_transcript_effect([first], dialog_hook=lambda d: False)
check("another variant in the gene under another symbol is listed too (found by position)",
      dialog.others_list.count() == 1 and dialog.others_list.item(0).text().startswith(f"{near[0]}:{near[1]:,}"),
      str(dialog.others_list.count()))
genes_none = TranscriptEffectDialog(w, "hg38", [first], lambda r: ([], ["NOT_IN_TABLE"]), lambda r, g: [], lambda r: req.genotype,
                                    lambda g, n: (1, "x"), "masked")
genes_none.also.click()
check("a variant whose gene isn't in the table: the dialog says so and the box doesn't fail",
      not genes_none.compute_button.isEnabled() and genes_none.others() == [])
if not comparison.agree:
    check("...a difference that the other variant explains is said so in the '?'",
          tab.vcf_why.isVisibleTo(tab) and "other variants were applied" in tab.vcf_why.toolTip())

# hg19: regression variants 1 and 2, one at a time.
rows19 = score([v for v in REGRESSION if v[0] == "hg19"], "hg19")
annotator = LiveScorer(str(helpers.HG19), "hg19", "masked")._annotator
for k, (build, chrom, pos, ref, alt) in enumerate(REGRESSION[:2], 1):
    row = next(r for r in rows19 if (r.chrom, r.pos) == (chrom, pos))
    dialog, ran = send(select({row.variant_id}), lambda d: True)
    req, comparison = tab.result.from_vcf
    check(f"regression variant {k} (hg19): prefilled, run, the engine's DS / DP within tolerance of the VCF tab's",
          ran and tab.build_combo.currentText() == "hg19" and tab.coords_edit.toPlainText() == vb.event_line(
              build, chrom, pos, ref, alt) and comparison.agree, f"{comparison.lines()} {comparison.problems}")
    report(f"regression variant {k}", tab.result)

# The context menu reaches the same dialog, for the row clicked.
import spliceai_gui.main_window as main_window  # noqa: E402
from PySide6.QtWidgets import QMenu  # noqa: E402


class ChoosingMenu(QMenu):
    """The menu, with its item chosen instead of shown (PySide6 doesn't let QMenu.exec itself be replaced)."""
    def exec(self, *_):
        return next(a for a in self.actions() if a.text().startswith("Analyse transcript effect"))


opened = []
main_window.QMenu = ChoosingMenu
TranscriptEffectDialog.exec = lambda self: opened.append(self) or QDialog.Rejected
w.table_view.clearSelection()
w._on_results_context_menu(w.table_view.visualRect(w.page_table_model.index(0, 0)).center())
check("the results table's right-click menu has 'Analyse transcript effect...' for the row clicked",
      len(opened) == 1 and opened[0].rows[0].variant_id == w.page_table_model.row_at(0).variant_id)
w.tabs.setCurrentIndex(0)
check("'Analyse transcript effect...' is a button above the results table", w.transcript_effect_button.isVisibleTo(w)
      and w.transcript_effect_button.text() == "Analyse transcript effect...")
w.table_view.clearSelection()
shown.clear()
check("with nothing selected it says to select a variant; more than 10 are refused",
      w.analyse_transcript_effect() is None and "Select a variant" in shown[-1][1]
      and w.analyse_transcript_effect([rows19[0]] * 11) is None and "at most 10" in shown[-1][1])

print("\nWhat the transcripts show (reported, not checked as truth):")
for line in reports:
    print("  " + line)
w.close()
helpers.finish()
