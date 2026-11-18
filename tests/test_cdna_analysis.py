"""Phase 7, the engine under the tab: "Check compatibility" and "Run" for a
synthetic patient (cdna_analysis.py), the results table and the warnings,
the genes without a usable reference transcript, the results bundle saved
and loaded again (bundle.py: identical tables, coordinates and scores), and
bedGraph export read back through tabix (bedgraph.py).

    .venv\\Scripts\\python tests\\test_cdna_analysis.py
"""
import os
import random
import tempfile
import time
import zipfile

import helpers
from helpers import check

helpers.watchdog(900)
helpers.setup()

import numpy as np  # noqa: E402

import sys  # noqa: E402

from spliceai_pipeline import align, bedgraph, bundle, cdna_analysis as ca, haplotype, refdb, transcripts  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.tabix import TabixFile  # noqa: E402

rng = random.Random(3)
REFSEQ, MANE = str(helpers.THIRDPARTY / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
tmp = tempfile.mkdtemp(prefix="spliceai_cdna_analysis_")

real_db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
reader = refdb.Reader(real_db) if os.path.exists(os.path.join(real_db, "manifest.json")) else None
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=reader, cache_dir=os.path.join(tmp, "cache"),
                               refseq_dir=REFSEQ, mane_dir=str(helpers.MANE_DIR))
fasta = scores.predictor().annotator.ref_fasta
table = ca.GeneTable("hg38", REFSEQ, MANE)

# --- where the RefSeq tables are found (spec 14 A3) ---------------------------------
import hashlib  # noqa: E402

PINNED = {"hg19": "daf41b6c1a7465f71e0de1e680b973249b9345bab6adad15b48acd620ea30062",
          "hg38": "67090326b606658bdabdfe5f8dcb97b7943ac802c7945399522be7174ec021ca"}
found_dir = transcripts.default_refseq_dir()
check("from source (.venv): the RefSeq tables are found in the repository's refseq folder, both builds, the pinned "
      "files (installer/fetch-thirdparty.ps1)", os.path.normcase(found_dir) == os.path.normcase(
          transcripts.DEFAULT_REFSEQ_DIR) and all(hashlib.sha256(open(transcripts.table_path(b), "rb").read())
                                                  .hexdigest() == PINNED[b] for b in PINNED), found_dir)
frozen = os.path.join(tmp, "frozen")                     # the PyInstaller onedir layout
for folder in ("_internal/refseq",):
    os.makedirs(os.path.join(frozen, folder))
    for b in ("hg19", "hg38"):
        open(os.path.join(frozen, folder, transcripts.TABLE_NAME.format(build=b)), "w").close()
saved = sys.executable, getattr(sys, "_MEIPASS", None), transcripts.DEFAULT_REFSEQ_DIR
sys.frozen, sys.executable, sys._MEIPASS = True, os.path.join(frozen, "SpliceAI-VariantScoring.exe"),     os.path.join(frozen, "_internal")
transcripts.DEFAULT_REFSEQ_DIR = os.path.join(frozen, "refseq")      # as when imported frozen
in_bundle = transcripts.default_refseq_dir()
os.makedirs(os.path.join(frozen, "refseq"))
for b in ("hg19", "hg38"):
    open(os.path.join(frozen, "refseq", transcripts.TABLE_NAME.format(build=b)), "w").close()
beside = transcripts.default_refseq_dir()
del sys.frozen
sys.executable, transcripts.DEFAULT_REFSEQ_DIR = saved[0], saved[2]
if saved[1] is None:
    del sys._MEIPASS
check("frozen: the tables bundled by the .spec are found inside it (_internal/refseq); a 'refseq' folder next to "
      "the program comes first", os.path.normcase(in_bundle) == os.path.normcase(os.path.join(frozen, "_internal", "refseq"))
      and os.path.normcase(beside) == os.path.normcase(os.path.join(frozen, "refseq")), f"{in_bundle} / {beside}")
spec = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "SpliceAI-VariantScoring.spec"),
            encoding="utf-8").read()
check("the PyInstaller .spec bundles both tables (the repository's refseq/) into 'refseq', and the validation sets",
      all(f'("refseq/{b}.ncbiRefSeqCurated.txt.gz", "refseq")' in spec
          and f'("validation-sets/validation-set-{b}.h5", "validation-sets")' in spec for b in ("hg19", "hg38")))
empty = tempfile.mkdtemp(dir=tmp)
try:
    ca.GeneTable("hg38", empty, MANE)
    message = ""
except FileNotFoundError as exc:
    message = str(exc)
check("a missing table is an error that says which file, where, and what to do -- never an empty gene list",
      "hg38.ncbiRefSeqCurated.txt.gz" in message and empty in message and "fetch-thirdparty.ps1" in message, message)

# --- finding genes -----------------------------------------------------------------
check("genes by SpliceAI's symbol and by RefSeq's (FAM132A is C1QTNF12)",
      [g.name for g in table.find("brca2")] == ["BRCA2"] and [g.name for g in table.find("C1QTNF12")] == ["FAM132A"])
check("an unknown symbol finds nothing", table.find("NOT_A_GENE") == [])
table19 = ca.GeneTable("hg19", REFSEQ, MANE)
kind, text = table19.reference_note(table19.find("BOP1")[0])
check("hg19 BOP1: no NM_ at the locus -> explained, no anchors", kind == "none" and "NM_015201" in text, text[:120])
kind, text = table19.reference_note(table19.find("RHBG")[0])
FALLBACK_NOTE = ("No MANE Select transcript is available for this gene on hg19. The reference cDNA and the anchors "
                 "(transcript start, end and ATG) are taken from NM_001412175.1.")
check("hg19 RHBG: no MANE Select transcript on hg19 -> the note names the fallback NM_ it is taken from",
      kind == "fallback" and text == FALLBACK_NOTE, text)
check("hg38 BRCA2 has a reference cDNA (no note)", table.reference_note(table.find("BRCA2")[0]) == (None, None))
check("exons named in the canonical numbering; a second exon within the same canonical one by its coordinates",
      ca.exon_labels([(100, 200), (300, 400), (500, 600), (700, 900), (1200, 1300), None],
                     [(100, 200), (300, 400), (500, 900)]) ==
      ["Exon 1", "Exon 2", "Exon 3", "its extra exon (700-900)", "its extra exon (1,200-1,300)",
       "its exon in inserted sequence"])

# --- a synthetic patient through compatibility and run ---------------------------------
gene = table.find("SAMSN1")[0]
record = fasta["chr" + gene.chrom]
fetch = lambda a, b: record[a - 1:b].seq
c = gene.canonical
exons = list(reversed(c.exons))         # minus strand: transcript order
a3, b3 = exons[2]
region = (gene.mask[0] - 5000, gene.mask[1] + 5000)
ref_region = fetch(*region).upper()
donor = a3 - 1                          # donor +1 of exon 3 (minus strand: below the exon)
edits = [Edit(donor, donor, "A" if ref_region[donor - region[0]] != "A" else "C", "SNV"),
         Edit(exons[5][1] + 900, exons[5][1] + 899, "".join(rng.choice("ACGT") for _ in range(70)), "insertion")]
patient_seq, _ = haplotype.build(ref_region, region[0], edits)
paste = ">patient\n" + "\n".join(patient_seq[i:i + 60] for i in range(4000, len(patient_seq) - 4000, 60))
compat = ca.check_compatibility(paste, gene, fetch, len(record))
kinds = sorted(e.kind for e in compat.alignment.events)
check("check: both edits found; the 70 bp insertion applied by default, the SNV not",
      kinds == ["SNV", "insertion"] and [e.apply for e in sorted(compat.alignment.events, key=lambda e: e.kind)]
      == [False, True], f"{kinds}")
notes = [text for _, text, _ in ca.check_notices(compat)]
check("check: the SNV left out next to the insertion is named in a notice",
      any("1 of 2 differences are not applied" in n for n in notes), str(notes))
lone = [align.Event("SNV", 100, "A", "C", 100, 100, "C", 1), align.Event("deletion", 200, "AT", "A", 201, 201, "", 1)]
align.set_default_apply(lone, 3)
alone = [e.apply for e in lone]
big = lone + [align.Event("deletion", 300, "ATTTT", "A", 301, 304, "", 4)]
align.set_default_apply(big, 3)
check("a lone SNV is applied by default, a 1 bp deletion is not; next to a 4 bp deletion neither is",
      alone == [True, False] and [e.apply for e in big] == [False, False, True], f"{alone} {[e.apply for e in big]}")
check("check: the report gives identity and the covered interval; the paste covers the whole gene, so no warning",
      any(line.startswith("Identity") for line in compat.report["lines"]) and not compat.report["problems"]
      and not compat.report["warnings"], str(compat.report)[:200])
short = ca.check_compatibility(paste[:30_000], gene, fetch, len(record))
check("check: a paste of part of the gene is warned about", any("doesn't cover the whole gene" in w
                                                                for w in short.report["warnings"]))
for e in compat.alignment.events:
    e.apply = True
check("check: the chosen events pass (no overlaps, inside the region)", not ca.chosen_check(compat)["problems"])
steps = []
result = ca.run_analysis("hg38", gene, compat, scores, fetch, table, on_progress=lambda f, t: steps.append((f, t)))
rows = ca.table_rows(result)
check("run: progress goes from start to done", steps[0][0] < 0.1 and steps[-1] == (1.0, "Done"))
check("run: the table has the reference, the wild type's 5 and the patient's 5 transcripts",
      [r["Type"] for r in rows] == ["Reference"] + ["WT predicted"] * 5 + ["Sample"] * 5
      and rows[0]["Change vs Reference"] == "No change", str([r["Type"] for r in rows]))
check("run: the wild type's top transcript is the reference's", rows[1]["Change vs Reference"] == "No change",
      rows[1]["Change vs Reference"])
p1 = rows[6]
check("run: the patient's top transcript skips exon 3 (in frame), said briefly with r. and p.",
      p1["Change vs Reference"].startswith("Exon 3 skipped") and "in frame" in p1["Change vs Reference"]
      and "r.(" in p1["Change vs Reference"] and "p.(" in p1["Change vs Reference"]
      and p1["Change vs WT predicted"].startswith("Exon 3 skipped"), p1["Change vs Reference"])
check("run: relative scores are shown as numbers, not percentages",
      all("%" not in r["Relative score"] for r in rows) and rows[6]["Relative score"].startswith("0."))
check("run: the patient's top transcript (its exon-3 donor lost) is not 'weak evidence': max delta >= 0.5",
      p1["Evidence vs WT"].startswith("max delta") and float(p1["Evidence vs WT"].split()[-1]) >= 0.5,
      p1["Evidence vs WT"])
check("run: how clearly rank 1 leads rank 2 is said", ca.lead_text(result).startswith("Sample rank 1 leads rank 2 by "),
      ca.lead_text(result))
sample1 = next(t for t in result.transcripts if t.kind == "Patient" and t.rank == 1)
rna = ca.rna_check(result, sample1)
exon3_len = exons[2][1] - exons[2][0] + 1
check("RNA check (spec 15 item 10): exon 3 skipped -> primers in exons 2 and 4; the predicted product is exon 3 "
      "shorter than the reference's, both measured from the cDNAs", rna and (rna["up"], rna["down"]) == (2, 4)
      and rna["reference"] - rna["predicted"] == exon3_len
      and rna["text"] == f"exon 2 -> exon 4: reference {rna['reference']:,} bp, predicted {rna['predicted']:,} bp"
      and p1["RNA check (RT-PCR)"] == rna["text"], str(rna))
wt_top = next(t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1)
check("RNA check: none for the WT prediction, nor for a sample candidate that is the WT's transcript",
      ca.rna_check(result, wt_top) is None
      and all(ca.rna_check(result, t) is None for t in result.transcripts if t.kind == "Patient" and t.max_delta is None))
local = ca.local_concordance(result)
check("local WT concordance (spec 14 C2): the SNV in intron 3 -> exons 2-5, the insertion in intron 5 -> exons 4-7; "
      "the WT matches the reference for all of them", sorted((d["where"], d["exons"], d["differing"]) for d in local)
      == [("intron 3", [2, 3, 4, 5], {}), ("intron 5", [4, 5, 6, 7], {})]
      and any(d["text"].endswith("matches the reference for exons 2-5 (4 of 4).") for d in local),
      [d["text"] for d in local])
import copy  # noqa: E402

fake = copy.copy(result)
wt1 = copy.deepcopy(next(t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1))
wt1.vs_ref = (wt1.vs_ref[0], [dict(r, status="end shifted +12") if r["ref"] == 4 else r for r in wt1.vs_ref[1]],
              *wt1.vs_ref[2:])
fake.transcripts = [wt1 if t.kind == "WT predicted" and t.rank == 1 else t for t in result.transcripts]
ex5, last = exons[4], exons[-1]
fake.edits = [Edit(ex5[0] + 10, ex5[0] + 30, "", "deletion"), Edit(last[0] - 50, last[0] - 50, "A", "SNV"),
              Edit(c.end + 200, c.end + 199, "ACGT", "insertion")]
fl = ca.local_concordance(fake)
check("local WT concordance: an edit in exon 5 -> exons 3-7; one after the transcript -> its last two exons; one "
      "before it -> its first two; a WT exon that differs is named with how",
      [(d["where"], d["exons"]) for d in fl] == [("exon 5", [3, 4, 5, 6, 7]), ("after the transcript", [7, 8]),
                                                 ("before the transcript", [1, 2])]
      and fl[0]["differing"] == {4: "end shifted +12"} and "exon 4 (end shifted +12)" in fl[0]["text"]
      and "(4 of 5 exons identical)" in fl[0]["text"], [d["text"] for d in fl])
# The DMPK-like case (spec 15 item 9): the WT prediction has an exon of its own in intron 1.
fake2 = copy.copy(result)
wt_x = copy.deepcopy(next(t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1))
mid1 = (exons[1][1] + exons[0][0]) // 2                 # minus strand: intron 1 lies between exon 2 and exon 1
wt_x.exons = list(wt_x.exons) + [(mid1, mid1 + 71)]
wt_x.ref_exons = list(wt_x.ref_exons) + [(mid1, mid1 + 71)]
wt_x.vs_ref = (wt_x.vs_ref[0], wt_x.vs_ref[1] + [{"ref": None, "other": len(wt_x.exons), "status": "extra (novel)"}],
               *wt_x.vs_ref[2:])
fake2.transcripts = [wt_x if t.kind == "WT predicted" and t.rank == 1 else t for t in result.transcripts]
fake2.edits = [Edit(mid1 - 500, mid1 - 451, "", "deletion")]
regions = ca.wt_regions(fake2)
l2 = ca.local_concordance(fake2)[0]
check("a WT-only exon in intron 1 marks intron 1 as discordant -- in the regions the schematic shades and in the "
      "local concordance of an event there", regions == {("intron", 1): "an exon of 72 nt in the WT prediction only"}
      and l2["where"] == "intron 1" and l2["introns"] == {1: "an exon of 72 nt in the WT prediction only"}
      and "intron 1 (an exon of 72 nt in the WT prediction only)" in l2["text"], f"{regions} | {l2['text']}")

# The Phase 8 control (decided after it): a random 65 bp insertion deep in SAMSN1's intron 1. SpliceAI raises a
# weak donor 32 nt away (DG 0.22), above exon 1's own (0.158): a model effect, labelled, ranking unchanged.
RANDOM_65 = "CAGGCGATGTGTCTACACCGAATGCTCCTTTTAAGAAAAGCTCACACGTAGGGGATCAACCGTTA"
weak_seq, _ = haplotype.build(ref_region, region[0], [Edit(14533713, 14533712, RANDOM_65, "insertion")])
weak_compat = ca.check_compatibility(">sample\n" + weak_seq, gene, fetch, len(record))
for e in weak_compat.alignment.events:
    e.apply = True
weak = ca.run_analysis("hg38", gene, weak_compat, scores, fetch, table)
samples = [t for t in weak.transcripts if t.kind == "Patient"]
labels = [r["Evidence vs WT"] for r in ca.table_rows(weak) if r["Type"] == "Sample"]
check("SAMSN1 random insertion: the top candidate is 'weak evidence (delta < 0.5)', max delta 0.22; rank 2 (the WT's "
      "transcript) is the same as the WT prediction", labels[0] == "weak evidence (delta < 0.5); max delta 0.22"
      and labels[1] == "same as the WT prediction", str(labels[:2]))
check("SAMSN1 random insertion: the label leaves the ranking as it was; the weak-evidence setting moves it",
      [t.relative for t in samples] == sorted((t.relative for t in samples), reverse=True)
      and "weak" not in ca.evidence(samples[0], ca.Settings(weak_delta=0.2)), ca.lead_text(weak))

# --- a single position in the Coordinates field (spec 14 C5) ----------------------------
h = 14533713 - region[0] + 1                    # the random insertion's place on that haplotype
piece = weak_seq[h - 2000:h + 2065]
near = ca.check_compatibility(piece, gene, fetch, len(record), position=14533713)
span = 2 * len(piece) + 10_000
check("a single position: the paste is looked for around it and found; the range found is given",
      near.found == (14531714, 14535713)
      and near.report["lines"][0] == f"Found within {span:,} bases of 14,533,713: chr21:{near.found[0]:,}-"
      f"{near.found[1]:,}." and [e.kind for e in near.alignment.events] == ["insertion"], near.report["lines"][0])
far = ca.check_compatibility(piece, gene, fetch, len(record), position=gene.mask[1] - 100)
check("a single position the paste isn't near: the whole gene is searched, and it is found there",
      far.found == near.found and far.report["lines"][0].startswith(f"Not found within {span:,} bases of ")
      and "; found in the whole gene: chr21:" in far.report["lines"][0], far.report["lines"][0])
check("the genes at a position, for a missing gene name: SAMSN1 (MANE Select) first at its insertion; none in a "
      "gene desert",
      [g.name for g in table.at("chr21", 14533713)][0] == "SAMSN1" and table.at("21", 1, 1000) == [],
      [g.name for g in table.at("chr21", 14533713)])

# --- the sequence-specificity control (spec 14 C3) -------------------------------------
from spliceai_pipeline import specificity  # noqa: E402

shuffles = []
for k in range(300):
    s = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 200)))
    shuffles.append((s, specificity.dinucleotide_shuffle(s, rng)))
check("dinucleotide shuffle: same length, first and last base and dinucleotide counts (300 random sequences), "
      "and not just the input back", all(len(a) == len(b) and a[0] == b[0] and a[-1] == b[-1]
                                         and specificity.counts(a) == specificity.counts(b) for a, b in shuffles)
      and sum(a != b for a, b in shuffles) > 280, "")
lo_i, hi_i = exons[3][1] + 1, exons[2][0] - 1      # intron 3 (minus strand: between exon 4 and exon 3 above it)
deletion = Edit((lo_i + hi_i) // 2, (lo_i + hi_i) // 2 + 79, "", "deletion")
stand_ins = [specificity.replacement(result, deletion, rng) for _ in range(50)]
check("deletion stand-ins: the same size, at random positions inside the same intron",
      all(d.end - d.start + 1 == 80 and lo_i <= d.start and d.end <= hi_i for d in stand_ins)
      and len({d.start for d in stand_ins}) > 40, f"{lo_i}..{hi_i}")
t0 = time.time()
controls = specificity.run(weak, scores, n=10)
ctl = controls[0]
check("SAMSN1 random insertion: observed top-1 changed with max delta 0.22 (SpliceAI's own DS), next to 10 random "
      "stand-ins of the same composition at the same place", len(controls) == 1 and ctl.observed_changed
      and abs(ctl.observed_delta - 0.22) < 0.01 and ctl.n == 10 and 0 <= ctl.changed <= 10
      and ctl.text().startswith("The insertion at 14,533,713 (65 bp): observed: top-1 changed, max delta 0.22; "
                                "random: "), f"{ctl.text()} in {time.time() - t0:.0f} s")
snv_only = specificity.run(result, scores, n=3)
ins_ctl = next(c for c in snv_only if c.kind == "insertion")
check("two events: the SNV has no random stand-in (said so, left out); the 70 bp insertion is measured without the "
      "SNV's effect -- the top-1 doesn't change because of it, and its delta stays small",
      sorted(c.n for c in snv_only) == [0, 3] and not ins_ctl.observed_changed and ins_ctl.observed_delta < 0.5
      and max(ins_ctl.deltas) < 0.5, [c.text() for c in snv_only])
weak.specificity = controls
saved_path = bundle.save(weak, os.path.join(tmp, "specificity"))
check("the control's result is saved in the bundle and comes back with it",
      bundle.load(saved_path)[0].specificity == controls, "")
notices = ca.result_notices(result)
check("notices (spec 15 item 5): what is always true is no notice -- the 70 bp insertion's 'outside typical use', "
      "a live computation; what is left is a caution or critical, each next to its rows",
      not any("typical use" in t or "computed live" in t for _, t, _ in notices)
      and all(level in (ca.CAUTION, ca.CRITICAL) and rows for level, _, rows in notices)
      and result.warnings == [t for _, t, _ in notices], str(notices)[:300])
import copy as _copy  # noqa: E402

from spliceai_pipeline import align  # noqa: E402

partial = ca.check_notices(short)
overlap = _copy.deepcopy(compat)
first = next(e for e in overlap.alignment.events if e.apply)
overlap.alignment.events.append(_copy.deepcopy(first))
at_atg = _copy.deepcopy(compat)
atg = at_atg.anchors["ATG"]
at_atg.alignment.events.append(align.Event("SNV", atg, "A", "C", atg, atg, "C", 1, apply=True))
check("notice classes: a paste of part of the gene is a caution; overlapping events are critical; a change "
      "touching the ATG is critical but can be run after confirming",
      [(lv, conf) for lv, t, conf in partial if "doesn't cover the whole gene" in t] == [(ca.CAUTION, False)]
      and any(lv == ca.CRITICAL and not conf and "overlap" in t for lv, t, conf in ca.check_notices(overlap))
      and any(lv == ca.CRITICAL and conf and "ATG" in t for lv, t, conf in ca.check_notices(at_atg)))

# --- the bundle -----------------------------------------------------------------------------
path = bundle.save(result, os.path.join(tmp, "result"), versions={"spliceai": "1.3.1"})
check("bundle: one file with the program's extension", path.endswith(bundle.EXTENSION) and os.path.isfile(path))
with zipfile.ZipFile(path) as z:
    names = set(z.namelist())
    summary = z.read("summary.tsv").decode()
expected = {"manifest.json", "summary.tsv", "concordance.tsv", "transcripts.tsv", "edits.tsv", "reference_scores.tsv",
            "sample_scores.tsv", "cdna.fasta", "protein.fasta", "alignments.txt", "analysis.json", "region.fa",
            "haplotype.fa", "wt_values.npy", "patient_values.npy"}
check("bundle: the human-readable tables, FASTA, alignments and the manifest are inside", expected <= names,
      str(sorted(expected - names)))
check("bundle: summary.tsv is the results table", summary.splitlines()[0].split("\t") == ca.COLUMNS
      and len(summary.splitlines()) == len(rows) + 1)
loaded, manifest, notes = bundle.load(path)
check("bundle: loaded again, the table is identical", ca.table_rows(loaded) == rows,
      str([(a, b) for a, b in zip(ca.table_rows(loaded), rows) if a != b][:1])[:300])
check("bundle: ...and so are coordinates, exons and per-base scores",
      [t.exons for t in loaded.transcripts] == [t.exons for t in result.transcripts]
      and [t.ref_exons for t in loaded.transcripts] == [t.ref_exons for t in result.transcripts]
      and np.array_equal(loaded.wt_values, result.wt_values) and np.array_equal(loaded.patient.values, result.patient.values)
      and loaded.hap == result.hap and loaded.cmap.segments == result.cmap.segments
      and [t.cdna.protein for t in loaded.transcripts] == [t.cdna.protein for t in result.transcripts])
check("bundle: the manifest records build, gene, transcript, parameters, versions and checksums",
      manifest["build"] == "hg38" and manifest["gene"] == "SAMSN1" and manifest["canonical"] == c.accession
      and manifest["parameters"]["threshold"] == 0.05 and manifest["spliceai"] == "1.3.1"
      and set(manifest["files"]) == names - {"manifest.json"} and not notes)
check("bundle: the manifest records each sample candidate's evidence and the lead of rank 1 over rank 2",
      manifest["evidence"]["Sample 1"]["label"] == rows[6]["Evidence vs WT"]
      and manifest["research_use"] == ca.RESEARCH_NOTE and manifest["strand"] == "-"
      and manifest["rna_check"]["Sample 1"]["up"] == 2
      and manifest["parameters"]["weak_delta"] == 0.5 and manifest["sample_lead"] == ca.lead_text(result)
      and manifest["local_concordance"] == [d["text"] for d in ca.local_concordance(result)],
      str(manifest.get("evidence"))[:200])
with zipfile.ZipFile(path) as z:
    edits_tsv = z.read("edits.tsv").decode().splitlines()
summaries = ca.event_summary(result)
check("bundle (spec 16): the manifest describes each applied event (type, size, place, nearest sites, origin, what "
      "sample 1 does there), edits.tsv has each event's label and the origin of its inserted bases",
      manifest["events"] == [d["line"] for d in summaries] and len(summaries) == 2
      and "origin of the inserted bases" in edits_tsv[0].split("\t") and "label" in edits_tsv[0].split("\t")
      and ca.event_summary(loaded) == [dict(d, event=x) for d, x in zip(summaries, [e for e in loaded.compat.alignment.events if e.apply])],
      str(manifest.get("events")))
damaged = os.path.join(tmp, "damaged" + bundle.EXTENSION)
with zipfile.ZipFile(path) as src, zipfile.ZipFile(damaged, "w") as dst:
    for item in src.namelist():
        body = src.read(item)
        dst.writestr(item, body.replace(b"Exon 3", b"Exon 4") if item == "summary.tsv" else body)
for what, target, expect in (("a damaged file", damaged, "doesn't match its checksum"),
                             ("not a zip", os.path.join(helpers.REPO, "README.md"), "not a readable zip")):
    try:
        bundle.load(target)
        check(f"bundle: {what} is refused", False)
    except bundle.BundleError as exc:
        check(f"bundle: {what} is refused, saying why", expect in str(exc), str(exc)[:120])

# --- the genes without their MANE Select transcript (hg19) -----------------------------------
real19 = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg19")
reader19 = refdb.Reader(real19) if os.path.exists(os.path.join(real19, "manifest.json")) else None
scores19 = refdb.ReferenceScores("hg19", str(helpers.HG19), reader=reader19, cache_dir=os.path.join(tmp, "cache19"),
                                 refseq_dir=REFSEQ, mane_dir=str(helpers.MANE_DIR))
record19 = scores19.predictor().annotator.ref_fasta
rhbg = table19.find("RHBG")[0]
rec = record19["chr" + rhbg.chrom]
fetch19 = lambda a, b: rec[a - 1:b].seq
seq19 = fetch19(rhbg.mask[0] - 3000, rhbg.mask[1] + 3000)
compat19 = ca.check_compatibility(seq19, rhbg, fetch19, len(rec))
r19 = ca.run_analysis("hg19", rhbg, compat19, scores19, fetch19, table19)
rows19 = ca.table_rows(r19)
check("RHBG run: the reference row is the fallback NM_, labelled with its accession, with its own cDNA and protein",
      rows19[0]["Type"] == "Reference (NM_001412175.1)" and rows19[0]["cDNA length"]
      and r19.transcripts[0].cdna.sequence and r19.transcripts[0].cdna.protein, str(rows19[0])[:200])
check("RHBG run: concordance and changes are against the fallback transcript (the patient is the reference here)",
      r19.transcripts[0].vs_ref[0].endswith("identical") and rows19[0]["Change vs Reference"] == "No change"
      and r19.reference_note == ("fallback", FALLBACK_NOTE), r19.transcripts[0].vs_ref[0])
path19 = bundle.save(r19, os.path.join(tmp, "rhbg"))
loaded19, manifest19, _ = bundle.load(path19)
check("RHBG: the results bundle's manifest records the note, and it comes back with the results",
      manifest19["reference_note"] == FALLBACK_NOTE and loaded19.reference_note == ("fallback", FALLBACK_NOTE)
      and ca.table_rows(loaded19) == rows19)
bop1 = table19.find("BOP1")[0]
try:
    ca.run_analysis("hg19", bop1, compat19, scores19, fetch19, table19)
    check("BOP1 (no NM_ at all): nothing is predicted", False)
except ValueError as exc:
    check("BOP1 (no NM_ at all): nothing is predicted, the explanation says why", "Nothing can be predicted" in str(exc))

# --- bedGraph ----------------------------------------------------------------------------
values = np.zeros(70_000)
values[100:103] = 0.5
values[40_000:40_010] = [0.1] * 5 + [0.12345] * 5
bg = os.path.join(tmp, "track.bedgraph.gz")
n = bedgraph.write(bg, "chr21", 1_000_001, values, description="test track")
check("bedGraph: equal runs merged, zero runs left out", n == 3 and bedgraph.runs(values, 1_000_001) ==
      [(1_000_100, 1_000_103, 0.5), (1_040_000, 1_040_005, 0.1), (1_040_005, 1_040_010, 0.1234)])
with TabixFile(bg) as tf:
    near = [line for line in tf.fetch("chr21", 1_040_001, 1_040_003)]
    far = [line for line in tf.fetch("chr21", 1_000_101, 1_000_101)]
    none = list(tf.fetch("chr22", 1, 10))
check("bedGraph: read back through its tabix index, by region", any(line.startswith("chr21\t1040000\t1040005\t0.1000")
      for line in near) and any("1000100\t1000103\t0.5000" in line for line in far) and not none, f"{near} {far}")
plan = result.plan
written = bedgraph.write(os.path.join(tmp, "acc.bedgraph.gz"), "chr" + gene.chrom, plan["start"], result.wt_values[:, 0])
with TabixFile(os.path.join(tmp, "acc.bedgraph.gz")) as tf:
    back = [line.split("\t") for line in tf.fetch("chr" + gene.chrom, plan["start"], plan["end"])]
rebuilt = np.zeros(plan["end"] - plan["start"] + 1)
for _, beg, end, value in back:
    rebuilt[int(beg) - plan["start"] + 1:int(end) - plan["start"] + 1] = float(value)
check(f"bedGraph: a real track ({written:,} lines) reads back to its values at 4 decimals",
      np.allclose(rebuilt, np.round(result.wt_values[:, 0].astype(np.float64), 4), atol=1e-9))


# Spec 17: an exon wholly inside inserted bases, a copy of a canonical exon or not.
from types import SimpleNamespace  # noqa: E402

from spliceai_pipeline.align import reverse_complement  # noqa: E402

dice = random.Random(17)
bases = "".join(dice.choice("ACGT") for _ in range(400))


def inserted(alt, strand="+"):
    hap, cmap = haplotype.build(bases, 1000, [Edit(1200, 1199, alt, "insertion")])
    gene = SimpleNamespace(canonical=SimpleNamespace(strand=strand, exons=[(1100, 1130), (1300, 1330)]))
    compat = SimpleNamespace(ref_seq=bases, region=(1000, 1399))
    r = SimpleNamespace(cmap=cmap, hap=hap, gene=gene, compat=compat)
    return ca.in_insertion(r, SimpleNamespace(exons=[(211, 241), (150, 160)], ref_exons=[None, (1150, 1160)]))


check("an exon in inserted bases that copy a canonical exon names it, in the transcript's own numbering",
      inserted(bases[90:180]) == {1: (1, 1200, 31)} and inserted(bases[90:180], "-") == {1: (2, 1200, 31)})
check("...and one in random bases, or in a copy of an intron, is no copy; an exon on the reference isn't listed",
      inserted(reverse_complement(bases)[:60]) == {1: (0, 1200, 31)} and inserted(bases[140:200]) == {1: (0, 1200, 31)})
helpers.finish()
