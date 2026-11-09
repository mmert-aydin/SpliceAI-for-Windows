"""Spec 14 C6: events given as coordinates (coord_events.py, the tab's "From
coordinates" input) -- every form read, HGVS g. with its build, the reference
checks, the sample built from them equal to the haplotype of the same edits
and to what the aligner finds in it, and the breakpoint sensitivity run.

    .venv\\Scripts\\python tests\\test_coord_events.py
"""
import os

import helpers
from helpers import check

helpers.setup()

from spliceai_pipeline import cdna_analysis as ca, coord_events as ce, haplotype, refdb, specificity  # noqa: E402
from spliceai_pipeline.align import reverse_complement  # noqa: E402


def raises(fn, text):
    try:
        fn()
    except ValueError as exc:
        return text in str(exc)
    return False


# --- reading ---------------------------------------------------------------------------------
plain = ce.parse("""# a comment, then a blank line

deletion chr17:17,228,330-17,228,462
tandem duplication chr17:100-200
inverted tandem duplication 17:300-400
tandem duplication chr17:500-600 inverted
inversion chr17:700-800
insertion chr17:900 acgtn
deletion chr17:1000-1003 ACGT
""", "hg38")
check("plain forms: deletion, tandem duplication (direct, inverted two ways), inversion, insertion, deleted bases; "
      "comments and blank lines skipped", [(e.kind, e.chrom, e.start, e.end, e.sequence, e.ref) for e in plain] == [
          ("deletion", "17", 17228330, 17228462, "", ""), ("duplication", "17", 100, 200, "", ""),
          ("inverted duplication", "17", 300, 400, "", ""), ("inverted duplication", "17", 500, 600, "", ""),
          ("inversion", "17", 700, 800, "", ""), ("insertion", "17", 900, 900, "ACGTN", ""),
          ("deletion", "17", 1000, 1003, "", "ACGT")], str(plain))
hgvs = ce.parse("""NC_000017.11:g.17228330_17228462del
NC_000017.11:g.100_103delACGT
NC_000017.11:g.100_200dup
NC_000017.11:g.300_400inv
NC_000017.11:g.900_901insACGT
NC_000017.11:g.1000_1005delinsTTA
NC_000017.11:g.2000A>G""", "hg38")
check("HGVS g.: del (with and without bases), dup, inv, ins between two bases, delins, substitution",
      [(e.kind, e.start, e.end, e.sequence, e.ref) for e in hgvs] == [
          ("deletion", 17228330, 17228462, "", ""), ("deletion", 100, 103, "", "ACGT"), ("duplication", 100, 200, "", ""),
          ("inversion", 300, 400, "", ""), ("insertion", 900, 900, "ACGT", ""), ("delins", 1000, 1005, "TTA", ""),
          ("SNV", 2000, 2000, "G", "A")], str(hgvs))
check("HGVS g.: the accession names the build -- GRCh37's on hg38 is refused, an unknown one too",
      raises(lambda: ce.parse("NC_000017.10:g.100_200del", "hg38"), "GRCh37 (hg19) chromosome, but hg38")
      and ce.parse("NC_000017.10:g.100_200del", "hg19")[0].chrom == "17"
      and raises(lambda: ce.parse("NC_000099.1:g.1del", "hg38"), "isn't a GRCh37 or GRCh38"))
check("the chromosome accessions (UCSC chromAlias): chr1, chrX, chrM of both builds",
      ce.ACCESSIONS["NC_000001.11"] == ("hg38", "1") and ce.ACCESSIONS["NC_000001.10"] == ("hg19", "1")
      and ce.ACCESSIONS["NC_000023.11"] == ("hg38", "X") and ce.ACCESSIONS["NC_012920.1"] == ("hg38", "M")
      and ce.ACCESSIONS["NC_001807.4"] == ("hg19", "M") and len(ce.ACCESSIONS) == 50)
check("unreadable lines are refused, naming the line",
      raises(lambda: ce.parse("deletion chr17:100-200\nfoo chr17:1-2", "hg38"), "Line 2")
      and raises(lambda: ce.parse("insertion chr17:100", "hg38"), "the inserted bases")
      and raises(lambda: ce.parse("deletion chr17:200-100", "hg38"), "not after the end")
      and raises(lambda: ce.parse("NC_000017.11:g.100_105insA", "hg38"), "two bases it lies between"))

# --- against the reference (SAMSN1, hg38) ----------------------------------------------------------
REFSEQ, MANE = str(helpers.THIRDPARTY / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
table = ca.GeneTable("hg38", REFSEQ, MANE)
gene = table.find("SAMSN1")[0]
from pyfaidx import Fasta  # noqa: E402

record = Fasta(str(helpers.HG38))["chr21"]
fetch = lambda a, b: record[a - 1:b].seq
region = (gene.mask[0] - 5000, gene.mask[1] + 5000)
ref = lambda a, b: fetch(a, b).upper()
exons = list(reversed(gene.canonical.exons))
e5 = exons[4]
x = 14533713                                          # deep in intron 1
lines = [f"deletion chr21:{e5[1] - 19}-{e5[1] + 60}",                     # exon 5's acceptor (minus strand: its end)
         f"tandem duplication chr21:{x}-{x + 99}",
         f"inverted tandem duplication chr21:{x + 3000}-{x + 3079}",
         f"inversion chr21:{exons[2][0] - 100}-{exons[2][1] + 100}",
         f"insertion chr21:{x + 6000} ACGTTGCAACGTTGCA",
         f"NC_000021.9:g.{x + 9000}_{x + 9004}delinsGG",
         f"NC_000021.9:g.{x + 9500}{ref(x + 9500, x + 9500)}>{'A' if ref(x + 9500, x + 9500) != 'A' else 'C'}"]
compat = ca.from_coordinates("\n".join(lines), gene, "hg38", fetch, len(record))
kinds = [e.kind for e in compat.alignment.events]
check("from coordinates: every event is a difference, all applied, in reference order; the report says no "
      "alignment was used (the short-read breakpoint note is general information now, in the introduction)",
      sorted(kinds) == sorted(["deletion", "duplication", "duplication", "inversion", "insertion", "complex", "SNV"])
      and all(e.apply for e in compat.alignment.events) and "with no alignment" in compat.report["lines"][0]
      and compat.report["warnings"] == [] and not ca.chosen_check(compat)["problems"], kinds)
edits = [haplotype.Edit(e5[1] - 19, e5[1] + 60, "", "deletion"),
         haplotype.Edit(x + 100, x + 99, ref(x, x + 99), "duplication"),
         haplotype.Edit(x + 3080, x + 3079, reverse_complement(ref(x + 3000, x + 3079)), "insertion"),
         haplotype.Edit(exons[2][0] - 100, exons[2][1] + 100, "", "inversion"),
         haplotype.Edit(x + 6001, x + 6000, "ACGTTGCAACGTTGCA", "insertion"),
         haplotype.Edit(x + 9000, x + 9004, "GG", "complex"),
         haplotype.Edit(x + 9500, x + 9500, compat.alignment.events[-1].inserted, "SNV")]
expected, _ = haplotype.build(ref(*region), region[0], edits)
check("the sample built is the reference region with exactly those edits (a duplication copied after itself, an "
      "inverted one reverse-complemented, an insertion after its position)", compat.sequence == expected)
again = ca.check_compatibility(compat.sequence, gene, fetch, len(record))
for e in again.alignment.events:
    e.apply = True
rebuilt, _ = haplotype.apply_events(again.ref_seq, again.region[0], again.alignment.events)
check("...and aligning that sample back finds differences that rebuild it exactly", rebuilt == compat.sequence,
      [e.kind for e in again.alignment.events])
origin = lambda e: ca.event_origin(e, "chr21", compat.ref_seq, compat.region[0])
by_kind = {}
for e in compat.alignment.events:
    by_kind.setdefault(e.kind, []).append(e)
tandem, inverted = sorted(by_kind["duplication"], key=lambda e: e.start)
check("spec 16 item 2: events given as coordinates are labelled like aligned ones, and say where inserted bases "
      "come from -- a tandem copy, an inverted copy (in place for an inversion), short ones too short to trace; "
      "deletions and substitutions put none in",
      origin(tandem) == (f"tandem copy of chr21:{x:,}-{x + 99:,}", [(x, x + 99, 1)])
      and origin(inverted) == (f"inverted copy of chr21:{x + 3000:,}-{x + 3079:,}", [(x + 3000, x + 3079, -1)])
      and origin(by_kind["inversion"][0])[0].startswith("inverted copy of chr21:")
      and origin(by_kind["inversion"][0])[0].endswith("(in place)")
      and origin(by_kind["insertion"][0])[0] == "too short to trace (under 20 bases)"
      and origin(by_kind["deletion"][0]) == ("", []) == origin(by_kind["SNV"][0])
      and tandem.label == "tandem duplication" and inverted.label == "inverted duplication",
      [(e.kind, e.label, origin(e)[0]) for e in compat.alignment.events])

bad = lambda line: lambda: ca.from_coordinates(line, gene, "hg38", fetch, len(record))
wrong = "A" if ref(x, x + 3) != "AAAA" else "C"
check("the reference is checked: deleted bases that differ, a substitution's wrong reference base, a stretch "
      "outside the gene's region or on another chromosome, events that overlap -- each refused, saying why",
      raises(bad(f"deletion chr21:{x}-{x + 3} {wrong * 4}"), f"is {ref(x, x + 3)}, not {wrong * 4}")
      and raises(bad(f"NC_000021.9:g.{x}{'T' if ref(x, x) != 'T' else 'G'}>A"), "the reference at chr21:")
      and raises(bad(f"deletion chr21:{region[1] + 10}-{region[1] + 20}"), "outside the region analysed")
      and raises(bad("deletion chr17:100-200"), "but the gene is on chr21")
      and raises(bad(f"deletion chr21:{x}-{x + 50}\ninversion chr21:{x + 40}-{x + 90}"), "overlap")
      and raises(bad("# nothing"), "No event given"))

# --- breakpoint sensitivity -----------------------------------------------------------------------
db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=refdb.Reader(db) if os.path.isdir(db) else None,
                               refseq_dir=REFSEQ, mane_dir=str(helpers.MANE_DIR))
one = ca.from_coordinates(lines[0], gene, "hg38", fetch, len(record))
result = ca.run_analysis("hg38", gene, one, scores, fetch, table)
checked = specificity.breakpoint_sensitivity(result, scores)
b = checked[0]
check("breakpoint sensitivity: both breakpoints of the deletion moved by -50, -10, +10 and +50, the top transcript "
      "compared with the sample's own each time", len(checked) == 1 and len(b.tried) == 8
      and {(side, s) for side, s, _ in b.tried} == {(side, s) for side in ("start", "end") for s in specificity.SHIFTS}
      and b.text().startswith(f"The deletion at {e5[1] - 19:,} (80 bp): the top transcript stays the same for "),
      b.text())
r0 = one.region[0]
ins = haplotype.Edit(x + 101, x + 100, "ACGT", "insertion")
dup = haplotype.Edit(x + 100, x + 99, ref(x, x + 99), "duplication")
check("shifted breakpoints: an insertion has one (moved, same bases), a duplication's copied stretch follows its "
      "breakpoints, a substitution has none",
      specificity.shifted(ins, "start", 10, one.ref_seq, r0) == haplotype.Edit(x + 111, x + 110, "ACGT", "insertion")
      and specificity.shifted(ins, "end", 10, one.ref_seq, r0) is None
      and specificity.shifted(dup, "end", 10, one.ref_seq, r0) == haplotype.Edit(x + 110, x + 109, ref(x, x + 109),
                                                                                  "duplication")
      and specificity.shifted(haplotype.Edit(x, x, "A", "SNV"), "start", 10, one.ref_seq, r0) is None)
helpers.finish()
