"""The three levels of benchmark.py's comparison, on exon chains made up here, both strands."""
import benchmark as b
from helpers import check, finish

PLUS = [(100, 199), (300, 399), (500, 599), (700, 799)]
MINUS = list(reversed(PLUS))


def g(text):
    return int(text[2:])


def compare(outcome, chain, canonical=PLUS, plus=True, marks=None):
    return b.match(outcome, chain, marks or {}, canonical, plus, g)


check("the canonical chain has no changes", not b.changes(PLUS, PLUS, True) and not b.changes(MINUS, MINUS, False))

skip = [(100, 199), (500, 599), (700, 799)]
check("an exon missing is a skip of that exon", b.changes(skip, PLUS, True) == {("skip", 2)})
check("the same skip matches, event and junction", compare({"skip": [2]}, skip) == (True, True))
check("a skip given by a position inside the exon", compare({"skip": ["g.350"]}, skip) == (True, True))
check("another exon's skip doesn't match", compare({"skip": [3]}, skip) == (False, False))
check("no change doesn't match", compare({"skip": [2]}, PLUS) == (False, False))

check("a range of exons, in either order", compare({"skip_range": ["g.250", "g.650"]}, [(100, 199), (700, 799)]) == (True, True)
      and compare({"skip_range": ["g.650", "g.250"]}, [(100, 199), (700, 799)]) == (True, True))

pe = [(100, 199), (300, 399), (440, 470), (500, 599), (700, 799)]
check("an exon in intron 2 is a pseudoexon there", b.changes(pe, PLUS, True) == {("pseudoexon", 2)})
check("the pseudoexon, exact", compare({"exons": [["g.440", "g.470"]]}, pe) == (True, True))
check("the right intron, other edges: event yes, junctions no", compare({"exons": [["g.440", "g.480"]]}, pe) == (True, False))
check("another intron: neither", compare({"exons": [["g.640", "g.670"]]}, pe) == (False, False))

pe_minus = list(reversed(pe))
check("minus strand: the intron counts in transcript order", b.changes(pe_minus, MINUS, False) == {("pseudoexon", 2)})
pe_early = [(700, 799), (620, 650), (500, 599), (300, 399), (100, 199)]
check("minus strand: an exon after the first exon is in intron 1", b.changes(pe_early, MINUS, False) == {("pseudoexon", 1)})
check("minus strand, exact, first and last in transcript order",
      compare({"exons": [["g.470", "g.440"]]}, pe_minus, MINUS, False) == (True, True))

donor = [(100, 199), (300, 420), (500, 599), (700, 799)]
check("a longer exon end is a cryptic donor", b.changes(donor, PLUS, True) == {("donor", 2)})
check("the donor, exact", compare({"donor": ["g.420"]}, donor) == (True, True))
check("a donor elsewhere in that intron: event yes, junctions no", compare({"donor": ["g.410"]}, donor) == (True, False))
short = [(100, 199), (300, 380), (500, 599), (700, 799)]
check("a donor inside the exon, exact", compare({"donor": ["g.380"]}, short) == (True, True))

acc_minus = [(700, 799), (500, 599), (300, 420), (100, 199)]
check("minus strand: the exon's higher end is its first base, an acceptor",
      b.changes(acc_minus, MINUS, False) == {("acceptor", 3)})
check("minus strand acceptor, exact", compare({"acceptor": ["g.420"]}, acc_minus, MINUS, False) == (True, True))

retained = [(100, 199), (300, 599), (700, 799)]
check("two exons in one is intron retention", b.changes(retained, PLUS, True) == {("retain", 2)})
check("retention, exact", compare({"retain": [2]}, retained) == (True, True))

dup = [(100, 199), (300, 399), None, (500, 599), (700, 799)]
marks = {("dup", 2): 100}
check("a copy of exon 2, whole", compare({"dup": [2]}, dup, marks=marks) == (True, True))
check("part of a copy: event yes, junctions no", compare({"dup": [2]}, dup, marks={("dup", 2): 60}) == (True, False))
check("no copy: neither", compare({"dup": [2]}, PLUS) == (False, False))

ins = [(100, 199), (300, 399), None, (500, 599), (700, 799)]
check("an exon of inserted bases, any size", compare({"inserted": 0}, ins, marks={("inserted", 0): [85]}) == (True, True))
check("an exon of inserted bases, its size", compare({"inserted": 85}, ins, marks={("inserted", 0): [85]}) == (True, True))
check("another size: event yes, junctions no", compare({"inserted": 90}, ins, marks={("inserted", 0): [85]}) == (True, False))

check("a pseudoexon known by its place: event yes, junctions not applicable",
      compare({"site": [["pseudoexon", "g.450"]]}, pe) == (True, None))
check("a place in another intron: no", compare({"site": [["pseudoexon", "g.650"]]}, pe)[0] is False)
check("a donor known by its place, past the exon", compare({"site": [["donor", "g.410"]]}, donor) == (True, None))
check("an acceptor known by its place, minus strand",
      compare({"site": [["acceptor", "g.410"]]}, acc_minus, MINUS, False) == (True, None))
out = b.compare([{"site": [["pseudoexon", "g.450"]]}], [(pe, {})], PLUS, True, g)
check("a case with no junction to compare stays None", out["exact_rank1"] is None and out["event_rank1"], str(out))
out = b.compare([{"site": [["pseudoexon", "g.450"]]}, {"skip": [2]}], [(skip, {})], PLUS, True, g)
check("a case with one junction known compares that one", out["exact_rank1"] is True, str(out))

check("a pseudoexon known by its last base", compare({"pe_end": "g.470"}, pe) == (True, True))
check("known by its first base", compare({"pe_start": "g.440"}, pe) == (True, True))
check("the wrong end: event yes, base no", compare({"pe_end": "g.460"}, pe) == (True, False))
check("minus strand: the last base is the low end", compare({"pe_end": "g.440"}, pe_minus, MINUS, False) == (True, True))

check("a change the WT prediction has already isn't the variant's",
      b.match({"exons": [["g.440", "g.470"]]}, pe, {}, PLUS, True, g, wild={("pseudoexon", 2)})[0] is False)
check("a WT change elsewhere doesn't matter",
      b.match({"exons": [["g.440", "g.470"]]}, pe, {}, PLUS, True, g, wild={("skip", 4)}) == (True, True))

predicted = [(PLUS, {}), (skip, {})]
out = b.compare([{"skip": [2]}], predicted, PLUS, True, g)
check("found at rank 2: the top counts, rank 1 doesn't",
      out["exact_top"] and out["event_top"] and not out["exact_rank1"] and not out["event_rank1"], str(out))
out = b.compare([{"skip": [3]}, {"skip": [2]}], [(skip, {})], PLUS, True, g)
check("any documented outcome counts", out["exact_rank1"] and out["outcome"] == 1, str(out))

import benchmark_report as r  # noqa: E402

lo, hi = r.clopper_pearson(7, 10)
check("exact interval of 7 of 10: 34.8-93.3%", round(lo, 3) == 0.348 and round(hi, 3) == 0.933, f"{lo}, {hi}")
check("exact interval of 0 of 10 starts at 0", r.clopper_pearson(0, 10)[0] == 0 and round(r.clopper_pearson(0, 10)[1], 3) == 0.308)
check("exact interval of 10 of 10 ends at 1", r.clopper_pearson(10, 10)[1] == 1 and round(r.clopper_pearson(10, 10)[0], 3) == 0.692)

finish()
