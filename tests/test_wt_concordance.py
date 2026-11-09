"""Spec 14 C1: the genome-wide wild-type concordance (wt_concordance.py) --
the categories on synthetic genes whose answer is known, the per-gene table
written and read back, the shipped tables, and a real chromosome from the
database.

    .venv\\Scripts\\python tests\\test_wt_concordance.py
"""
import os
import tempfile

import helpers
from helpers import check

helpers.setup()

import numpy as np  # noqa: E402

from spliceai_pipeline import refdb, transcripts, wt_concordance as wc  # noqa: E402
from spliceai_pipeline.mane import find_cached_summary  # noqa: E402

START = 1000
EXONS = ((1100, 1200), (1500, 1600), (2000, 2100), (2500, 2600))


def gene(other_exons=None):
    """A plus-strand gene with 4 canonical exons (and another NM_ if given)."""
    tx = lambda acc, exons: transcripts.Transcript(acc, "SYN", "chr1", "+", exons[0][0], exons[-1][1],
                                                   (1150, 2550), tuple(exons))
    c = tx("NM_000001.1", EXONS)
    ts = [c] + ([tx("NM_000002.1", other_exons)] if other_exons else [])
    return transcripts.Gene("SYN", "1", "+", 1100, 2600, EXONS, canonical=c, rule=transcripts.RULE_MANE, transcripts=ts)


def values(changes=()):
    """Every canonical site at 0.9, then {position: (kind, p)} on top."""
    v = np.zeros((2000, 2), np.float32)
    for a, b in EXONS[1:]:
        v[a - START, 0] = 0.9
    for a, b in EXONS[:-1]:
        v[b - START, 1] = 0.9
    for x, (kind, p) in dict(changes).items():
        v[x - START, kind] = p
    return v


same = wc.gene_concordance(gene(), values(), START)
check("synthetic: every site strong -> 4 of 4 exons identical, the whole transcript exact, every exon annotated",
      (same["identical"], same["whole_exact"], same["predicted_annotated"], same["categories"]) == (4, 1, 4, {}), same)
weak = wc.gene_concordance(gene(), values({2000: (0, 0.01), 2020: (0, 0.9)}), START)
check("synthetic: exon 3's acceptor below the site threshold -> 'canonical site below the site threshold'",
      weak["identical"] == 3 and weak["categories"] == {"canonical site below the site threshold": 1}, weak)
alt_changes = {1600: (1, 0.3), 1650: (1, 0.95)}
alt = wc.gene_concordance(gene(((1100, 1200), (1500, 1650), (2000, 2100), (2500, 2600))), values(alt_changes), START)
check("synthetic: exon 2's end moved to a donor another NM_ uses -> 'alternative exon', and that exon annotated",
      alt["categories"] == {"alternative exon": 1} and alt["predicted_annotated"] == 4, alt)
other = wc.gene_concordance(gene(), values(alt_changes), START)
check("synthetic: the same move with no NM_ using that donor -> 'other', that exon not annotated",
      other["categories"] == {"other": 1} and other["predicted_annotated"] == 3, other)
skip = wc.gene_concordance(gene(((1100, 1200), (2000, 2100), (2500, 2600))),
                           values({1500: (0, 0.06), 1600: (1, 0.06)}), START)
check("synthetic: exon 2 skipped, as another NM_ skips it -> 'alternative exon'",
      skip["categories"] == {"alternative exon": 1} and skip["predicted_exons"] == 3, skip)
edge = wc.gene_concordance(gene(), values({2500: (0, 0.3), 2550: (0, 0.95)}), START)
check("synthetic: the last exon starts elsewhere -> 'first or last exon edge'; first and last counted apart",
      edge["categories"] == {"first or last exon edge": 1} and edge["last_identical"] == 0
      and edge["first_identical"] == 1, edge)

tmp = tempfile.mkdtemp(prefix="spliceai_wt_concordance_")
path = os.path.join(tmp, "t.tsv.gz")
rows = [same, weak, alt]
for r, name in zip(rows, ("A", "B", "C")):
    r["gene"] = name
wc.write_table(rows, path)
back = wc.read_table(path)
check("the per-gene table is written and read back unchanged", [back[r["gene"]] for r in rows] == rows)
check("for_gene: a gene's row from the table of its build; None for an unknown gene or build",
      wc.for_gene("x", "B", tmp) is None and wc.for_gene("hg38", "NOT_A_GENE") is None)

for build in ("hg19", "hg38"):
    shipped = os.path.join(wc.DATA_DIR, wc.TABLE_NAME.format(build=build))
    table = wc.read_table(shipped) if os.path.isfile(shipped) else {}
    check(f"{build}: the shipped per-gene table is there, with every gene of a canonical transcript",
          len(table) > 17000 and table.get("BRCA2", {}).get("exons", 0) >= 27, f"{len(table)} genes")

spec = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "SpliceAI-VariantScoring.spec"),
            encoding="utf-8").read()
check("the PyInstaller .spec bundles both tables where for_gene looks (spliceai_pipeline/data)",
      all(f'("spliceai_pipeline/data/wt_concordance_{b}.tsv.gz", "spliceai_pipeline/data")' in spec
          for b in ("hg19", "hg38")))

db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
if os.path.exists(os.path.join(db, "manifest.json")):
    every = transcripts.load_genes("hg38", str(helpers.THIRDPARTY / "refseq"), find_cached_summary(str(helpers.MANE_DIR)))
    genes = [g for g in every if g.chrom == "21"]
    found, skipped = wc.run("hg38", refdb.Reader(db), genes)
    multi = [r for r in found if r["exons"] > 1]
    rate = sum(r["identical"] for r in multi) / sum(r["exons"] for r in multi)
    check("hg38 chr21 from the database: every gene with a canonical transcript, none skipped; > 90% of exons "
          "identical", len(found) == sum(g.canonical is not None for g in genes) and not skipped and rate > 0.9,
          f"{len(found)} genes, {rate:.1%}")
    shipped = wc.read_table(os.path.join(wc.DATA_DIR, wc.TABLE_NAME.format(build="hg38")))
    check("...and the shipped table has the same rows for them", all(shipped.get(r["gene"]) == {
        **r, "categories": r["categories"]} for r in found), "")
    # Spec 17 item 4: the same concordance at each minimum exon length.
    tp53 = next(g for g in every if g.name == "TP53")
    totals, shorter = wc.sweep(refdb.Reader(db), [tp53], values=(25, 22, 3))
    check("the sweep: TP53's exon 3 has 22 nt, so 9 of 11 exons are identical at a minimum of 25 and all 11 at 22 or 3",
          [totals[v]["identical"] for v in (25, 22, 3)] == [9, 11, 11] and totals[3]["short_found"][25] == 1
          and shorter[25]["internal"] == 1 and shorter[22]["internal"] == 0,
          str([(v, totals[v]["identical"]) for v in totals]))
    some = [g for g in genes if g.canonical and len(g.canonical.exons) > 1][:15]
    names = {g.name for g in some}
    totals, shorter = wc.sweep(refdb.Reader(db), some, values=(3, 25))
    multi = [r for r in found if r["gene"] in names]
    table = wc.sweep_report("hg38", totals, shorter, 60)
    check("...at the default it equals the plain concordance, and the report has a row per value",
          totals[3]["identical"] == sum(r["identical"] for r in multi) and "| 25 nt |" in table and "| 3 nt |" in table,
          f"{totals[3]['identical']} of {totals[3]['exons']}")
else:
    print("SKIP the real chromosome: no hg38 database here")
helpers.finish()
