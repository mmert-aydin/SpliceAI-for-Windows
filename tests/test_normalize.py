"""Normalization regression: a deletion and an insertion in repeats (hg19).

An earlier version of normalize.py moved indels onto a different variant while
left-aligning them, and two such variants then scored 0.91 and 0.99. The two
here are synthetic, of the same shape (a 4-base unit deleted from a repeat, a
C inserted before a CCC run), deep in an ANKRD26 intron -- not patient data.

Each variant is written in every equivalent form found within 30 bases --
shifted right through its repeat, and with one shared base appended, as some
callers write them -- and every form must give the same sequence, normalize to
the same left-aligned variant, and score masked max 0.00 through the pipeline.

    .venv\\Scripts\\python tests\\test_normalize.py
"""
import helpers
from helpers import check

helpers.watchdog(900)
helpers.setup()

from pyfaidx import Fasta  # noqa: E402

from spliceai_pipeline.cli import run_pipeline_core  # noqa: E402
from spliceai_pipeline.normalize import make_base_fetcher, normalize_variant  # noqa: E402

CASES = [
    ("ANKRD26 deletion", "chr10", 27298042, "ACAGG", "A"),
    ("ANKRD26 insertion", "chr10", 27297926, "G", "GC"),
]

hg19 = Fasta(str(helpers.HG19), rebuild=False)
fetch = make_base_fetcher(hg19)


def equivalent_forms(chrom, pos, ref, alt, reach=30):
    """(pos, ref, alt) forms within `reach` bases that spell the same haplotype
    as the given one, read off the reference (so none is made up)."""
    lo, hi = pos - 10, pos + len(ref) + reach + 10
    window = fetch(chrom, lo, hi)

    def apply(p, r, a):
        return window[:p - lo] + a + window[p - lo + len(r):]

    target = apply(pos, ref, alt)
    forms = []
    for p in range(pos, pos + reach):
        r = window[p - lo:p - lo + len(ref)]
        a = target[p - lo:p - lo + len(alt)]
        if r[0] == a[0] and apply(p, r, a) == target:
            forms.append((p, r, a))
            nxt = window[p - lo + len(r)]
            forms.append((p, r + nxt, a + nxt))
    return forms, apply, target


lines = []
for gene, chrom, pos, ref, alt in CASES:
    check(f"{gene}: REF {ref} matches hg19 at {chrom}:{pos}", fetch(chrom, pos, pos + len(ref) - 1) == ref,
          fetch(chrom, pos - 3, pos + len(ref) + 3))
    forms, apply, target = equivalent_forms(chrom, pos, ref, alt)
    check(f"{gene}: equivalent forms found (shifted and padded)", len(forms) >= 2, f"{len(forms)}, e.g. {forms[1:4]}")
    check(f"{gene}: every form spells the same sequence", all(apply(*f) == target for f in forms))
    normalized = {normalize_variant(fetch, chrom, *f) for f in forms}
    check(f"{gene}: every form normalizes to {chrom}:{pos} {ref}>{alt}", normalized == {(pos, ref, alt)}, str(normalized))
    lines.extend(f"{chrom}\t{p}\t{gene}\t{r}\t{a}\t.\t.\t.\n" for p, r, a in forms)

rows = run_pipeline_core(helpers.write_vcf(lines), "hg19", "masked", str(helpers.HG19), on_progress=lambda *a: None)
for gene, chrom, pos, ref, alt in CASES:
    mine = [r for r in rows if r.gene == gene.split()[0] and (len(r.ref) > len(r.alt)) == (len(ref) > len(alt))]
    check(f"{gene}: every form scores masked max 0.00 (the old bug scored such indels 0.91 / 0.99)",
          mine and all(r.max_score == 0.0 and (r.pos, r.ref, r.alt) == (pos, ref, alt) for r in mine),
          f"{len(mine)} forms -> {sorted({(r.pos, r.ref, r.alt, r.max_score) for r in mine})}")

helpers.finish()
