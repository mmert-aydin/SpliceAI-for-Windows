"""Fast, structural-only VCF compatibility check.

Runs against whatever text is currently in the input box (loaded or pasted)
without invoking any part of the scoring pipeline, so the user gets a quick
yes/no before committing to a run that might take a long time. Checks the
same structural properties the pipeline itself depends on (a #CHROM header
with the mandatory columns, at least one parseable data record) plus a
sampled per-record sanity pass, without re-implementing normalize.py/vcfio.py
-- this is deliberately a cheaper, independent pre-check, not a dry run of
the real parser.
"""
from dataclasses import dataclass, field

from .i18n import t

REQUIRED_HEADER_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT")
SAMPLE_LIMIT = 2000
VALID_BASE_CHARS = set("ACGTNacgtn")


@dataclass
class VcfCheckResult:
    ok: bool
    summary: str
    details: list = field(default_factory=list)


def data_lines(lines):
    """Every non-empty, non-header line of a VCF -- i.e. one variant record.
    Factored out so anything else that needs a variant count (e.g.
    vcf_info.summarize_vcf_text's "About this VCF" popup) agrees with this
    check on exactly what counts as one, rather than recomputing its own
    slightly-different filter."""
    return [ln for ln in lines if ln.strip() and not ln.startswith("#")]


def check_vcf_text(text):
    if not text or not text.strip():
        return VcfCheckResult(False, t("Input is empty -- nothing to check."))

    lines = text.splitlines()

    header_line = next((ln for ln in lines if ln.startswith("#CHROM")), None)
    if header_line is None:
        return VcfCheckResult(False, t("No #CHROM header line found."))

    header_cols = header_line.split("\t")
    missing = [c for c in REQUIRED_HEADER_COLUMNS if c not in header_cols]
    if missing:
        return VcfCheckResult(False, t("Header is missing required column(s): {columns}", columns=", ".join(missing)))

    data = data_lines(lines)
    if not data:
        return VcfCheckResult(False, t("No variant records found (0 data lines after the header)."))

    sample = data[:SAMPLE_LIMIT]
    bad = []
    chrom_styles = set()
    for i, line in enumerate(sample, 1):
        fields = line.split("\t")
        if len(fields) < 5:
            bad.append(t("line {i}: fewer than 5 tab-separated fields ({n} found)", i=i, n=len(fields)))
            continue
        chrom, pos, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        if not pos.isdigit() or int(pos) < 1:
            bad.append(t("line {i}: POS is not a positive integer: {pos}", i=i, pos=repr(pos)))
            continue
        if not ref or any(c not in VALID_BASE_CHARS for c in ref):
            bad.append(t("line {i}: REF has unexpected characters: {ref}", i=i, ref=repr(ref)))
        if not alt:
            bad.append(t("line {i}: ALT is empty", i=i))
        chrom_styles.add(t("chr-prefixed") if chrom.lower().startswith("chr") else t("bare"))

    bad_ratio = len(bad) / len(sample)
    style_desc = t(" and ").join(sorted(chrom_styles)) if chrom_styles else t("unrecognized")
    sample_note = t(" (sampled first {n} of {total})", n=len(sample), total=len(data)) if len(data) > len(sample) else ""

    if bad_ratio > 0.5:
        return VcfCheckResult(
            False,
            t("{bad}/{n} sampled records failed basic parsing{note} -- this doesn't look like well-formed VCF data.",
              bad=len(bad), n=len(sample), note=sample_note),
            details=bad[:20],
        )

    summary = t("{n} variant record(s) found, header OK, chrom naming: {style}.", n=len(data), style=style_desc)
    if bad:
        summary += " " + t("{bad}/{n} sampled line(s) had minor issues{note} (see details).", bad=len(bad),
                           n=len(sample), note=sample_note)
    return VcfCheckResult(True, summary, details=bad[:20])
