"""Sequential VCF reading.

Only a linear pass over the input VCF is needed (every record has to be
normalized and scored), so this uses the stdlib gzip module -- bgzip output is
valid multi-member gzip and Python's gzip reader handles concatenated members
transparently. Indexed random access (needed only for the huge precomputed
score files) lives in tabix.py.
"""
import gzip

from .allele_fraction import extract_allele_fractions


def _opener(path):
    # utf-8-sig: a VCF saved "UTF-8 with BOM" (some Windows editors) would
    # otherwise start with an invisible character that turns its first header
    # line into a broken data record and stops the run. errors="replace":
    # header text in some other encoding can't stop a run either.
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "rt", encoding="utf-8-sig", errors="replace")


def read_vcf_records(path):
    """Yields (chrom, pos, ref, [alts], [allele_fractions]) for each data line
    of a VCF. allele_fractions is aligned 1:1 with alts (see
    allele_fraction.extract_allele_fractions) -- None for any allele whose
    fraction couldn't be determined from the record's INFO/FORMAT fields."""
    with _opener(path) as fh:
        for line in fh:
            # Blank lines (e.g. pasted text ending in an empty line) aren't
            # records -- skipped rather than failing the whole run.
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            alts = alt.split(",")
            info_field = fields[7] if len(fields) > 7 else ""
            format_field = fields[8] if len(fields) > 8 else ""
            sample_field = fields[9] if len(fields) > 9 else ""
            allele_fractions = extract_allele_fractions(info_field, format_field, sample_field, len(alts))
            yield chrom, pos, ref, alts, allele_fractions


def count_vcf_records(path):
    n = 0
    with _opener(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n


def strip_vcf_extension(path):
    """Strips a trailing .vcf.gz/.vcf/.bgz (case-insensitive) from a path or
    filename -- used to derive an output filename from an input VCF's name
    (e.g. "my_sample.vcf.gz" -> "my_sample") without leaving the compound
    ".vcf.gz" extension half-stripped (a plain path.rsplit(".vcf", 1) does
    the wrong thing on ".vcf.gz" if a caller isn't careful about ordering, so
    this checks the longest/compound suffix first). Returns path unchanged if
    none of these extensions match.
    """
    for suffix in (".vcf.gz", ".vcf.bgz", ".vcf"):
        if path.lower().endswith(suffix):
            return path[:-len(suffix)]
    return path
