"""Sequential VCF reading.

Only a linear pass over the input VCF is needed (every record has to be
normalized and scored), so this uses the stdlib gzip module -- bgzip output is
valid multi-member gzip and Python's gzip reader handles concatenated members
transparently. Indexed random access (needed only for the huge precomputed
score files) lives in tabix.py.
"""
import gzip

from .allele_fraction import extract_allele_fractions


def is_gzip(path):
    """By content, not by name: real VCFs are sometimes bgzipped under a plain
    .vcf name (and a .vcf.gz may not be compressed at all). The GUI's
    main_window._is_gzip checks the same two bytes."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _opener(path):
    # utf-8-sig: a VCF saved "UTF-8 with BOM" (some Windows editors) would
    # otherwise start with an invisible character that turns its first header
    # line into a broken data record and stops the run. errors="replace":
    # header text in some other encoding can't stop a run either.
    if is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "rt", encoding="utf-8-sig", errors="replace")


def zygosities(format_field, sample_field, n_alts):
    """Per ALT allele (1:1 with alts), what the first sample's GT says:
    "het" (one copy), "hom" (two), "hemi" (a haploid call carrying it, e.g.
    chrX in a male), or None (no GT, a missing call, or not carried). "1/2"
    is "het" for both alleles. vcf_bridge.genotype() spells the same rule
    out in words for one variant."""
    keys = format_field.split(":")
    if "GT" not in keys:
        return [None] * n_alts
    values = sample_field.split(":")
    i = keys.index("GT")
    calls = values[i].replace("|", "/").split("/") if i < len(values) else []
    result = []
    for allele in range(1, n_alts + 1):
        mine = calls.count(str(allele))
        result.append(None if not mine else "hemi" if len(calls) == 1 else "hom" if mine >= 2 else "het")
    return result


def read_vcf_records(path):
    """Yields (chrom, pos, ref, [alts], [allele_fractions], [zygosities]) for
    each data line of a VCF. allele_fractions is aligned 1:1 with alts (see
    allele_fraction.extract_allele_fractions) -- None for any allele whose
    fraction couldn't be determined from the record's INFO/FORMAT fields;
    zygosities likewise (see zygosities())."""
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
            yield chrom, pos, ref, alts, allele_fractions, zygosities(format_field, sample_field, len(alts))


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
