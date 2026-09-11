"""VCF record normalization: left-alignment of indels and chrom-name standardization.

Implements the standard parsimony + left-alignment algorithm (Tan et al. 2015,
"Unified representation of genetic variants", Bioinformatics) -- the same
algorithm underlying `bcftools norm` / `vt normalize`. It is reimplemented here
in pure Python (against the reference FASTA via pyfaidx) because cyvcf2/htslib
provide no way to do this on Windows (see tabix.py docstring).

Alleles that can't be scored at all -- their chromosome isn't in the FASTA, or
their REF doesn't match the FASTA at that position (almost always a wrong
genome build) -- aren't dropped or allowed to crash the run: they come back
un-normalized with Variant.skip_reason set, and end up in the output as a row
whose source column says why.
"""
from dataclasses import dataclass

SKIP_UNKNOWN_CHROM = "skipped: chromosome not in reference FASTA"
SKIP_REF_MISMATCH = "skipped: REF doesn't match reference FASTA"


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    orig_chrom: str
    orig_pos: int
    orig_ref: str
    orig_alt: str
    allele_fraction: float = None
    skip_reason: str = None


def fasta_uses_chr_prefix(fasta):
    return any(str(k).lower().startswith("chr") for k in fasta.keys())


def standardize_chrom(chrom, want_chr_prefix):
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return "chr" + bare if want_chr_prefix else bare


# The mitochondrion goes by all of these (UCSC hg19/hg38: chrM; Ensembl/b37: MT).
_MITO_NAMES = ("chrM", "chrMT", "MT", "M")


def resolve_contig(chrom, contigs):
    """The FASTA's own name for a VCF chromosome, or None if the FASTA doesn't
    have it. contigs: set of the FASTA's sequence names. Accepts either naming
    style ("1" / "chr1") and every name of the mitochondrion."""
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    candidates = [chrom, bare, "chr" + bare]
    if bare.upper() in ("M", "MT"):
        candidates.extend(_MITO_NAMES)
    return next((c for c in candidates if c in contigs), None)


def make_base_fetcher(fasta):
    """Returns fetch(chrom, start_1based, end_1based) -> uppercase str, reading from pyfaidx."""

    def fetch(chrom, start, end):
        return fasta[chrom][start - 1:end].seq.upper()

    return fetch


def normalize_variant(fetch_bases, chrom, pos, ref, alt):
    """Left-align and minimize a single (chrom, pos, ref, alt) allele.

    chrom must already be in the naming convention used by `fetch_bases`
    (i.e. matching the FASTA's contig names). pos is 1-based.
    """
    ref = ref.upper()
    alt = alt.upper()

    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt

    # Right-trim shared trailing bases; if an allele empties out, extend both
    # alleles one base to the left from the reference.
    while True:
        if len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif len(ref) == 0 or len(alt) == 0:
            pos -= 1
            lb = fetch_bases(chrom, pos, pos)
            ref, alt = lb + ref, lb + alt
        else:
            break

    # Left-trim shared leading bases, keeping at least one base in each allele.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1

    # Indels: shift left through repeated sequence while both alleles end in
    # the same base -- Tan et al.'s "trim the shared last base, extend one
    # base left" step. (Comparing the shorter allele's first base with the
    # longer one's last instead is the same test for a plain insertion or
    # deletion, but moves a complex indel such as AT>GCA onto a different
    # variant.)
    if len(ref) != len(alt):
        while pos > 1 and ref[-1] == alt[-1]:
            prev_base = fetch_bases(chrom, pos - 1, pos - 1)
            pos -= 1
            ref = prev_base + ref[:-1]
            alt = prev_base + alt[:-1]

    return pos, ref, alt


def normalize_record(fasta, chrom, pos, ref, alts, allele_fractions=None, contigs=None):
    """Normalize one VCF data line, splitting multi-allelic ALTs into one Variant per allele.

    allele_fractions, if given, must be aligned 1:1 with alts (see
    vcfio.read_vcf_records / allele_fraction.extract_allele_fractions) -- each
    resulting Variant gets the fraction for the specific allele it came from.

    contigs: set of the FASTA's sequence names; computed when not given, but
    pass it in when normalizing many records.

    If the record can't be scored -- chromosome not in the FASTA, or REF not
    matching the FASTA there -- its alleles are returned un-normalized with
    skip_reason set (see module docstring).
    """
    if contigs is None:
        contigs = set(fasta.keys())
    if allele_fractions is None:
        allele_fractions = [None] * len(alts)

    fetch_bases = make_base_fetcher(fasta)
    fasta_chrom = resolve_contig(chrom, contigs)
    if fasta_chrom is None:
        skip_reason = SKIP_UNKNOWN_CHROM
    elif fetch_bases(fasta_chrom, pos, pos + len(ref) - 1) != ref.upper():
        skip_reason = SKIP_REF_MISMATCH
    else:
        skip_reason = None

    variants = []
    for alt, allele_fraction in zip(alts, allele_fractions):
        if alt in (".", "*") or alt.startswith("<"):
            continue
        if skip_reason:
            variants.append(Variant(
                chrom=fasta_chrom or chrom, pos=pos, ref=ref, alt=alt,
                orig_chrom=chrom, orig_pos=pos, orig_ref=ref, orig_alt=alt,
                allele_fraction=allele_fraction, skip_reason=skip_reason,
            ))
            continue
        npos, nref, nalt = normalize_variant(fetch_bases, fasta_chrom, pos, ref, alt)
        variants.append(Variant(
            chrom=fasta_chrom, pos=npos, ref=nref, alt=nalt,
            orig_chrom=chrom, orig_pos=pos, orig_ref=ref, orig_alt=alt,
            allele_fraction=allele_fraction,
        ))
    return variants
