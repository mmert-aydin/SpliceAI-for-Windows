"""Parses a quick, informational summary out of a VCF's own header + content,
for the "About this VCF" popup (see main_window._on_about_vcf).

Independent of check_vcf_text's structural pass/fail check, but reuses its
exact variant-record counting logic (vcf_check.data_lines) rather than
recomputing a slightly different filter, so the two can never silently
disagree on what "N variants" means.

Purely informational -- summarize_vcf_text() never changes any GUI setting
itself (e.g. the Build dropdown), even when it detects a build mismatch;
main_window._on_about_vcf only *reports* one, leaving the fix to the user.
"""
import re
from dataclasses import dataclass, field

from .vcf_check import data_lines

# Maps a free-form "##reference=" value to one of the Build dropdown's
# builds; the first pattern that matches wins.
#
# (The original comments weren't recoverable from the compiled program;
# blank/comment lines like these keep every line where it was in the
# original, which Python 3.13 bakes into class definitions.)
_BUILD_ALIASES = [
    (re.compile(r"hg19|grch37|hs37d5|g1k_v37|\bb37\b", re.IGNORECASE), "hg19"),
    (re.compile(r"hg38|grch38|hs38|\bb38\b", re.IGNORECASE), "hg38"),
]

# chr1's length differs between the two builds, and most variant callers write
# every ##contig line with its length even when there's no ##reference line.
_CHR1_LENGTH_TO_BUILD = {"249250621": "hg19", "248956422": "hg38"}


def _match_build(value):
    for pattern, normalized in _BUILD_ALIASES:
        if pattern.search(value):
            return normalized
    return None


def _header_lines(text):
    """Just the header lines (up to #CHROM) -- cheap even for a huge VCF,
    since the build is re-detected whenever the input text changes."""
    end = text.find("\n#CHROM")
    return (text[:end] if end != -1 else text[:1_000_000]).splitlines()


def detect_vcf_build(text):
    """(build, evidence) for the build a VCF's own header names, or (None,
    None). Checks ##reference and ##assembly first, then chr1's length in the
    ##contig lines. evidence is a short, user-readable note of what matched."""
    lines = _header_lines(text)
    for key in ("##reference", "##assembly"):
        value = _header_value(lines, key)
        if value:
            build = _match_build(value)
            if build:
                shown = value if len(value) <= 40 else "..." + value[-37:]
                return build, f"{key[2:]}={shown}"
    for ln in lines:
        if ln.startswith("##contig=<"):
            fields = dict(kv.split("=", 1) for kv in ln[len("##contig=<"):].rstrip(">").split(",") if "=" in kv)
            if fields.get("ID") in ("1", "chr1"):
                build = _CHR1_LENGTH_TO_BUILD.get(fields.get("length", ""))
                if build:
                    return build, "chr1 length in its ##contig lines"
    return None, None


@dataclass
class VcfInfo:
    variant_count: int
    detected_build: str = None
    detected_build_raw: str = None
    file_date: str = None
    samples: list = field(default_factory=list)
    source: str = None
    contig_count: int = None


def _header_value(lines, key):
    """First "##key=value" header line's value (key includes the leading "##",
    e.g. "##reference"), or None if absent. Stops at #CHROM -- the header
    block is over by then, no need to scan into the (possibly huge) variant
    data that follows."""
    prefix = f"{key}="
    for ln in lines:
        if ln.startswith("#CHROM"):
            break
        if ln.startswith(prefix):
            return ln[len(prefix):].strip()
    return None


def summarize_vcf_text(text):
    lines = text.splitlines()

    raw_build = _header_value(lines, "##reference")
    detected_build = _match_build(raw_build) if raw_build else None

    header_line = next((ln for ln in lines if ln.startswith("#CHROM")), None)
    samples = []
    if header_line:
        cols = header_line.split("\t")
        if "FORMAT" in cols:
            samples = cols[cols.index("FORMAT") + 1:]

    return VcfInfo(
        variant_count=len(data_lines(lines)),
        detected_build=detected_build,
        detected_build_raw=raw_build,
        file_date=_header_value(lines, "##fileDate"),
        samples=samples,
        source=_header_value(lines, "##source"),
        contig_count=sum(1 for ln in lines if ln.startswith("##contig=")) or None,
    )
