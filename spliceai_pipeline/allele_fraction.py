"""Extracts a per-allele allele fraction (AF) from a VCF record's INFO and/or
FORMAT/sample fields, for the results table's "Allele Fraction" column.

Not every caller reports AF the same way. Checked against the two real VCFs
available during development (both produced by
the same "Genomize" pipeline) declare AF (and AD) as valid FORMAT keys in
their header, but every actual data line's FORMAT column is just "GT:DP" --
AF/AD are only ever populated at the INFO level in practice for these files.
FORMAT/sample-level AF is still preferred over INFO whenever a record
genuinely has it (some other caller's output might), since it's the actual
per-sample measured value rather than a potentially cohort/site-level INFO
annotation -- it just happens not to change anything for the two files this
was verified against, since their FORMAT column never has it.

Falls back to computing AF from AD ("ref_depth,alt1_depth[,alt2_depth,...]",
the standard VCF convention: one more value than there are ALT alleles) as
alt_depth / (ref_depth + alt_depth), with the same FORMAT-preferred-over-INFO
rule, when no AF is present at all. A header's declared Number= for AD/AF is
not always trustworthy in practice (Genomize's own header declares AD as
Number=1, but the real field is a comma-separated pair) -- so this parses
whatever comma-separated content is actually there rather than trusting it.

Returns None (never 0.0) for an allele where neither AF nor a usable AD can
be found, so a VCF without this data just leaves the results column blank
rather than erroring or reporting a misleading zero.
"""


def _parse_info(info_field):
    """INFO field ('AD=1,2;AF=0.5;FLAG') -> {key: value_string}. Flag-only
    keys (no '=') are skipped -- never relevant here."""
    info = {}
    for kv in (info_field or "").split(";"):
        if "=" in kv:
            key, _, value = kv.partition("=")
            info[key] = value
    return info


def _parse_format_sample(format_field, sample_field):
    """FORMAT ('GT:DP:AF') + the first sample column's values ('0/1:877:0.447')
    -> {key: value_string}. Uses the first sample column only, consistent with
    this project's single-sample assumption elsewhere (e.g. vcf_info.py's
    "About this VCF" popup)."""
    if not format_field or not sample_field:
        return {}
    return dict(zip(format_field.split(":"), sample_field.split(":")))


def _per_allele_floats(raw_value, n_alts):
    """'0.447' or '0.1,0.2' -> a list of n_alts floats-or-None. Padded/
    truncated to n_alts rather than trusting the header's Number=A/R -- a
    caller that reports one value regardless of ALT count, or a header/data
    mismatch, shouldn't crash this, just leave the unmatched allele(s) None."""
    if not raw_value or raw_value == ".":
        return [None] * n_alts
    parts = [p.strip() for p in raw_value.split(",")]
    out = []
    for i in range(n_alts):
        if i < len(parts) and parts[i] not in ("", "."):
            try:
                out.append(float(parts[i]))
            except ValueError:
                out.append(None)
        else:
            out.append(None)
    return out


def _ad_to_allele_fractions(raw_ad, n_alts):
    """'ref_depth,alt1_depth[,alt2_depth,...]' -> a list of n_alts allele
    fractions. Standard VCF AD convention: the first value is the reference
    allele's depth, the rest are one per ALT allele in order. A single number
    with no comma isn't usable (no way to tell ref depth from alt depth), so
    that's treated the same as AD being absent entirely."""
    if not raw_ad or raw_ad == ".":
        return [None] * n_alts
    parts = [p.strip() for p in raw_ad.split(",")]
    if len(parts) < 2:
        return [None] * n_alts
    try:
        depths = [float(p) if p not in ("", ".") else None for p in parts]
    except ValueError:
        return [None] * n_alts
    ref_depth, alt_depths = depths[0], depths[1:]
    out = []
    for i in range(n_alts):
        if i >= len(alt_depths) or ref_depth is None or alt_depths[i] is None:
            out.append(None)
            continue
        total = ref_depth + alt_depths[i]
        out.append(alt_depths[i] / total if total > 0 else None)
    return out


def _prefer_format(format_values, info_values):
    """Combines two same-length per-allele lists, taking the FORMAT/sample
    value for each allele when it's present and falling back to INFO's value
    for that same allele otherwise."""
    return [f if f is not None else i for f, i in zip(format_values, info_values)]


def extract_allele_fractions(info_field, format_field, sample_field, n_alts):
    """Returns a list of n_alts Optional[float] allele fractions, aligned 1:1
    with the VCF record's (unsplit) ALT alleles in order -- matching the
    order normalize_record() later splits them into individual Variants."""
    info = _parse_info(info_field)
    fmt = _parse_format_sample(format_field, sample_field)

    af = _prefer_format(
        _per_allele_floats(fmt.get("AF"), n_alts),
        _per_allele_floats(info.get("AF"), n_alts),
    )
    if any(v is not None for v in af):
        return af

    return _prefer_format(
        _ad_to_allele_fractions(fmt.get("AD"), n_alts),
        _ad_to_allele_fractions(info.get("AD"), n_alts),
    )
