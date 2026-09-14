"""Combine precomputed-lookup and live-scoring results into unified output rows."""
from dataclasses import dataclass
from typing import Optional


@dataclass
class ScoreRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    DS_AG: Optional[float]
    DP_AG: Optional[int]
    DS_AL: Optional[float]
    DP_AL: Optional[int]
    DS_DG: Optional[float]
    DP_DG: Optional[int]
    DS_DL: Optional[float]
    DP_DL: Optional[int]
    max_score: Optional[float]
    source: str
    snpeff_transcript: Optional[str] = None
    snpeff_hgvs_c: Optional[str] = None
    # Where in the gene this variant falls, from the same SnpEff entry the
    # transcript came from: "Intron 5/26", "Exon 3/27", "5' UTR", "Intergenic".
    snpeff_region: Optional[str] = None
    allele_fraction: Optional[float] = None

    @property
    def variant_id(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    @property
    def alteration(self):
        """Display-only ref>alt shorthand (e.g. "T>C") for the GUI's combined
        "Alteration" column -- like variant_id, this is derived, never a stored
        field, so it doesn't affect CSV export/import schemas (ref and alt stay
        available as their own columns everywhere)."""
        return f"{self.ref}>{self.alt}"


def parse_annotation(annotation):
    """'ALT|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL' -> (alt, gene, [4 DS floats-or-None], [4 DP ints-or-None]).

    Both live scoring (spliceai.utils.get_delta_scores) and the official
    precomputed score files use this same 10-field format -- verified against
    a real downloaded precomputed VCF's SpliceAI= INFO field, not assumed.
    Older or malformed annotations with only 6 fields still parse fine, just
    with all-None DP values.
    """
    parts = annotation.split("|")
    alt, gene = parts[0], parts[1]
    ds = [None if v == "." else float(v) for v in parts[2:6]]
    dp = [None if v == "." else int(v) for v in parts[6:10]] if len(parts) >= 10 else [None, None, None, None]
    return alt, gene, ds, dp


def rows_from_annotations(chrom, pos, ref, alt, annotations, source):
    rows = []
    for annotation in annotations:
        a_alt, gene, (ds_ag, ds_al, ds_dg, ds_dl), (dp_ag, dp_al, dp_dg, dp_dl) = parse_annotation(annotation)
        scored = [v for v in (ds_ag, ds_al, ds_dg, ds_dl) if v is not None]
        max_score = max(scored) if scored else None
        rows.append(ScoreRow(
            chrom, pos, ref, a_alt, gene,
            ds_ag, dp_ag, ds_al, dp_al, ds_dg, dp_dg, ds_dl, dp_dl,
            max_score, source,
        ))
    return rows


_SCORE_FIELDS = ("DS_AG", "DP_AG", "DS_AL", "DP_AL", "DS_DG", "DP_DG", "DS_DL", "DP_DL")


def _same_scores(row_a, row_b):
    return all(getattr(row_a, f) == getattr(row_b, f) for f in _SCORE_FIELDS)


def _snpeff_for_gene(snpeff_ann, gene, n_genes):
    """The SnpEff (transcript_id, hgvs_c, region) for one output row, or None.

    snpeff_ann is snpeff.run_snpeff_annotation()'s per-variant value --
    {gene_name: (transcript_id, hgvs_c, region), ..., None: overall pick} -- or
    a plain tuple, which applies to every row. A row gets its
    own gene's annotation. If SnpEff has none under that gene name, a variant
    with a single gene still gets the overall pick (SpliceAI's and SnpEff's gene
    sets don't always use the same symbol); with several genes that row is
    left blank rather than showing another gene's transcript.
    """
    if isinstance(snpeff_ann, tuple):
        return snpeff_ann
    if gene and gene in snpeff_ann:
        return snpeff_ann[gene]
    if n_genes <= 1:
        return snpeff_ann.get(None)
    return None


def merge_variant(chrom, pos, ref, alt, precomputed_annotations, live_annotations, snpeff_ann=None,
                  allele_fraction=None, skip_reason=None):
    """precomputed_annotations: list[str] if found in the precomputed file, else None.
    live_annotations: list[str] if live scoring was run for this variant, else None
    (note: an empty list means live scoring ran but found no overlapping gene).
    snpeff_ann: optional SnpEff annotation for this (chrom, pos, ref, alt) --
    see _snpeff_for_gene for how it's matched to each row's gene.
    skip_reason: why this variant couldn't be scored (see normalize.py and
    score.SKIP_REASONS), shown as the placeholder row's source instead of
    "none" -- which then only ever means "scored, but no gene overlaps".
    allele_fraction: optional float, stamped onto every row the same way as
    snpeff_ann -- it's a property of the input variant call itself (from the
    VCF's own AF/AD), not of which gene/source produced a given row.

    A variant found in the precomputed data still gets live-scored too whenever
    skip_precomputed is off (the default -- see cli.py), specifically so the two
    can be compared for QC. But when both scored the same gene and landed on
    identical DS/DP values -- the overwhelmingly common case, since it's the
    same deterministic model computing the same variant -- emitting both as
    separate rows isn't a comparison, it's just the same result twice; only the
    live-sourced row is dropped in that case, per gene. If they disagree for a
    given gene (the actual QC-relevant case this dual-scoring exists to catch),
    both rows are kept so the mismatch stays visible rather than being hidden.

    Returns the list of ScoreRow rows to emit for this variant -- one row per
    (source, gene) combination that isn't a redundant precomputed/live
    duplicate, or a single placeholder row with source='none' if neither source
    produced anything, so every input variant is still represented in the output.
    """
    precomputed_rows = rows_from_annotations(chrom, pos, ref, alt, precomputed_annotations, "precomputed") if precomputed_annotations else []
    live_rows = rows_from_annotations(chrom, pos, ref, alt, live_annotations, "live") if live_annotations else []

    precomputed_by_gene = {row.gene: row for row in precomputed_rows}
    rows = list(precomputed_rows)
    for live_row in live_rows:
        matching_precomputed = precomputed_by_gene.get(live_row.gene)
        if matching_precomputed is not None and _same_scores(matching_precomputed, live_row):
            continue
        rows.append(live_row)

    if not rows:
        rows.append(ScoreRow(chrom, pos, ref, alt, "", None, None, None, None, None, None, None, None, None,
                             skip_reason or "none"))
    if snpeff_ann is not None:
        n_genes = len({row.gene for row in rows if row.gene})
        for row in rows:
            picked = _snpeff_for_gene(snpeff_ann, row.gene, n_genes)
            if picked is not None:
                # A 2-tuple is accepted so a caller (or an older saved result)
                # that has no region still merges.
                row.snpeff_transcript, row.snpeff_hgvs_c = picked[0], picked[1]
                row.snpeff_region = picked[2] if len(picked) > 2 else None
    if allele_fraction is not None:
        for row in rows:
            row.allele_fraction = allele_fraction
    return rows
