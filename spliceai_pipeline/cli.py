"""Command-line entry point wiring normalize -> lookup -> score -> merge -> write.

run_pipeline_core() is the reusable core (no file I/O for the output table, just
returns rows, and reports progress via a callback instead of print()) so both
this CLI and spliceai_gui can drive the exact same pipeline logic without any
duplication. run_pipeline() is the thin CLI-specific wrapper that adds
print-based progress and writes the output file.
"""
import argparse
import sys
import time

from pyfaidx import Fasta

from .lookup import PrecomputedScores
from .merge import merge_variant
from .normalize import SKIP_REF_MISMATCH, SKIP_UNKNOWN_CHROM, normalize_record
from .vcfio import read_vcf_records, strip_vcf_extension
from .writer import delimiter_for_path, write_rows

BUILDS = ("hg19", "hg38")
MODES = ("raw", "masked")

# Above this share of variants whose REF doesn't match the FASTA, the run
# stops before scoring instead of producing a table of skipped rows: the
# right build gives ~0%, the wrong one ~70% (measured on a real hg38 VCF
# checked against hg19).
MAX_REF_MISMATCH_FRACTION = 0.10


class ReferenceMismatchError(ValueError):
    """The VCF and the reference FASTA don't belong together -- wrong genome
    build, or a FASTA whose chromosome names match none of the VCF's. Raised
    before any scoring, with a message fit to show directly to a user."""


def check_reference_fit(normalized, build, fasta_path, progress):
    """Raises ReferenceMismatchError if the normalized variants show the VCF
    doesn't fit the FASTA (see MAX_REF_MISMATCH_FRACTION); otherwise reports
    how many alleles will be listed as skipped, and why."""
    unknown = [v for v in normalized if v.skip_reason == SKIP_UNKNOWN_CHROM]
    mismatched = [v for v in normalized if v.skip_reason == SKIP_REF_MISMATCH]
    unknown_chroms = sorted({v.orig_chrom for v in unknown})
    chrom_list = ", ".join(unknown_chroms[:5]) + (", ..." if len(unknown_chroms) > 5 else "")
    on_reference = len(normalized) - len(unknown)

    if normalized and on_reference == 0:
        raise ReferenceMismatchError(
            f"None of the VCF's chromosomes ({chrom_list}) exist in the reference FASTA {fasta_path}. "
            "It probably names them differently (e.g. RefSeq 'NC_000001.11' instead of 'chr1' or '1'). "
            "Use a UCSC-style FASTA such as hg19.fa / hg38.fa (\"Download reference...\"). Nothing was scored."
        )
    if on_reference and len(mismatched) / on_reference > MAX_REF_MISMATCH_FRACTION:
        raise ReferenceMismatchError(
            f"{len(mismatched)} of {on_reference} variants ({len(mismatched) / on_reference:.0%}) have a REF "
            f"allele that doesn't match the reference FASTA at their position. That almost always means the "
            f"wrong genome build: {build} is selected, but this VCF is probably not {build}. Check the Build "
            f"setting and that the FASTA is {build}, then run again. Nothing was scored."
        )
    if unknown:
        progress(f"  {len(unknown)} allele(s) on chromosomes the FASTA doesn't have ({chrom_list}) -- listed as skipped")
    if mismatched:
        progress(f"  {len(mismatched)} allele(s) whose REF doesn't match the FASTA -- listed as skipped")


def format_elapsed(seconds):
    """H:MM:SS, whole seconds only (no decimals) -- e.g. 1574s -> "0:26:14".
    Shared by every elapsed-time display (CLI/GUI progress messages, the GUI's
    Summary dialog run duration) so they all read the same way."""
    total = int(round(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Score VCF variants with SpliceAI (precomputed lookup + live fallback).")
    p.add_argument("vcf", help="Input VCF file path (.vcf or .vcf.gz)")
    p.add_argument("--build", choices=BUILDS, required=True)
    p.add_argument("--mode", choices=MODES, required=True)
    p.add_argument("--precomputed-dir", default=None,
                   help="Directory containing spliceai_scores.*.vcf.gz(.tbi). Omit to skip precomputed "
                        "lookup entirely and live-score everything.")
    p.add_argument("--fasta", required=True, help="Reference FASTA (faidx-indexed)")
    p.add_argument("--use-snpeff", action="store_true",
                   help="Also annotate each variant with its transcript (NM_ ID) and coding position "
                        "(HGVS.c) via SnpEff. Off by default -- SnpEff is optional, unlike --fasta. "
                        "Requires SnpEff + the database for --build to already be set up (see "
                        "--snpeff-dir).")
    p.add_argument("--snpeff-dir", default=None,
                   help="SnpEff install root (contains a 'snpEff' subfolder with snpEff.jar and its data/ "
                        "databases). Default: <project root>/snpeff")
    p.add_argument("--mane-dir", default=None,
                   help="Folder containing a downloaded NCBI MANE Select summary file (see "
                        "spliceai_pipeline/mane.py), used to prefer the canonical transcript when a "
                        "variant overlaps several. Optional -- if not found, canonical selection falls "
                        "back to curated (NM_/NR_) over predicted (XM_/XR_) only. Default: <project "
                        "root>/mane")
    p.add_argument("--skip-precomputed", action="store_true",
                   help="Only live-score variants NOT found in the precomputed file. Default: live-score "
                        "everything, for QC/comparison against precomputed. Has no effect if "
                        "--precomputed-dir is omitted (everything is live-scored either way).")
    p.add_argument("-o", "--output", default=None,
                   help="Output table path (.csv or .tsv). Default: <vcf, extension stripped>-calculated.tsv")
    return p.parse_args(argv)


def run_pipeline_core(vcf_path, build, mode, fasta_path, precomputed_dir=None, skip_precomputed=False,
                      use_snpeff=False, snpeff_dir=None, mane_dir=None,
                      on_progress=None, control=None):
    """Normalize -> lookup (if precomputed_dir given) -> live-score -> merge.

    control: optional control.RunControl -- pause/end from the GUI, checked
    between variants; an ended run raises control.RunCancelled.

    Returns a list of merge.ScoreRow. Writing the output table is left to the
    caller (see write_rows) so this can be reused by a GUI that wants the rows
    in memory rather than a file on disk.

    on_progress(message, current=None, total=None): called for every status
    update. current/total are only given together, when there's a concrete
    progress fraction for the current phase (e.g. "2000/120888 looked up");
    otherwise this is just a phase-boundary announcement. Defaults to print().
    """
    def progress(message, current=None, total=None):
        if on_progress:
            on_progress(message, current, total)
        else:
            print(message, flush=True)

    def checkpoint():
        if control is not None:
            control.check()

    # Seconds spent in each phase, reported in the final "Done" message.
    # (Original comments weren't recoverable from the compiled program; blank and
    # comment lines like these keep line numbers where they were.)
    phase_seconds = {}

    t0 = time.time()
    progress(f"Loading reference FASTA: {fasta_path}")
    fasta = Fasta(fasta_path, rebuild=False)
    phase_seconds["fasta_load"] = time.time() - t0

    precomputed = None
    if precomputed_dir:
        progress(f"Opening precomputed scores ({build}/{mode}): {precomputed_dir}")
        precomputed = PrecomputedScores(precomputed_dir, build, mode)
    else:
        progress("No precomputed-data directory given -- every variant will be live-scored.")

    t0 = time.time()
    progress(f"Reading and normalizing variants from {vcf_path}")
    normalized = []
    contigs = set(fasta.keys())
    for chrom, pos, ref, alts, allele_fractions in read_vcf_records(vcf_path):
        checkpoint()
        normalized.extend(normalize_record(fasta, chrom, pos, ref, alts, allele_fractions, contigs=contigs))
    phase_seconds["normalize"] = time.time() - t0
    progress(f"  {len(normalized)} variant alleles after splitting multi-allelic records")
    check_reference_fit(normalized, build, fasta_path, progress)

    snpeff_annotations = {}
    if use_snpeff:
        # The snpeff module is imported here, only when SnpEff annotation was
        # requested, rather than at the top of the file.
        # (Original comments weren't recoverable from the compiled program; blank and
        # comment lines like these keep line numbers where they were.)
        checkpoint()
        progress("Running SnpEff annotation...")
        from .snpeff import JavaNotFoundError, SnpEffNotSetUpError, run_snpeff_annotation
        t0 = time.time()
        try:
            unique_variants = {(v.chrom, v.pos, v.ref, v.alt) for v in normalized if not v.skip_reason}
            snpeff_annotations = run_snpeff_annotation(
                unique_variants, build, snpeff_dir=snpeff_dir, mane_dir=mane_dir,
                on_progress=lambda m: progress(f"  {m}"),
                should_cancel=(lambda: control.cancelled) if control is not None else None,
            )
        except (JavaNotFoundError, SnpEffNotSetUpError, RuntimeError) as exc:
            # On these errors the run continues without SnpEff annotations
            # (SnpEff is optional, unlike --fasta).
            # (Original comments weren't recoverable from the compiled program; blank and
            # comment lines like these keep line numbers where they were.)
            progress(f"  SnpEff annotation skipped: {exc}")
            snpeff_annotations = {}
        phase_seconds["snpeff"] = time.time() - t0
        progress(f"  SnpEff annotation took {format_elapsed(phase_seconds['snpeff'])}")

    precomputed_hits = [None] * len(normalized)
    if precomputed is not None:
        progress("Looking up precomputed scores...")
        n_found = 0
        t_start = time.time()
        for n, v in enumerate(normalized, 1):
            checkpoint()
            hit = None if v.skip_reason else precomputed.lookup(v.chrom, v.pos, v.ref, v.alt)
            if hit is not None:
                n_found += 1
            precomputed_hits[n - 1] = hit
            if n % 2000 == 0 or n == len(normalized):
                elapsed = time.time() - t_start
                progress(f"  looked up {n}/{len(normalized)} ({n_found} found, {format_elapsed(elapsed)} elapsed)",
                         n, len(normalized))
        phase_seconds["precomputed_lookup"] = time.time() - t_start
        progress(f"  {n_found}/{len(normalized)} found in precomputed scores")
        precomputed.close()

    # (Original comments weren't recoverable from the compiled program; blank and
    # comment lines like these keep line numbers where they were.)
    to_live_score = [
        i for i in range(len(normalized))
        if not normalized[i].skip_reason and (precomputed_hits[i] is None or not skip_precomputed)
    ]

    live_results = {}
    skip_reasons = {}
    if to_live_score:
        progress(f"Running live SpliceAI scoring on {len(to_live_score)} variant(s) (CPU)...")
        from .score import LiveScorer
        t_start = time.time()
        scorer = LiveScorer(fasta_path, build, mode)
        phase_seconds["live_scoring_model_load"] = time.time() - t_start
        progress(f"Live scorer ready (TensorFlow imported, models loaded, {format_elapsed(phase_seconds['live_scoring_model_load'])}).")
        t_start = time.time()
        for n, idx in enumerate(to_live_score, 1):
            v = normalized[idx]
            checkpoint()
            try:
                live_results[idx], reason = scorer.score_with_reason(v.chrom, v.pos, v.ref, v.alt)
            except Exception as exc:
                # One unscorable variant must never cost a whole (possibly
                # hour-long) run: it's listed with the error instead.
                live_results[idx], reason = [], f"skipped: scoring error ({type(exc).__name__}: {exc})"
            if reason:
                skip_reasons[idx] = reason
            if n % 10 == 0 or n == len(to_live_score):
                elapsed = time.time() - t_start
                progress(f"  live-scored {n}/{len(to_live_score)} ({format_elapsed(elapsed)} elapsed)",
                         n, len(to_live_score))
        phase_seconds["live_scoring"] = time.time() - t_start
        if skip_reasons:
            progress(f"  SpliceAI skipped {len(skip_reasons)} variant(s) -- the reason is in the source column")
    else:
        progress("No variants require live scoring.")

    t0 = time.time()
    progress("Merging results...")
    rows = []
    for i, v in enumerate(normalized):
        snpeff_ann = snpeff_annotations.get((v.chrom, v.pos, v.ref, v.alt))
        rows.extend(merge_variant(
            v.chrom, v.pos, v.ref, v.alt, precomputed_hits[i], live_results.get(i), snpeff_ann=snpeff_ann,
            allele_fraction=v.allele_fraction, skip_reason=v.skip_reason or skip_reasons.get(i),
        ))
    phase_seconds["merge"] = time.time() - t0

    total = sum(phase_seconds.values())
    breakdown = ", ".join(f"{name}={format_elapsed(seconds)}" for name, seconds in phase_seconds.items())
    progress(f"Done ({len(rows)} rows). Phase timing: {breakdown} (total {format_elapsed(total)})")
    return rows


def run_pipeline(args):
    out_path = args.output or strip_vcf_extension(args.vcf) + "-calculated.tsv"
    delimiter = delimiter_for_path(out_path)

    rows = run_pipeline_core(
        args.vcf, args.build, args.mode, args.fasta,
        precomputed_dir=args.precomputed_dir,
        skip_precomputed=args.skip_precomputed,
        use_snpeff=args.use_snpeff,
        snpeff_dir=args.snpeff_dir,
        mane_dir=args.mane_dir,
    )

    print(f"Writing {len(rows)} rows to {out_path}", flush=True)
    write_rows(rows, out_path, delimiter=delimiter)
    return out_path


def main(argv=None):
    args = parse_args(argv)
    run_pipeline(args)


if __name__ == "__main__":
    sys.exit(main())
