"""NCBI MANE Select transcript data -- used by snpeff.py to prefer the one
clinically-referenced transcript per gene over an arbitrary/predicted one when
a variant overlaps several annotated transcripts.

MANE (Matched Annotation from NCBI and EBI) Select is a single, expert-curated
RefSeq/Ensembl transcript pair per protein-coding gene, agreed upon by NCBI and
Ensembl as *the* representative transcript. It's defined once against GRCh38
coordinates, but the RefSeq accession it names is a property of the transcript
itself (the same mRNA/NM_ entry), not of the genome build -- so the same MANE
Select list applies equally to hg19-coordinate annotations, matched by base
accession only (ignoring the version suffix, since SnpEff's GRCh37.p13 vs
GRCh38.p14 databases can reference different versions of the identical
transcript, e.g. NM_007194.3 vs NM_007194.4).

Deliberately excludes "MANE Plus Clinical" (a secondary list of additional
clinically-relevant transcripts beyond the one Select pick per gene) -- only
"MANE Select" rows are used, per spec.

Not vendored in this repo (see .gitignore's "mane/" entry) -- same reasoning
as the reference FASTA and SnpEff: this is a separately-downloaded, standard
public reference file, not project code.
"""
import csv
import gzip
import os
import re

from . import net
from .app_paths import app_root_dir

MANE_CURRENT_DIR_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"

# Default download location: a "mane" folder next to the app (see
# app_paths.app_root_dir()), which is where the pipeline reads it from.
#
# (Original comments weren't recoverable from the compiled program; blank and
# comment lines like these keep line numbers where they were.)
DEFAULT_MANE_DIR = os.path.join(app_root_dir(), "mane")

_SUMMARY_FILENAME_RE = re.compile(r"MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz")


class ManeDownloadError(Exception):
    pass


def find_current_summary_filename(opener=None):
    """Discovers the exact current summary filename (e.g.
    "MANE.GRCh38.v1.5.summary.txt.gz") by reading NCBI's own directory listing,
    rather than hardcoding a version number that will silently go stale the
    next time MANE releases a new version."""
    opener = opener or net.urlopen
    try:
        with opener(MANE_CURRENT_DIR_URL) as response:
            listing = response.read().decode("utf-8", errors="replace")
    except Exception as exc:
        raise ManeDownloadError(net.download_error_message(exc, f"Could not list {MANE_CURRENT_DIR_URL}.")) from exc
    match = _SUMMARY_FILENAME_RE.search(listing)
    if not match:
        raise ManeDownloadError(
            "Could not find a MANE.GRCh38.vX.X.summary.txt.gz filename in the directory listing at "
            f"{MANE_CURRENT_DIR_URL} -- NCBI may have changed their naming scheme."
        )
    return match.group(0)


def fetch_mane_summary(dest_dir=None, opener=None):
    """Downloads the current MANE Select summary file into dest_dir (default
    DEFAULT_MANE_DIR). Small (~1MB gzipped, ~19k gene rows) compared to the
    FASTA/SnpEff downloads, so this is a plain one-shot download, no chunked
    progress reporting needed. Returns the local file path.
    """
    dest_dir = dest_dir or DEFAULT_MANE_DIR
    opener = opener or net.urlopen
    filename = find_current_summary_filename(opener=opener)
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, filename)
    tmp_path = dest_path + ".part"
    url = MANE_CURRENT_DIR_URL + filename
    try:
        with opener(url) as response, open(tmp_path, "wb") as out:
            out.write(response.read())
        os.replace(tmp_path, dest_path)
    except Exception as exc:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise ManeDownloadError(net.download_error_message(exc, f"Failed to download {url}.")) from exc
    return dest_path


def is_summary_filename(path):
    """True for a MANE summary file name, e.g. MANE.GRCh38.v1.5.summary.txt.gz
    -- the only kind find_cached_summary() picks up."""
    return bool(_SUMMARY_FILENAME_RE.fullmatch(os.path.basename(path)))


def sibling_mane_dir(snpeff_dir):
    """The "mane" folder beside a given SnpEff install, or None.

    Setup puts snpeff, java and mane next to each other. When the app is
    running from somewhere else than the install -- from source, or after the
    folder is moved -- this is where the summary actually is, so it is worth
    looking here before reporting that MANE isn't set up. Same reasoning as
    snpeff.sibling_java_dir.
    """
    if not snpeff_dir:
        return None
    candidate = os.path.join(os.path.dirname(os.path.abspath(snpeff_dir)), "mane")
    return candidate if os.path.isdir(candidate) else None


def find_cached_summary(mane_dir=None):
    """Returns the path to an already-downloaded summary file under mane_dir,
    or None if none is present. Doesn't download anything."""
    mane_dir = mane_dir or DEFAULT_MANE_DIR
    if not os.path.isdir(mane_dir):
        return None
    candidates = [
        f for f in os.listdir(mane_dir)
        if _SUMMARY_FILENAME_RE.fullmatch(f)
    ]
    if not candidates:
        return None
    # Version numbers sort correctly as strings here, so the last one is the
    # newest.
    candidates.sort()
    return os.path.join(mane_dir, candidates[-1])


def load_mane_select(path):
    """Parses a MANE summary file into {base_refseq_accession: (full_accession,
    gene_symbol)} for MANE Select rows only (MANE Plus Clinical is excluded --
    see module docstring). base_refseq_accession has no version suffix (e.g.
    "NM_007194", not "NM_007194.4"), matching how SnpEff's ANN transcript IDs
    are compared against this table."""
    mane = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("MANE_status") != "MANE Select":
                continue
            full_accession = row.get("RefSeq_nuc") or ""
            if not full_accession:
                continue
            base_accession = full_accession.split(".")[0]
            mane[base_accession] = (full_accession, row.get("symbol", ""))
    return mane


def load_cached_mane_select(mane_dir=None):
    """Convenience: find_cached_summary() + load_mane_select(), or {} if no
    cached file is present. This is the "gracefully degrade" entry point
    snpeff.py uses -- MANE preference is skipped (falls back to curated-over-
    predicted only) rather than failing the whole SnpEff step, consistent with
    how a missing SnpEff database itself is handled."""
    path = find_cached_summary(mane_dir)
    if not path:
        return {}
    return load_mane_select(path)
