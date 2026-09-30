"""SnpEff annotation step: runs SnpEff on the normalized variant set to get a
transcript ID + coding-position (HGVS.c) annotation per variant, merged into
the same chrom/pos/ref/alt row as SpliceAI's own scores.

When a variant overlaps multiple annotated transcripts, the transcript/HGVS
reported is chosen by priority (see _select_transcript_annotation): (1) NCBI
MANE Select, matched by base RefSeq accession -- version-agnostic, since
SnpEff's GRCh37.p13 database can reference an older version of the identical
transcript than GRCh38.p14 does, and MANE Select's chosen accession is a
property of the transcript itself, not of genome build coordinates, so the
same MANE table applies to both hg19 and hg38 -- (2) any curated NM_/NR_
transcript over a predicted XM_/XR_ one, (3) otherwise whichever entry SnpEff
lists first (unchanged fallback). See spliceai_pipeline/mane.py.

Uses RefSeq-transcript SnpEff databases (GRCh37.p13 / GRCh38.p14), not the
more commonly-referenced Ensembl-based GRCh3x.NN databases -- those annotate
with Ensembl transcript IDs (ENST...), not the NM_ RefSeq IDs this pipeline
needs to report.

SnpEff itself is a standalone Java tool, not a Python package -- it is
deliberately NOT vendored in this git repo and NOT listed in requirements.txt
(same reasoning as the reference FASTA and precomputed scores: multi-hundred-
MB-to-GB downloads don't belong in version control). Its default install
location is a plain "snpeff/" folder directly under the project root
(DEFAULT_SNPEFF_DIR), colocated with the project the same way data/ holds the
FASTA/precomputed files -- gitignored (see .gitignore), never committed.
Configurable via spliceai_gui.config's "snpeff_dir" the same way
fasta_path/precomputed_dir are. Setup (downloading SnpEff + the two
databases) is a one-time step, either manual (see check_snpeff_setup()'s
error message for the exact commands) or via the GUI's "Download SnpEff..."
button (spliceai_gui.snpeff_download_dialog).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

from .app_paths import app_root_dir
from .control import RunCancelled

# Default install root: a "snpeff" folder next to the app (see
# app_paths.app_root_dir()); the jar and databases live in its "snpEff"
# subfolder (see snpeff_jar_path()).
DEFAULT_SNPEFF_DIR = os.path.join(app_root_dir(), "snpeff")

# Oldest Java major version SnpEff will run on. find_java() skips anything
# older -- e.g. a Java 8 install that shadows a newer JDK on PATH -- instead
# of letting SnpEff fail later with an unhelpful Java error.
#
MIN_JAVA_VERSION = 21

SNPEFF_DB_BY_BUILD = {
    "hg19": "GRCh37.p13",
    "hg38": "GRCh38.p14",
}

# Where "Download SnpEff..." puts a private Java runtime when the PC has no
# Java 21+ of its own (see spliceai_gui/java_download.py): a "java" folder next
# to the app, like "snpeff" and "mane". find_java() looks here first.
DEFAULT_JAVA_DIR = os.path.join(app_root_dir(), "java")

# SnpEff upper-cases type names with the JVM's default locale; under Turkish
# rules "string" becomes "STRİNG" and its type lookup fails. Every java launch
# pins an English locale.
JAVA_LOCALE_ARGS = ["-Duser.language=en", "-Duser.country=US"]

# java.exe is a console program: without this, every launch from the windowed
# app would flash a console window.
SUBPROCESS_NO_WINDOW = {"creationflags": subprocess.CREATE_NO_WINDOW} if sys.platform == "win32" else {}


# Folders searched (newest-named subfolder first) for bin\java.exe when the
# java on PATH is missing or too old: Eclipse Adoptium (Temurin), Oracle, and
# Microsoft builds of OpenJDK.
_KNOWN_JAVA_LOCATIONS = [
    r"C:\Program Files\Eclipse Adoptium",
    r"C:\Program Files\Java",
    r"C:\Program Files\Microsoft",
]


class JavaNotFoundError(RuntimeError):
    """No Java installation new enough for SnpEff could be found."""


class SnpEffNotSetUpError(RuntimeError):
    """SnpEff itself (jar) or the database for the requested build isn't installed."""


_java_version_cache = {}


def java_major_version(java_exe):
    """Returns the major version as an int (handles both the old '1.8.0_x' and
    current '21.0.x' java -version formats), or None if it couldn't be determined.

    Cached per executable: the GUI re-checks SnpEff's status on every input
    change, and starting a JVM each time would make typing lag."""
    if java_exe in _java_version_cache:
        return _java_version_cache[java_exe]
    try:
        result = subprocess.run(
            [java_exe, "-version"], capture_output=True, text=True, timeout=15, **SUBPROCESS_NO_WINDOW
        )
    except (OSError, subprocess.SubprocessError):
        return None
    text = result.stderr or result.stdout or ""
    version = None
    m = re.search(r'version "(\d+)\.(\d+)', text)
    if m:
        major, minor = int(m.group(1)), int(m.group(2))
        version = minor if major == 1 else major
    else:
        m = re.search(r'version "(\d+)"', text)
        version = int(m.group(1)) if m else None
    if version is not None:
        _java_version_cache[java_exe] = version
    return version


def private_java_exe(java_dir=None):
    """java.exe of the app's own Java runtime (DEFAULT_JAVA_DIR unless given),
    newest folder first, or None if there isn't one. Doesn't check the version."""
    java_dir = java_dir or DEFAULT_JAVA_DIR
    if not os.path.isdir(java_dir):
        return None
    for entry in sorted(os.listdir(java_dir), reverse=True):
        exe = os.path.join(java_dir, entry, "bin", "java.exe")
        if os.path.isfile(exe):
            return exe
    return None


def sibling_java_dir(snpeff_dir):
    """The "java" folder beside a given SnpEff install, or None.

    Setup puts snpeff, java and mane next to each other, so whoever pointed
    SnpEff at a folder almost always has the matching Java one level up from
    it. Looking there as well as next to the app is what makes SnpEff work
    when the app is run from somewhere else than the install -- from source,
    or after the install folder is moved -- instead of reporting no Java while
    a perfectly good one sits beside the jar.
    """
    if not snpeff_dir:
        return None
    candidate = os.path.join(os.path.dirname(os.path.abspath(snpeff_dir)), "java")
    return candidate if os.path.isdir(candidate) else None


def find_java(snpeff_dir=None):
    """Locates a java executable satisfying MIN_JAVA_VERSION, without modifying
    system PATH/environment. Checks, in order: the app's own runtime in
    DEFAULT_JAVA_DIR (installed by "Download SnpEff..." when needed), the java
    folder beside snpeff_dir when one is given (see sibling_java_dir), 'java'
    on PATH, then common Windows JDK install locations -- so this still works
    when an older Java shadows a newer one on PATH.

    Raises JavaNotFoundError, with a message fit to show directly to a user,
    if nothing suitable is found.
    """
    candidates = []
    private = private_java_exe()
    if private:
        candidates.append(private)
    beside = private_java_exe(sibling_java_dir(snpeff_dir))
    if beside:
        candidates.append(beside)
    on_path = shutil.which("java")
    if on_path:
        candidates.append(on_path)
    for base in _KNOWN_JAVA_LOCATIONS:
        if os.path.isdir(base):
            for entry in sorted(os.listdir(base), reverse=True):
                exe = os.path.join(base, entry, "bin", "java.exe")
                if os.path.isfile(exe):
                    candidates.append(exe)

    checked = []
    for exe in candidates:
        version = java_major_version(exe)
        checked.append(f"{exe} (version {version if version is not None else 'unknown'})")
        if version is not None and version >= MIN_JAVA_VERSION:
            return exe

    raise JavaNotFoundError(
        f"No Java {MIN_JAVA_VERSION}+ found -- SnpEff needs it. Click \"Download SnpEff...\" and the app "
        "downloads its own copy of Java automatically (nothing to install by hand).\nChecked: "
        f"{', '.join(checked) if checked else '(no java found on PATH or common install locations)'}"
    )


def snpeff_jar_path(snpeff_dir):
    """snpeff_dir is the configurable install root (DEFAULT_SNPEFF_DIR unless the
    user pointed it elsewhere) -- the actual jar sits one level down, in a
    "snpEff" subfolder, because that's the layout the official snpEff_latest_core.zip
    unpacks into; this mirrors it exactly rather than reshuffling files after download.

    Always returns an absolute path: run_snpeff_annotation() runs the SnpEff
    subprocess with cwd set to this path's own parent directory (so SnpEff can
    find its config file alongside the jar), and a relative jar_path would get
    silently re-resolved against that new cwd instead of the caller's original
    one -- e.g. --snpeff-dir snpeff from the project root would have Java look
    for snpeff/snpEff/snpeff/snpEff/snpEff.jar and fail with "Unable to access
    jarfile". Made absolute here once so every caller (check_snpeff_setup, the
    subprocess command line) is immune regardless of whether snpeff_dir itself
    was given as relative or absolute.
    """
    return os.path.abspath(os.path.join(snpeff_dir, "snpEff", "snpEff.jar"))


def snpeff_db_dir(snpeff_dir, build):
    db = SNPEFF_DB_BY_BUILD.get(build)
    return os.path.join(snpeff_dir, "snpEff", "data", db) if db else None


def check_snpeff_setup(snpeff_dir, build):
    """Raises SnpEffNotSetUpError with a clear message if snpEff.jar or the
    database for this build isn't installed under snpeff_dir. Doesn't check
    Java -- call find_java() separately for that."""
    db = SNPEFF_DB_BY_BUILD.get(build)
    if db is None:
        raise SnpEffNotSetUpError(f"No SnpEff database is configured for build {build!r}.")
    jar_path = snpeff_jar_path(snpeff_dir)
    if not os.path.isfile(jar_path):
        raise SnpEffNotSetUpError(
            f"SnpEff isn't installed under {snpeff_dir} -- expected {jar_path}. Download "
            "snpEff_latest_core.zip from https://pcingola.github.io/SnpEff/download/ and unzip it into "
            f"{snpeff_dir}."
        )
    db_dir = snpeff_db_dir(snpeff_dir, build)
    if not os.path.isdir(db_dir):
        raise SnpEffNotSetUpError(
            f"SnpEff database {db!r} for build {build!r} isn't installed -- expected {db_dir}. Run: java -jar "
            f"{jar_path} download {db}"
        )


_VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _write_variants_vcf(variants, out_path):
    """variants: iterable of (chrom, pos, ref, alt) tuples, already in the same
    normalized/left-aligned representation used everywhere else in this
    pipeline, so the merge key back in cli.py matches exactly."""
    with open(out_path, "w", newline="") as fh:
        fh.write(_VCF_HEADER)
        for chrom, pos, ref, alt in variants:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")


def _ann_entries(info_field):
    """Parses every comma-separated entry of an ANN= INFO field into its
    pipe-separated fields, dropping entries too short to be usable. Empty list
    if there's no ANN field at all.

    ANN entry format (SnpEff spec): Allele|Annotation|Annotation_Impact|
    Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|
    HGVS.p|cDNA.pos/cDNA.length|CDS.pos/CDS.length|AA.pos/AA.length|Distance|ERRORS
    """
    m = re.search(r"ANN=([^\t;]+)", info_field)
    if not m:
        return []
    entries = []
    for entry in m.group(1).split(","):
        fields = entry.split("|")
        if len(fields) >= 10:
            entries.append(fields)
    return entries


# Where in the gene's structure a variant falls, from SnpEff's effect term.
# Checked in order, first match wins: an effect is often several terms joined
# by "&" (e.g. "splice_region_variant&intron_variant"), and the structural part
# is what this column is for -- a splice-region variant sitting in an intron is
# reported as the intron it is in, with its number.
_REGION_RULES = (
    ("intron", "Intron"),
    ("splice_acceptor", "Intron"),
    ("splice_donor", "Intron"),
    ("5_prime_utr", "5' UTR"),
    ("3_prime_utr", "3' UTR"),
    ("upstream_gene", "Upstream"),
    ("downstream_gene", "Downstream"),
    ("intergenic", "Intergenic"),
    ("exon", "Exon"),
    # Coding consequences don't say "exon" in the term, but that is where they
    # are, and SnpEff's Rank then counts exons.
    ("missense", "Exon"),
    ("synonymous", "Exon"),
    ("stop_", "Exon"),
    ("start_", "Exon"),
    ("frameshift", "Exon"),
    ("inframe_", "Exon"),
    ("initiator_codon", "Exon"),
    ("coding_sequence", "Exon"),
)

# Regions where SnpEff's Rank ("3/27") counts something, so it is worth showing.
_RANKED_REGIONS = ("Intron", "Exon")


def describe_region(effect, rank=""):
    """A short "where in the gene" label from SnpEff's Annotation and Rank
    fields, e.g. "Intron 5/26", "Exon 3/27", "5' UTR", "Intergenic".

    effect is the raw ANN Annotation field, which may join several terms with
    "&"; rank is the raw ANN Rank field ("3/27" or empty). An effect SnpEff
    reports that isn't in _REGION_RULES is shown as itself, tidied up, rather
    than dropped -- an unexplained blank would be worse than an unfamiliar term.
    """
    effect = (effect or "").strip()
    if not effect:
        return None
    lowered = effect.lower()
    label = next((name for term, name in _REGION_RULES if term in lowered), None)
    if label is None:
        # e.g. "non_coding_transcript_variant" -> "Non coding transcript variant"
        first = lowered.split("&")[0].replace("_variant", "").replace("_", " ").strip()
        return first[:1].upper() + first[1:] if first else None
    rank = (rank or "").strip()
    if rank and label in _RANKED_REGIONS:
        return f"{label} {rank}"
    return label


def _select_transcript_annotation(info_field, mane_select, gene=None):
    """Extracts (transcript_id, hgvs_c, region) from an ANN= INFO field, picking among
    however many transcripts SnpEff annotated for this variant (only those of
    `gene`, if given) by priority:

      1. NCBI MANE Select, matched by base RefSeq accession (version-agnostic).
         mane_select: {base_accession: (full_accession, gene_symbol)}, e.g. from
         spliceai_pipeline.mane.load_cached_mane_select() -- pass {} (or None)
         to skip this and go straight to (2), e.g. when no MANE data has been
         downloaded yet; this is a graceful degrade, not a hard failure.
      2. Any curated transcript (NM_/NR_ prefix) over a predicted one (XM_/XR_).
      3. Whichever entry SnpEff lists first (the original, un-prioritized
         behavior), if neither of the above found anything better.

    region (e.g. "Intron 5/26") comes from the same chosen entry as the
    transcript, never from a different one -- an exon number that belonged to
    some other transcript would be worse than none.

    Returns (None, None, None) if there's no ANN field or no parseable entry
    at all.
    """
    entries = _ann_entries(info_field)
    if gene is not None:
        entries = [fields for fields in entries if fields[3] == gene]
    if not entries:
        return None, None, None

    def transcript_id_of(fields):
        return fields[6] or None

    def picked(fields):
        # fields: 1 = Annotation (effect), 6 = Feature_ID, 8 = Rank, 9 = HGVS.c
        return transcript_id_of(fields), fields[9] or None, describe_region(fields[1], fields[8])

    if mane_select:
        for fields in entries:
            transcript_id = transcript_id_of(fields)
            if transcript_id and transcript_id.split(".")[0] in mane_select:
                return picked(fields)

    for fields in entries:
        transcript_id = transcript_id_of(fields)
        if transcript_id and transcript_id.startswith(("NM_", "NR_")):
            return picked(fields)

    return picked(entries[0])


def _gene_annotations(info_field, mane_select):
    """{gene_name: (transcript_id, hgvs_c, region)} for every gene SnpEff annotated,
    plus the key None for the overall pick across all genes -- so a variant
    overlapping two genes gets each row's own gene's transcript, never the
    other gene's (see merge._snpeff_for_gene). Empty if nothing usable."""
    result = {}
    nothing = (None, None, None)
    overall = _select_transcript_annotation(info_field, mane_select)
    if overall != nothing:
        result[None] = overall
    for gene in {fields[3] for fields in _ann_entries(info_field) if fields[3]}:
        picked = _select_transcript_annotation(info_field, mane_select, gene=gene)
        if picked != nothing:
            result[gene] = picked
    return result


def run_snpeff_annotation(variants, build, snpeff_dir=None, java_exe=None, mane_dir=None, on_progress=None,
                          should_cancel=None):
    """variants: iterable of (chrom, pos, ref, alt) tuples (normalized
    representation, deduplicated by the caller if needed).
    snpeff_dir: install root (see module docstring); defaults to
    DEFAULT_SNPEFF_DIR, but callers pass through the GUI/CLI-configured path so
    this respects wherever the user actually has SnpEff installed.
    mane_dir: where to look for an already-downloaded MANE Select summary (see
    spliceai_pipeline.mane); defaults to mane.DEFAULT_MANE_DIR. If none is
    found, canonical-transcript selection just skips the MANE-preference step
    and falls back to curated-over-predicted -- this is a soft dependency, not
    a hard requirement, same philosophy as SnpEff itself being optional.

    Returns {(chrom, pos, ref, alt): {gene_name: (transcript_id, hgvs_c), ...,
    None: (transcript_id, hgvs_c)}} -- one pick per annotated gene plus the
    overall pick under None (see _gene_annotations) -- for every variant
    SnpEff produced a usable ANN= entry for; variants with no annotation (e.g.
    intergenic, or outside any transcript in this database) are simply absent
    from the returned dict rather than mapped to (None, None). When a variant
    overlaps multiple transcripts, which one is reported follows the priority
    described in _select_transcript_annotation's docstring.

    Raises JavaNotFoundError / SnpEffNotSetUpError if the prerequisites aren't
    in place, or RuntimeError if the SnpEff subprocess itself fails.
    """
    def progress(msg):
        if on_progress:
            on_progress(msg)

    snpeff_dir = snpeff_dir or DEFAULT_SNPEFF_DIR
    check_snpeff_setup(snpeff_dir, build)
    db = SNPEFF_DB_BY_BUILD[build]
    jar_path = snpeff_jar_path(snpeff_dir)
    java_exe = java_exe or find_java(snpeff_dir)

    from .mane import load_cached_mane_select
    mane_select = load_cached_mane_select(mane_dir)
    progress(
        f"Using MANE Select data ({len(mane_select)} genes) for canonical transcript preference."
        if mane_select else
        "No MANE Select data found -- canonical transcript preference will fall back to curated-over-predicted only."
    )

    variants = list(variants)
    if not variants:
        return {}

    with tempfile.TemporaryDirectory(prefix="spliceai_snpeff_") as tmpdir:
        in_vcf = os.path.join(tmpdir, "variants.vcf")
        out_vcf = os.path.join(tmpdir, "annotated.vcf")
        _write_variants_vcf(variants, in_vcf)

        progress(f"Running SnpEff ({db}) on {len(variants)} variant(s)...")
        cmd = [java_exe, *JAVA_LOCALE_ARGS, "-Xmx4g", "-jar", jar_path, "-noStats", "-noLog", db, in_vcf]
        err_path = os.path.join(tmpdir, "snpeff_stderr.txt")
        # Polled rather than subprocess.run(), so ending the run (should_cancel)
        # stops SnpEff at once; stderr goes to a file, so a chatty SnpEff can
        # never fill a pipe and hang.
        with open(out_vcf, "w", encoding="utf-8") as out_fh, open(err_path, "w", encoding="utf-8") as err_fh:
            proc = subprocess.Popen(
                cmd, stdout=out_fh, stderr=err_fh, cwd=os.path.dirname(jar_path), **SUBPROCESS_NO_WINDOW,
            )
            while True:
                try:
                    proc.wait(timeout=0.3)
                    break
                except subprocess.TimeoutExpired:
                    if should_cancel and should_cancel():
                        proc.kill()
                        proc.wait()
                        raise RunCancelled("Run ended by the user.")
        if proc.returncode != 0:
            with open(err_path, encoding="utf-8", errors="replace") as fh:
                stderr = fh.read()
            raise RuntimeError(f"SnpEff failed (exit {proc.returncode}): {stderr[-2000:]}")

        annotations = {}
        with open(out_vcf, encoding="utf-8") as fh:
            for line in fh:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                chrom, pos, ref, alt, info = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
                gene_annotations = _gene_annotations(info, mane_select)
                if gene_annotations:
                    annotations[(chrom, pos, ref, alt)] = gene_annotations

        progress(f"SnpEff annotated {len(annotations)}/{len(variants)} variant(s).")
        return annotations
