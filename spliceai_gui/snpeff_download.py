"""Fetch and set up SnpEff: download the core distribution, unzip it, then
download the RefSeq-transcript database for a build -- mirrors
reference_download.py's download/progress/cancel pattern (reusing its
download_file() directly for the zip phase) but for a Java tool + its own
`download` subcommand instead of a plain reference FASTA.

SnpEff is optional (unlike the reference FASTA): this module is only ever
invoked from the GUI when the user has explicitly opted in via the "I will
use SnpEff" checkbox.
"""
import os
import shutil
import subprocess
import tempfile
import time
import zipfile

from spliceai_pipeline.snpeff import (
    JAVA_LOCALE_ARGS, SNPEFF_DB_BY_BUILD, SUBPROCESS_NO_WINDOW, JavaNotFoundError, find_java, snpeff_db_dir,
    snpeff_jar_path,
)

from .java_download import fetch_java
from .reference_download import DownloadCancelled, DownloadError, download_file

SNPEFF_ZIP_URL = "https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"


def _remove_if_exists(path):
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


def fetch_snpeff(dest_dir, build, on_progress=None, opener=None, should_cancel=None):
    """Full pipeline: Java (only if needed) -> download snpEff_latest_core.zip ->
    unzip -> download the RefSeq-transcript database for `build`
    (GRCh37.p13/GRCh38.p14 -- see spliceai_pipeline.snpeff's module docstring
    for why RefSeq and not the more commonly-referenced Ensembl GRCh3x.NN
    databases). Returns the jar path.

    on_progress(phase, done, total) with phase in {"java", "java_unzip",
    "download", "unzip", "database"}; total is None for indeterminate progress
    (unzipping and the database-download subcommand don't expose a byte-accurate
    total, same as fetch_reference's "index" phase).

    SnpEff needs Java 21+. If find_java() finds none -- no Java at all, or only
    an old one such as Java 8 -- a private Java runtime is downloaded first
    (java_download.fetch_java), so the user never installs Java by hand.
    Steps that are already done are skipped: an existing snpEff.jar isn't
    downloaded again, and neither is a database that's already installed for
    this build.

    The SnpEff zip is downloaded into a scratch folder under the system temp
    directory and unpacked from there (see java_download.fetch_java for why:
    antivirus ransomware protection can block creating a .zip in a Desktop
    folder). The scratch folder is always removed; a partial unzip/database
    download is left in place rather than guessed-and-deleted, since re-running
    this function safely overwrites/completes it either way.
    """
    if build not in SNPEFF_DB_BY_BUILD:
        raise DownloadError(f"Unknown build: {build!r} (expected one of {list(SNPEFF_DB_BY_BUILD)})")
    db = SNPEFF_DB_BY_BUILD[build]

    try:
        java_exe = find_java()
    except JavaNotFoundError:
        java_exe = fetch_java(on_progress=on_progress, opener=opener, should_cancel=should_cancel)

    os.makedirs(dest_dir, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix="spliceai_snpeff_")
    zip_path = os.path.join(work_dir, "snpEff_latest_core.zip")
    jar_path = snpeff_jar_path(dest_dir)

    try:
        if not os.path.isfile(jar_path):
            download_file(
                SNPEFF_ZIP_URL, zip_path,
                on_progress=lambda done, total: on_progress and on_progress("download", done, total),
                opener=opener, should_cancel=should_cancel, direct=True,
            )

            if should_cancel and should_cancel():
                raise DownloadCancelled("Cancelled before unzip.")
            if on_progress:
                on_progress("unzip", 0, None)
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(dest_dir)
            if on_progress:
                on_progress("unzip", 1, 1)

            if not os.path.isfile(jar_path):
                raise DownloadError(f"Unzip completed but {jar_path} is missing -- unexpected archive layout.")

        if os.path.isdir(snpeff_db_dir(dest_dir, build)):
            return jar_path

        if should_cancel and should_cancel():
            raise DownloadCancelled("Cancelled before database download.")
        if on_progress:
            on_progress("database", 0, None)
        log_fd, log_path = tempfile.mkstemp(prefix="snpeff_db_download_", suffix=".log")
        os.close(log_fd)
        try:
            with open(log_path, "w", encoding="utf-8") as log_fh:
                proc = subprocess.Popen(
                    [java_exe, *JAVA_LOCALE_ARGS, "-jar", jar_path, "download", "-v", db],
                    cwd=os.path.dirname(jar_path), stdout=log_fh, stderr=subprocess.STDOUT,
                    **SUBPROCESS_NO_WINDOW,
                )
                while proc.poll() is None:
                    if should_cancel and should_cancel():
                        proc.terminate()
                        try:
                            proc.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            proc.kill()
                        raise DownloadCancelled("Database download cancelled.")
                    time.sleep(0.3)
            if proc.returncode != 0:
                with open(log_path, encoding="utf-8", errors="replace") as log_fh:
                    output = log_fh.read()
                raise DownloadError(f"Database download failed (exit {proc.returncode}): {output[-2000:]}")
        finally:
            _remove_if_exists(log_path)
        if on_progress:
            on_progress("database", 1, 1)
    except (DownloadCancelled, DownloadError):
        raise
    except Exception as exc:
        raise DownloadError(f"SnpEff setup failed: {exc}") from exc
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return jar_path
