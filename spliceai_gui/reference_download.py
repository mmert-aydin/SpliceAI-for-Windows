"""Fetch and prepare a reference FASTA from UCSC's goldenPath: download,
decompress, and build a .fai index -- so a user without a local reference
doesn't have to find/build one manually. samtools may not be available on
Windows (a real setup issue hit earlier in this project), so indexing uses
pyfaidx (pure Python, pip-installable) instead of shelling out to `samtools
faidx`.

Split into small, independently-testable functions (download_file /
decompress_gz / build_fai_index / fetch_reference) so tests can substitute a
tiny fake HTTP response instead of the real ~1GB download, while still
exercising the real URL construction, progress reporting, decompression, and
indexing logic end to end.
"""

import gzip
import os
import time
import urllib.request


REFERENCE_URLS = {
    "hg19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
}

CHUNK_SIZE = 1024 * 1024


class DownloadError(Exception):
    pass


class DownloadCancelled(DownloadError):
    pass


def _remove_if_exists(path):
    if path and os.path.exists(path):
        try:
            os.remove(path)
        except OSError:
            pass


# ~10 s in total: long enough for an antivirus on-access scan of a freshly
# downloaded file, short enough that a genuinely read-only folder still fails
# promptly.
RENAME_ATTEMPTS = 50
RENAME_DELAY_S = 0.2


def _replace_with_retry(src, dst, attempts=RENAME_ATTEMPTS, delay=RENAME_DELAY_S):
    """os.replace(), retried briefly. On Windows, antivirus software can still
    be scanning a file that was just closed and refuse the rename for a moment
    ([WinError 5] Access denied / [WinError 32]) -- seen with "Download
    SnpEff..." on a PC running Avast."""
    for attempt in range(attempts):
        try:
            os.replace(src, dst)
            return
        except PermissionError as exc:
            if attempt == attempts - 1:
                raise PermissionError(
                    f"{exc}\n\nThe download finished, but Windows wouldn't let the file be renamed "
                    "into place. This is usually antivirus software still scanning it -- try again, "
                    "or add this folder to your antivirus exclusions."
                ) from exc
            time.sleep(delay)


def download_file(url, dest_path, on_progress=None, opener=None, should_cancel=None, direct=False):
    """Streams url to dest_path via a .part temp file, renamed on success.

    direct=True writes straight to dest_path instead (no .part, no rename) --
    for scratch files that are deleted after use anyway, e.g. the Java/SnpEff
    zips. That avoids the rename, which antivirus sandboxing (Avast's
    auto-sandbox runs a new exe virtualized on its first launch) can refuse.

    on_progress(bytes_done, total_bytes_or_None). opener defaults to
    urllib.request.urlopen; tests can inject a fake to avoid a real network call.
    """
    opener = opener or urllib.request.urlopen
    tmp_path = dest_path if direct else dest_path + ".part"
    try:
        with opener(url) as response:
            total = response.getheader("Content-Length") if hasattr(response, "getheader") else None
            total = int(total) if total else None
            done = 0
            with open(tmp_path, "wb") as out:
                while True:
                    if should_cancel and should_cancel():
                        raise DownloadCancelled("Download cancelled.")
                    chunk = response.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if on_progress:
                        on_progress(done, total)
        if not direct:
            _replace_with_retry(tmp_path, dest_path)
    except DownloadCancelled:
        _remove_if_exists(tmp_path)
        raise
    except Exception as exc:
        _remove_if_exists(tmp_path)
        raise DownloadError(f"Download failed: {exc}") from exc


def decompress_gz(gz_path, dest_path, on_progress=None, should_cancel=None):
    """Decompresses gz_path to dest_path. Progress is reported as compressed
    bytes consumed (the true uncompressed size isn't known upfront for
    multi-GB gzip members), which is a fine proxy for a progress bar."""
    tmp_path = dest_path + ".part"
    try:
        total = os.path.getsize(gz_path)
        with gzip.open(gz_path, "rb") as src, open(tmp_path, "wb") as out:
            while True:
                if should_cancel and should_cancel():
                    raise DownloadCancelled("Decompression cancelled.")
                chunk = src.read(CHUNK_SIZE)
                if not chunk:
                    break
                out.write(chunk)
                if on_progress:
                    raw = getattr(src, "fileobj", None)
                    pos = raw.tell() if raw is not None else None
                    on_progress(pos, total)
        _replace_with_retry(tmp_path, dest_path)
    except DownloadCancelled:
        _remove_if_exists(tmp_path)
        raise
    except Exception as exc:
        _remove_if_exists(tmp_path)
        raise DownloadError(f"Decompression failed: {exc}") from exc


def build_fai_index(fasta_path):
    """Builds fasta_path + '.fai' using pyfaidx -- no samtools required."""
    from pyfaidx import Fasta

    try:
        Fasta(fasta_path, rebuild=True)
    except Exception as exc:
        raise DownloadError(f"Indexing failed: {exc}") from exc


def fetch_reference(build, dest_dir, on_progress=None, opener=None, should_cancel=None):
    """Full pipeline: download -> decompress -> index. Returns the final .fa path.

    on_progress(phase, done, total) with phase in {"download", "decompress", "index"};
    total is None for indeterminate progress (e.g. "index" has no meaningful fraction).
    Any failure (including cancellation) removes every file this call created --
    never leaves a half-downloaded .gz or half-built .fa/.fai behind.
    """
    if build not in REFERENCE_URLS:
        raise DownloadError(f"Unknown build: {build!r} (expected one of {list(REFERENCE_URLS)})")

    url = REFERENCE_URLS[build]
    os.makedirs(dest_dir, exist_ok=True)
    gz_path = os.path.join(dest_dir, f"{build}.fa.gz")
    fa_path = os.path.join(dest_dir, f"{build}.fa")
    fai_path = fa_path + ".fai"

    try:
        download_file(
            url, gz_path,
            on_progress=lambda done, total: on_progress and on_progress("download", done, total),
            opener=opener, should_cancel=should_cancel,
        )
        if on_progress:
            on_progress("decompress", 0, None)
        decompress_gz(
            gz_path, fa_path,
            on_progress=lambda done, total: on_progress and on_progress("decompress", done, total),
            should_cancel=should_cancel,
        )
        if on_progress:
            on_progress("index", 0, None)
        build_fai_index(fa_path)
        if on_progress:
            on_progress("index", 1, 1)
    except Exception:
        _remove_if_exists(gz_path)
        _remove_if_exists(fa_path)
        _remove_if_exists(fai_path)
        raise
    else:
        _remove_if_exists(gz_path)

        return fa_path
