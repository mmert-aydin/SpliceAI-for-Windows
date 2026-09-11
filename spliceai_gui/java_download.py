"""Fetch a private Java runtime for SnpEff, so nobody has to install Java by hand.

SnpEff is a Java program and needs Java 21+. Many PCs have no Java at all, or
only an old Java 8 on PATH. When spliceai_pipeline.snpeff.find_java() finds
nothing suitable, "Download SnpEff..." calls fetch_java() first: it downloads
the official Eclipse Temurin JRE zip from Adoptium, checks it against the
SHA-256 Adoptium publishes, and unzips it into DEFAULT_JAVA_DIR (a "java"
folder next to the app, alongside "snpeff" and "mane"). find_java() looks
there first. Nothing is installed system-wide and PATH is never changed.
"""
import hashlib
import json
import os
import platform
import shutil
import tempfile
import urllib.request
import zipfile

from spliceai_pipeline.snpeff import DEFAULT_JAVA_DIR, MIN_JAVA_VERSION, java_major_version, private_java_exe

from .reference_download import DownloadCancelled, DownloadError, download_file

ADOPTIUM_ASSETS_URL = (
    "https://api.adoptium.net/v3/assets/latest/{version}/hotspot"
    "?architecture={arch}&image_type=jre&os=windows&vendor=eclipse"
)

# Adoptium's API answers 403 to urllib's default "Python-urllib/x.y" agent.
USER_AGENT = "SpliceAI-VariantScoring (Java runtime for SnpEff)"


def _request(url):
    return urllib.request.Request(url, headers={"User-Agent": USER_AGENT})


def _arch():
    return "aarch64" if platform.machine().lower() in ("arm64", "aarch64") else "x64"


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_jre_package(opener=None):
    """Asks Adoptium for the latest Temurin JRE (MIN_JAVA_VERSION, Windows) and
    returns its package info: {"name", "link", "size", "checksum", ...}."""
    opener = opener or urllib.request.urlopen
    url = ADOPTIUM_ASSETS_URL.format(version=MIN_JAVA_VERSION, arch=_arch())
    try:
        with opener(_request(url)) as response:
            assets = json.load(response)
    except Exception as exc:
        raise DownloadError(f"Could not look up the Java download at {url}: {exc}") from exc
    for asset in assets:
        package = (asset.get("binary") or {}).get("package") or {}
        if package.get("name", "").endswith(".zip") and package.get("link"):
            return package
    raise DownloadError(f"Adoptium listed no Windows Java {MIN_JAVA_VERSION} zip to download ({url}).")


def fetch_java(dest_dir=None, on_progress=None, opener=None, should_cancel=None):
    """Downloads, verifies and unzips a Temurin JRE into dest_dir (default
    DEFAULT_JAVA_DIR). Returns the path of its java.exe.

    on_progress(phase, done, total) with phase "java" (download bytes) or
    "java_unzip" (indeterminate).

    The zip itself goes to a scratch folder under the system temp directory,
    never into dest_dir: the app folder is often on the Desktop, where antivirus
    ransomware protection can stop an unfamiliar program from creating a .zip
    (the final .part -> .zip rename then fails with "Access denied", however
    long it's retried). Only the unpacked runtime ends up in dest_dir, and the
    scratch folder is always removed.
    """
    dest_dir = dest_dir or DEFAULT_JAVA_DIR
    package = find_jre_package(opener)
    os.makedirs(dest_dir, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix="spliceai_java_")
    zip_path = os.path.join(work_dir, package["name"])

    try:
        download_file(
            _request(package["link"]), zip_path,
            on_progress=lambda done, total: on_progress and on_progress("java", done, total or package.get("size")),
            opener=opener, should_cancel=should_cancel, direct=True,
        )
        expected = (package.get("checksum") or "").lower()
        if expected and _sha256(zip_path) != expected:
            raise DownloadError("The downloaded Java package is damaged (checksum mismatch). Please try again.")

        if should_cancel and should_cancel():
            raise DownloadCancelled("Cancelled before unpacking Java.")
        if on_progress:
            on_progress("java_unzip", 0, None)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(dest_dir)

        java_exe = private_java_exe(dest_dir)
        if not java_exe:
            raise DownloadError(f"Java was unpacked into {dest_dir}, but no bin\\java.exe was found in it.")
        version = java_major_version(java_exe)
        if version is None or version < MIN_JAVA_VERSION:
            raise DownloadError(f"The downloaded Java ({java_exe}) reports version {version}, not {MIN_JAVA_VERSION}+.")
        if on_progress:
            on_progress("java_unzip", 1, 1)
    except (DownloadCancelled, DownloadError):
        raise
    except Exception as exc:
        raise DownloadError(f"Java setup failed: {exc}") from exc
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return java_exe
