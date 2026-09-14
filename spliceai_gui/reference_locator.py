"""Finds the reference genome FASTAs by itself, so a fresh install doesn't ask
the user to point at files that are already on the machine (or on the USB drive
Setup was run from).

Setup copies the drive's reference-data folder to ~\\SpliceAI_reference_data
when it is run from the USB, so on a normal install location 2 below already
holds both FASTAs and nothing is ever asked. The remaining locations cover the
cases where it doesn't: a user who declined the copy, an older install, files
put somewhere by hand, or the drive still plugged in.

Search order, first usable file per build wins:
  1. %SPLICEAI_REF_DIR%              -- deliberate override, always first
  2. ~\\SpliceAI_reference_data       -- where Setup and "Download..." put them
  3. reference-data next to the app  -- files copied into the install folder
  4. reference-data / SpliceAI-USB-* on each fixed and removable drive

Local disks are searched before removable drives, so a copy on the PC is
preferred over the same file on a USB drive that can be unplugged.

"Usable" means a FASTA whose name identifies the build and that has its .fa.fai
index next to it -- without the index, pyfaidx tries to build one, which takes
minutes and fails outright on a read-only drive. Nothing here downloads,
copies, writes or opens the files; it is only directory listings, so it is
cheap enough to run every time the window opens.
"""
import ctypes
import glob
import os
from pathlib import Path

from spliceai_pipeline.app_paths import app_root_dir

from . import config

BUILDS = ("hg19", "hg38")
FASTA_SUFFIXES = (".fa", ".fasta", ".fna")

# The folder "Download reference..." fills, and the one Setup copies the USB
# drive's reference-data into. Kept in sync with the installer's [Files] entry.
USER_REFERENCE_DIR = Path.home() / "SpliceAI_reference_data"

_DRIVE_REMOVABLE = 2
_DRIVE_FIXED = 3


def _drive_roots():
    """Roots of the fixed and removable drives, local disks first. Network
    drives are left out on purpose: listing one can block for a long time when
    the share is slow or gone, and reading a 3 GB FASTA over a network share
    isn't something to arrange behind the user's back anyway."""
    try:
        mask = ctypes.windll.kernel32.GetLogicalDrives()
        drive_type = ctypes.windll.kernel32.GetDriveTypeW
    except (AttributeError, OSError):  # not Windows -- used by the tests
        return []
    fixed, removable = [], []
    for i in range(26):
        if not mask & (1 << i):
            continue
        root = f"{chr(ord('A') + i)}:\\"
        kind = drive_type(ctypes.c_wchar_p(root))
        if kind == _DRIVE_FIXED:
            fixed.append(root)
        elif kind == _DRIVE_REMOVABLE:
            removable.append(root)
    return fixed + removable


def candidate_dirs():
    """Folders to look in, in the order they are preferred."""
    dirs = []
    override = os.environ.get("SPLICEAI_REF_DIR")
    if override:
        dirs.append(override)
    dirs.append(str(USER_REFERENCE_DIR))
    app_dir = app_root_dir()
    dirs += [os.path.join(app_dir, "reference-data"), os.path.join(app_dir, "reference_data")]
    for root in _drive_roots():
        dirs.append(os.path.join(root, "reference-data"))
        dirs.append(os.path.join(root, "SpliceAI_reference_data"))
        # The USB layout: Setup.exe and reference-data inside a versioned folder.
        dirs += sorted(glob.glob(os.path.join(root, "SpliceAI-USB-*", "reference-data")))
    # Keep the first occurrence of each folder.
    seen, unique = set(), []
    for path in dirs:
        key = os.path.normcase(os.path.normpath(path))
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def is_usable(path):
    """A FASTA this program can read straight away: the file exists and its
    .fai index is next to it."""
    return bool(path) and os.path.isfile(path) and os.path.isfile(path + ".fai")


def _fastas_in(directory):
    """{build: path} for the usable FASTAs in one folder. An exact hg19.fa /
    hg38.fa wins over any other name that merely mentions the build."""
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return {}
    found = {}
    for name in names:
        if not name.lower().endswith(FASTA_SUFFIXES):
            continue
        build = config.guess_build_from_name(name)
        if build not in BUILDS:
            continue
        path = os.path.join(directory, name)
        if not is_usable(path):
            continue
        preferred = name.lower() in {f"{build}{suffix}" for suffix in FASTA_SUFFIXES}
        if build not in found or preferred:
            found[build] = path
    return found


def find_references(dirs=None):
    """{build: path} for every build found, searching candidate_dirs()."""
    found = {}
    for directory in dirs if dirs is not None else candidate_dirs():
        for build, path in _fastas_in(directory).items():
            found.setdefault(build, path)
        if len(found) == len(BUILDS):
            break
    return found


def resolve(settings, dirs=None):
    """The FASTA paths to start the window with, and where they came from.

    Returns ({build: path}, {build: path}). A saved path that still works is
    always kept -- automatic detection only fills a build whose setting is
    empty or points at a file that is no longer there (moved, deleted, or a USB
    drive that isn't plugged in). The second dict holds only the builds that
    were filled in automatically, so the window can say where they came from
    (and stop saying it once the user picks a different file).
    """
    paths, detected = {}, {}
    missing = []
    for build in BUILDS:
        saved = (settings.get(config.fasta_key(build)) or "").strip()
        if saved and os.path.isfile(saved):
            paths[build] = saved
        else:
            paths[build] = ""
            missing.append(build)
    if missing:
        for build, path in find_references(dirs).items():
            if build in missing:
                paths[build] = path
                detected[build] = path
    return paths, detected
