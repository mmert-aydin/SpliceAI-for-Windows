"""Detects whether the `spliceai` package is importable, and -- since it is
never bundled with this app (see LICENSE-THIRD-PARTY.md: SpliceAI's source is
PolyForm Strict, which forbids redistribution, and its trained models are
CC BY-NC; every user has to obtain it directly from Illumina under their own
acceptance of those terms) -- installs it for the user.

The program finds SpliceAI the same way whether it is installed or running
from source: in PACKAGES_DIR, a per-user folder next to the settings file
(~/.spliceai_gui/python-packages), which activate_packages_dir() puts on
sys.path at startup. It is never pip-installed into a Python environment: a
copy in the build's .venv would be swept into the build, which
installer/build.ps1 refuses.

PACKAGES_DIR is deliberately NOT the app's own bundle directory
(sys._MEIPASS, the "_internal" folder): that folder is replaced whenever the
app is updated, re-copied or rebuilt, which silently deleted earlier installs
-- the app then said "SpliceAI isn't installed" again on the next launch. An
install in PACKAGES_DIR survives all of that. Older installs inside _internal
are still found, since _internal is always on sys.path.

Installed builds (installer/, 2026-09-11): Setup downloads the pinned wheel
below from PyPI and unpacks it into _internal, so this module's install path
is only a fallback (e.g. the PC was offline during Setup). It needs no Python
or pip on the PC: a pure-Python wheel with no scripts is just a zip file, so
install_spliceai_wheel() downloads it with the app's own Python, checks its
SHA-256 and unzips it into PACKAGES_DIR -- exactly what pip would do for it.
urllib uses the Windows certificate store and proxy settings, so it works
behind corporate proxies where pip's own certificate bundle fails.
"""

import hashlib
import importlib
import importlib.util
import os
import shutil
import sys
import tempfile
import zipfile

from .config import CONFIG_PATH
from .i18n import N_

SPLICEAI_REPO_URL = "https://github.com/Illumina/SpliceAI"

# The exact SpliceAI release this app is tested with, as published on PyPI.
# installer/build.ps1 reads these same values, so Setup and the app's own
# fallback always fetch the identical, hash-checked file.
SPLICEAI_VERSION = "1.3.1"
SPLICEAI_WHEEL_NAME = f"spliceai-{SPLICEAI_VERSION}-py2.py3-none-any.whl"
SPLICEAI_WHEEL_URL = (
    "https://files.pythonhosted.org/packages/d6/2b/"
    "9dbf72fdd948cd606c21826cc3735a5beea52633dab72d95d9936a9454d4/" + SPLICEAI_WHEEL_NAME
)
SPLICEAI_WHEEL_SHA256 = "63c633b8de6803d4ffb613e6ef62f1a896fec9fc13245c0f000ff7adea279776"

PACKAGES_DIR = CONFIG_PATH.parent / "python-packages"


def is_frozen():
    return bool(getattr(sys, "frozen", False))


def bundle_dir():
    """This frozen app's own package directory (already on sys.path), or
    None when not frozen -- see module docstring."""
    return getattr(sys, "_MEIPASS", None) if is_frozen() else None


def _register_with_pkg_resources(path):
    """spliceai/__init__.py calls pkg_resources.get_distribution("spliceai"),
    which only looks at pkg_resources' own list of installed packages. In the
    frozen app that list is built at startup by PyInstaller's pkg_resources
    runtime hook, before any of this app's code runs -- so a folder added to
    sys.path afterwards (or a package installed into it later) has to be
    added to that list explicitly, or importing spliceai fails with
    DistributionNotFound even though the files are there."""
    pkg_resources = sys.modules.get("pkg_resources")
    if pkg_resources is None:
        return
    try:
        pkg_resources.working_set.add_entry(path)
    except Exception:
        pass


def activate_packages_dir():
    """Makes PACKAGES_DIR importable (call once at startup, before any
    is_spliceai_installed() check): creates it if needed, appends it to
    sys.path, and registers it with pkg_resources. Safe to call again."""
    path = str(PACKAGES_DIR)
    try:
        PACKAGES_DIR.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    if path not in sys.path:
        sys.path.append(path)
        _register_with_pkg_resources(path)


def is_spliceai_installed():
    """Cheap, side-effect-free presence check.

    Deliberately uses importlib.util.find_spec rather than actually
    importing spliceai: spliceai/__init__.py calls signal.signal() at import
    time (see spliceai_pipeline/score.py's module docstring for why that's
    only safe from the main thread), which a mere "is it there" check has no
    business triggering, especially since this is called from startup code
    that runs before anything thread-related is set up.

    Returns False (never raises) for every "not usable" shape: a plain
    missing package (ImportError/ModuleNotFoundError, itself an ImportError
    subclass) or a malformed/partial install with a bad __spec__ (ValueError).
    """
    try:
        return importlib.util.find_spec("spliceai") is not None
    except (ImportError, ValueError):
        return False


def check_installation():
    """Details for the "SpliceAI found" button, without importing SpliceAI
    (see is_spliceai_installed): None if it isn't installed, else a dict with
    its version, folder, where it came from, and any missing model or gene
    annotation files -- a partial install can't live-score."""
    try:
        spec = importlib.util.find_spec("spliceai")
    except (ImportError, ValueError):
        spec = None
    if spec is None or not spec.origin:
        return None
    pkg_dir = os.path.dirname(spec.origin)
    parent = os.path.dirname(pkg_dir)
    version = None
    try:
        for name in os.listdir(parent):
            if name.lower().startswith("spliceai-") and name.endswith(".dist-info"):
                version = name[len("spliceai-"):-len(".dist-info")]
    except OSError:
        pass
    models = [f"spliceai{i}.h5" for i in range(1, 6)]
    annotations = ["grch37.txt", "grch38.txt"]
    here = os.path.normcase(os.path.abspath(pkg_dir))
    if bundle_dir() and here.startswith(os.path.normcase(os.path.abspath(bundle_dir()))):
        origin = N_("installed by Setup, with the program")
    elif here.startswith(os.path.normcase(os.path.abspath(PACKAGES_DIR))):
        origin = N_("installed from within the program")
    else:
        origin = N_("installed in this Python environment")
    return {
        "version": version,
        "location": pkg_dir,
        "origin": origin,
        "missing_models": [m for m in models if not os.path.isfile(os.path.join(pkg_dir, "models", m))],
        "missing_annotations": [a for a in annotations
                                if not os.path.isfile(os.path.join(pkg_dir, "annotations", a))],
        "n_models": len(models),
        "n_annotations": len(annotations),
    }


class WheelInstallError(Exception):
    pass


def sha256_of_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def unpack_wheel(wheel_path, target_dir):
    """Installs the pinned SpliceAI wheel into target_dir the way pip does for
    a pure-Python wheel without scripts or data files: by unzipping it.
    Refuses any file that isn't the exact published one."""
    digest = sha256_of_file(wheel_path)
    if digest != SPLICEAI_WHEEL_SHA256:
        raise WheelInstallError(
            f"{os.path.basename(wheel_path)} is not the published SpliceAI {SPLICEAI_VERSION} "
            f"package (SHA-256 {digest}, expected {SPLICEAI_WHEEL_SHA256})."
        )
    root = os.path.realpath(target_dir)
    os.makedirs(root, exist_ok=True)
    with zipfile.ZipFile(wheel_path) as zf:
        for name in zf.namelist():
            if not os.path.realpath(os.path.join(root, name)).startswith(root + os.sep):
                raise WheelInstallError(f"Unexpected path in wheel: {name}")
        zf.extractall(root)


def install_spliceai_wheel(on_line=None, wheel_path=None):
    """The install: downloads the pinned wheel from PyPI (or uses
    wheel_path, e.g. a copy the user picked from a USB drive), verifies it
    and unpacks it into PACKAGES_DIR. Raises on any failure."""
    from .reference_download import download_file

    log = on_line or (lambda _line: None)
    target = str(PACKAGES_DIR)
    if wheel_path:
        log(f"Installing from {wheel_path}")
        unpack_wheel(wheel_path, target)
    else:
        work_dir = tempfile.mkdtemp(prefix="spliceai_wheel_")
        try:
            path = os.path.join(work_dir, SPLICEAI_WHEEL_NAME)
            log(f"Downloading {SPLICEAI_WHEEL_URL}")
            download_file(SPLICEAI_WHEEL_URL, path, direct=True)
            log("Download finished, checking the file...")
            unpack_wheel(path, target)
        finally:
            shutil.rmtree(work_dir, ignore_errors=True)
    log(f"SpliceAI {SPLICEAI_VERSION} installed into {target}")


def refresh_after_install():
    """Re-checks importability after an install -- invalidate_caches() is
    needed because Python's import system caches "this module doesn't exist"
    lookups (via its finders' negative caches) from before the install, and
    pkg_resources has to re-scan PACKAGES_DIR for the newly installed
    package's metadata (see _register_with_pkg_resources)."""
    activate_packages_dir()
    importlib.invalidate_caches()
    _register_with_pkg_resources(str(PACKAGES_DIR))
    return is_spliceai_installed()
