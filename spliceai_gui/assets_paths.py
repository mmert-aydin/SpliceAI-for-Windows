"""Resolves paths to bundled static assets (logo etc), independent of the
current working directory the app happens to be launched from.

Works unmodified whether running from source or frozen (PyInstaller), with
no explicit frozen/not-frozen branch needed: __file__-relative resolution is
correct here specifically because these are bundled READ-ONLY resources
travelling *inside* the app's own package structure (spliceai_gui/assets/),
not user-writable external data -- PyInstaller's `datas` preserves the same
source-tree-relative layout inside the frozen bundle, so the same relative
walk up from this file lands in the right place either way. Contrast with
spliceai_pipeline/app_paths.py's DEFAULT_SNPEFF_DIR/DEFAULT_MANE_DIR, which
need to land *next to* the frozen exe (not inside its bundle) because that
data is user-downloaded, not shipped -- a __file__-relative path there
resolves into the wrong place when frozen. Don't copy this file's pattern
for anything the app writes to.
"""
from pathlib import Path

ASSETS_DIR = Path(__file__).resolve().parent / "assets"
LOGO_PNG = ASSETS_DIR / "logo.png"
LOGO_ICO = ASSETS_DIR / "logo.ico"

LICENSE_THIRD_PARTY_PATH = Path(__file__).resolve().parent.parent / "LICENSE-THIRD-PARTY.md"
