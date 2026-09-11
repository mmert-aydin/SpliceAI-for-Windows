"""Shared "where does user-downloaded data live" resolution for snpeff.py and
mane.py's DEFAULT_*_DIR constants.

Both SnpEff and MANE default to a plain folder colocated with "the app" --
when running from source, that's the actual project root (two directories
above this file). When frozen with PyInstaller (onedir), naively reusing the
same os.path.dirname(os.path.dirname(__file__)) trick resolves to the
_internal bundle directory instead (this module's __file__ then lives at
e.g. dist/SpliceAI/_internal/spliceai_pipeline/app_paths.py, so two dirname()
calls land on .../_internal, not on dist/SpliceAI) -- a real bug caught by
testing the packaged exe from outside the source tree: downloads would have
silently landed inside the frozen bundle's internals rather than next to the
.exe where a user would expect to find/manage them.

sys.frozen / sys.executable are the standard PyInstaller-set attributes for
detecting this at runtime (see PyInstaller's own docs on runtime information)
-- checked here instead of just trying "does _internal exist" or similar
guesswork.
"""

import os
import sys


def app_root_dir():
    """Directory to treat as "the app's own folder" for default download
    locations (SnpEff, MANE). Frozen: the folder containing the .exe (its
    sibling _internal/ holds the bundled Python runtime and packages, but
    user-downloaded data should sit next to the exe, not inside it, so it
    survives being obviously visible/manageable and isn't mistaken for part
    of the frozen bundle). Not frozen: the actual project root.
    """
    if getattr(sys, "frozen", False):
        return os.path.dirname(os.path.abspath(sys.executable))
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
