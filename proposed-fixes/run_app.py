"""PyInstaller entry point.

Deliberately a thin top-level script rather than pointing PyInstaller
directly at spliceai_gui/__main__.py: PyInstaller executes its entry script
without any package context, so __main__.py's `from .main import main`
(a relative import) would fail with "attempted relative import with no known
parent package" the moment the frozen exe ran. This file lives outside the
spliceai_gui package and uses a plain absolute import instead, sidestepping
that entirely. `python -m spliceai_gui` (via __main__.py) keeps working
unchanged for running from source -- this file is only for the frozen build.
"""

import sys

from spliceai_gui.main import main

if __name__ == "__main__":
    sys.exit(main())
