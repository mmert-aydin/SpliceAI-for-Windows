"""PyInstaller entry point.

Deliberately a thin top-level script rather than pointing PyInstaller
directly at spliceai_gui/__main__.py: PyInstaller executes its entry script
without any package context, so __main__.py's `from .main import main`
(a relative import) would fail with "attempted relative import with no known
parent package" the moment the frozen exe ran. This file lives outside the
spliceai_gui package and uses a plain absolute import instead, sidestepping
that entirely. `python -m spliceai_gui` (via __main__.py) keeps working
unchanged for running from source -- this file is only for the frozen build.

`SpliceAI-VariantScoring.exe --cli <spliceai_pipeline.cli arguments>` runs the
pipeline without the GUI, e.g. to test an installed copy from a script. The
exe is windowed (no console), so write results with -o; progress and any
error go to the file named by the SPLICEAI_CLI_LOG environment variable (if
set), and the exit code is 0 on success, 1 on failure.
"""

import os
import sys
import traceback


def run_cli(argv):
    log_path = os.environ.get("SPLICEAI_CLI_LOG")
    log = open(log_path, "a", encoding="utf-8") if log_path else open(os.devnull, "w")
    sys.stdout = sys.stderr = log
    try:
        from spliceai_gui.spliceai_setup import activate_packages_dir
        from spliceai_pipeline.cli import main as cli_main

        activate_packages_dir()
        cli_main(argv)
        return 0
    except BaseException:
        # An uncaught exception in a windowed exe would pop up an error box.
        traceback.print_exc()
        return 1
    finally:
        log.flush()


if __name__ == "__main__":
    if sys.argv[1:2] == ["--cli"]:
        sys.exit(run_cli(sys.argv[2:]))

    from spliceai_gui.main import main

    sys.exit(main())
