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

`SpliceAI-VariantScoring.exe --build-refdb <spliceai_pipeline.refdb_build
arguments>` builds the reference score database the same way, as a process of
its own: TensorFlow's thread limit is set before it starts, and a crash can't
take the window down.

`SpliceAI-VariantScoring.exe --check-refdb <folder> <hg19|hg38>`: whether a
reference score database there is usable (see check_refdb); Setup's question.
"""

import os
import sys
import traceback


def check_refdb(folder, build):
    """Setup asks this before copying a database from the USB drive: exit 0 if
    `folder` holds a database the program may use for `build`, 2 if there is
    none, 1 if there is one it refuses (the reason goes to SPLICEAI_CLI_LOG).
    Checked against the pinned model files' hashes: Setup may run this before
    SpliceAI is downloaded, and the program checks the installed ones itself."""
    log_path = os.environ.get("SPLICEAI_CLI_LOG")
    log = open(log_path, "a", encoding="utf-8") if log_path else open(os.devnull, "w")
    try:
        from spliceai_pipeline import refdb
        from spliceai_pipeline.predict import MODEL_SHA256
        return 0 if refdb.check(folder, build, MODEL_SHA256) is not None else 2
    except BaseException as exc:
        log.write(f"{type(exc).__name__}: {exc}\n")
        return 1
    finally:
        log.close()


def run_cli(argv, build_refdb=False):
    log_path = os.environ.get("SPLICEAI_CLI_LOG")
    log = open(log_path, "a", encoding="utf-8") if log_path else open(os.devnull, "w")
    sys.stdout = sys.stderr = log
    try:
        from spliceai_gui.spliceai_setup import activate_packages_dir
        if build_refdb:
            from spliceai_pipeline.refdb_build import main as cli_main
        else:
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
    if sys.argv[1:2] == ["--check-refdb"] and len(sys.argv) == 4:
        sys.exit(check_refdb(sys.argv[2], sys.argv[3]))
    if sys.argv[1:2] == ["--build-refdb"]:
        sys.exit(run_cli(sys.argv[2:], build_refdb=True))

    from spliceai_gui.main import main

    sys.exit(main())
