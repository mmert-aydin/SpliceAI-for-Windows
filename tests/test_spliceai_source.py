"""From source the program finds SpliceAI as the installed one does (spec 17):
in ~/.spliceai_gui/python-packages, with the same setup dialog if it's missing.
A separate process with an empty home, and nothing but the repository on its path.

    .venv\\Scripts\\python tests\\test_spliceai_source.py
"""
import json
import os
import subprocess
import sys
import tempfile

import helpers
from helpers import check

CHILD = r"""
import importlib.util, json, os, sys
os.environ["QT_QPA_PLATFORM"] = "offscreen"
out = {"in_environment": importlib.util.find_spec("spliceai") is not None}
from PySide6.QtWidgets import QApplication, QPushButton
app = QApplication(sys.argv)
from spliceai_gui import spliceai_setup as s, worker
from spliceai_gui.spliceai_setup_dialog import EXPLANATION_HTML, SpliceAISetupDialog
s.activate_packages_dir()
out["on_path"] = str(s.PACKAGES_DIR) in sys.path
out["before"] = s.is_spliceai_installed()
dialog = SpliceAISetupDialog()
out["buttons"] = [b.text() for b in dialog.findChildren(QPushButton)]
out["texts"] = EXPLANATION_HTML + worker.SPLICEAI_MISSING_MESSAGE
s.install_spliceai_wheel(wheel_path=sys.argv[1])
out["after"] = s.refresh_after_install()
info = s.check_installation()
out["info"] = {k: info[k] for k in ("version", "origin", "missing_models", "missing_annotations")}
import spliceai
out["imported_from"] = os.path.dirname(os.path.dirname(spliceai.__file__))
out["packages_dir"] = str(s.PACKAGES_DIR)
print("RESULT " + json.dumps(out), flush=True)
"""

helpers.watchdog(300)
helpers.setup()
from spliceai_gui.spliceai_setup import SPLICEAI_WHEEL_NAME  # noqa: E402

home = tempfile.mkdtemp(prefix="spliceai_source_")
env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
env.update(USERPROFILE=home, HOME=home, PYTHONPATH=str(helpers.REPO))
run = subprocess.run([sys.executable, "-W", "ignore", "-c", CHILD, str(helpers.CACHE / SPLICEAI_WHEEL_NAME)], env=env,
                     cwd=home, capture_output=True, text=True, timeout=240)
line = next((x for x in run.stdout.splitlines() if x.startswith("RESULT ")), None)
if line is None:
    check("the child process ran", False, (run.stdout + run.stderr)[-600:])
    helpers.finish()
got = json.loads(line[7:])
check(".venv doesn't hold SpliceAI (installer/build.ps1 refuses it)", not got["in_environment"])
check("from source python-packages is on the path, and SpliceAI found missing there",
      got["on_path"] and not got["before"])
check("the setup dialog is the installed program's: 'Download and install', 'Install from file...', no pip",
      {"Download and install", "Install from file..."} <= set(got["buttons"]) and "pip " not in got["texts"]
      and "Download and install" in got["texts"], str(got["buttons"]))
check("after installing the pinned package it is found, imported from python-packages, complete",
      got["after"] and os.path.samefile(got["imported_from"], got["packages_dir"]) and got["info"]["version"] == "1.3.1"
      and got["info"]["origin"] == "installed from within the program" and not got["info"]["missing_models"]
      and not got["info"]["missing_annotations"], str(got["info"]))
build = (helpers.REPO / "installer" / "build.ps1").read_text(encoding="utf-8")
check("build.ps1 says what to remove and why", "pip uninstall -y spliceai" in build and "forbids redistributing" in build)
helpers.finish()
