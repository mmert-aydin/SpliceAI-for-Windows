"""GUI entry point.

Startup sequence:
  1. app icon
  2. FirstLaunchDialog          -- licence/terms notice, every launch (Decline -> exit)
  3. SystemCheckDialog          -- first launch, or any launch where a REQUIRED
                                   prerequisite (SpliceAI, reference FASTA, a
                                   writable working folder) is currently missing.
                                   Offers in-app download/install for what's
                                   missing. (Quit -> exit)
  4. MainWindow

The old standalone "SpliceAI is not installed" notice
(show_spliceai_setup_notice_if_needed) is subsumed by step 3 -- SpliceAI is now
just one row in the system check, alongside the reference genome, Java/SnpEff,
MANE data and disk space.
"""

import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication, QDialog

from spliceai_gui.assets_paths import LOGO_PNG
from spliceai_gui.first_launch_dialog import FirstLaunchDialog
from spliceai_gui.main_window import MainWindow
from spliceai_gui.system_check_dialog import run_system_check, should_show_on_startup


def configure_app_icon(app):
    """Set the app-wide default window icon (inherited by any window/dialog
    that doesn't set its own). Split out from main() so it's callable/testable
    without invoking main()'s blocking app.exec()."""
    if LOGO_PNG.exists():
        app.setWindowIcon(QIcon(str(LOGO_PNG)))


def show_first_launch_notice(app):
    """Modal licence/terms notice. Returns True if accepted, False if declined.
    Keeps the app alive (no main window yet) while it's up."""
    previous = app.quitOnLastWindowClosed()
    app.setQuitOnLastWindowClosed(False)
    try:
        dialog = FirstLaunchDialog()
        result = dialog.exec()
    finally:
        app.setQuitOnLastWindowClosed(previous)
    return result == QDialog.Accepted


def build_and_show_main_window():
    window = MainWindow()
    window.show()
    return window


def main():
    app = QApplication(sys.argv)
    configure_app_icon(app)

    if not show_first_launch_notice(app):
        return 0

    if should_show_on_startup() and not run_system_check(app):
        return 0

    window = build_and_show_main_window()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
