import sys

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication, QDialog

from . import theme
from .assets_paths import LOGO_PNG
from .first_launch_dialog import FirstLaunchDialog
from .main_window import MainWindow
from .reference_note import show_reference_note_if_needed
from .spliceai_setup import activate_packages_dir, is_spliceai_installed
from .spliceai_setup_dialog import SpliceAISetupDialog


def configure_app_icon(app):
    """Sets the app-wide default window icon (inherited by any window/dialog
    that doesn't set its own). Split out from main() so it's callable/testable
    without invoking main()'s blocking app.exec()."""
    if LOGO_PNG.exists():
        app.setWindowIcon(QIcon(str(LOGO_PNG)))


def show_first_launch_notice(app):
    """Shows the mandatory startup notice and waits for it to be dismissed.

    Runs before the main window exists, so this is the one point where the
    about-to-close dialog could be Qt's "last window" -- quitOnLastWindowClose
    defaults to True, which would otherwise risk quitting the whole app the
    moment this dialog closes, before MainWindow ever gets a chance to open.
    Disabled for the duration of this call and restored right after, so the
    dialog closing can never be mistaken for the app's last window closing.

    Returns True if the user clicked "I have read and accept the terms",
    False if they clicked "Decline" -- callers must not build MainWindow at
    all when this returns False.
    """
    previous = app.quitOnLastWindowClosed()
    app.setQuitOnLastWindowClosed(False)
    try:
        dialog = FirstLaunchDialog()
        result = dialog.exec()
    finally:
        app.setQuitOnLastWindowClosed(previous)
    return result == QDialog.Accepted


def show_spliceai_setup_notice_if_needed(app):
    """Shows the SpliceAI setup dialog if the package isn't importable --
    same quitOnLastWindowClosed guard as show_first_launch_notice, and for
    the same reason (no MainWindow exists yet at this point).

    Unlike the first-launch notice, this is never a hard gate: SpliceAI is
    only required for live scoring, so both "Done" (installed successfully)
    and "Continue without SpliceAI" (declined/skipped) let startup proceed --
    callers don't need to check this call's return value to decide whether to
    build MainWindow. It exists purely to make the "what's missing and what
    do I do about it" message visible up front instead of only surfacing the
    first time a live-scoring run fails.
    """
    if is_spliceai_installed():
        return
    previous = app.quitOnLastWindowClosed()
    app.setQuitOnLastWindowClosed(False)
    try:
        SpliceAISetupDialog().exec()
    finally:
        app.setQuitOnLastWindowClosed(previous)


def build_and_show_main_window():
    window = MainWindow()
    window.show()
    return window


def main():
    # SpliceAI installed by this app lives in a per-user folder outside the
    # app (see spliceai_setup.PACKAGES_DIR); make it importable before anything
    # checks for it.
    activate_packages_dir()
    app = QApplication(sys.argv)
    theme.apply_theme(app)
    configure_app_icon(app)
    accepted = show_first_launch_notice(app)
    if not accepted:
        return 0
    show_spliceai_setup_notice_if_needed(app)
    window = build_and_show_main_window()
    # Non-blocking: points a fresh install at the reference data it still
    # needs. Skipped when the window found the files by itself.
    show_reference_note_if_needed(window, window.reference_paths())
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
