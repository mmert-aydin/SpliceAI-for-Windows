"""Background threads of the Sequence-to-cDNA tab, so the window never
blocks: reading a build's gene table, "Check compatibility" (reading and
aligning a paste of several Mb) and "Run". Like worker.PipelineWorker, only
Qt plumbing around spliceai_pipeline.cdna_analysis; Run's Cancel goes
through spliceai_pipeline.control.RunControl, checked between model calls.
"""
from PySide6.QtCore import QThread, Signal

from spliceai_pipeline.control import RunCancelled, RunControl


class TaskWorker(QThread):
    """Runs fn() once: finished_ok(result) or failed(message)."""

    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, parent=None):
        super().__init__(parent)
        self.fn = fn

    def run(self):
        try:
            result = self.fn()
        except Exception as exc:  # noqa: BLE001 -- shown to the user
            message = str(exc) if isinstance(exc, (ValueError, KeyError, FileNotFoundError)) else \
                f"{type(exc).__name__}: {exc}"
            self.failed.emit(message.strip("'\""))
            return
        self.finished_ok.emit(result)


class RunWorker(QThread):
    """cdna_analysis.run_analysis() with progress, pause and cancel."""

    progress = Signal(float, str)
    finished_ok = Signal(object)
    failed = Signal(str)
    cancelled = Signal()

    def __init__(self, run, parent=None):
        """run(checkpoint, on_progress) -> Analysis."""
        super().__init__(parent)
        self._run = run
        self.control = RunControl()

    def stop(self):
        self.control.cancel()

    def run(self):
        try:
            result = self._run(self.control.check, lambda fraction, text: self.progress.emit(fraction, text))
        except RunCancelled:
            self.cancelled.emit()
            return
        except ModuleNotFoundError as exc:
            if exc.name == "spliceai" or (exc.name or "").startswith("spliceai."):
                from .i18n import t

                self.failed.emit(t("SpliceAI is not installed, so the sample's sequence cannot be predicted. To install it, click the SpliceAI "
                                   "button at the top of the window."))
                return
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        except Exception as exc:  # noqa: BLE001 -- shown to the user
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        self.finished_ok.emit(result)
