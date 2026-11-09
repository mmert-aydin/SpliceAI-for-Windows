"""Background thread that drives the pipeline so the GUI never blocks.

Reuses spliceai_pipeline.cli.run_pipeline_core() directly -- no pipeline logic
is reimplemented here, just Qt signal plumbing around it.
"""

from PySide6.QtCore import QThread, Signal

from spliceai_pipeline.cli import ReferenceMismatchError, run_pipeline_core
from spliceai_pipeline.control import RunCancelled, RunControl



def _spliceai_missing_message():
    """What to tell a user whose run needed SpliceAI: no Python or command line
    for it (see spliceai_setup.install_spliceai_wheel)."""
    from .i18n import t

    return t("SpliceAI is not installed, so live scoring could not be completed in this run (variants not covered by "
             "precomputed data).\n"
             "\n"
             "To install it, click the red \"SpliceAI not installed\" button at the top of the window, then \"Download "
             "and install\" (an internet connection is needed once) or \"Install from file...\".")


SPLICEAI_MISSING_MESSAGE = _spliceai_missing_message()


class PipelineWorker(QThread):
    progress = Signal(str, object, object)
    finished_ok = Signal(list)
    failed = Signal(str)
    cancelled = Signal()

    def __init__(self, vcf_path, build, mode, fasta_path, precomputed_dir, skip_precomputed,
                 use_snpeff=False, snpeff_dir=None, mane_dir=None, use_mane=True, genes=None, parent=None):
        super().__init__(parent)
        self.vcf_path = vcf_path
        self.build = build
        self.mode = mode
        self.fasta_path = fasta_path
        self.precomputed_dir = precomputed_dir
        self.skip_precomputed = skip_precomputed
        self.use_snpeff = use_snpeff
        self.snpeff_dir = snpeff_dir
        self.mane_dir = mane_dir
        self.use_mane = use_mane
        self.genes = genes
        # The GUI's Pause/End buttons (see spliceai_pipeline.control).
        self.control = RunControl()

    def pause(self):
        self.control.pause()

    def resume(self):
        self.control.resume()

    def stop(self):
        self.control.cancel()

    def run(self):
        try:
            def on_progress(message, current=None, total=None):
                self.progress.emit(message, current, total)

            rows = run_pipeline_core(
                self.vcf_path, self.build, self.mode, self.fasta_path,
                precomputed_dir=self.precomputed_dir or None,
                skip_precomputed=self.skip_precomputed,
                use_snpeff=self.use_snpeff,
                snpeff_dir=self.snpeff_dir or None,
                mane_dir=self.mane_dir or None,
                use_mane=self.use_mane,
                on_progress=on_progress,
                control=self.control,
                genes=self.genes,
            )
            self.finished_ok.emit(rows)
        except RunCancelled:
            self.cancelled.emit()
        except ReferenceMismatchError as exc:
            # Already written for the user (wrong build / wrong FASTA).
            self.failed.emit(str(exc))
        except ModuleNotFoundError as exc:
            if exc.name == "spliceai" or (exc.name or "").startswith("spliceai."):
                self.failed.emit(SPLICEAI_MISSING_MESSAGE)
                return
            self.failed.emit(f"{type(exc).__name__}: {exc}")
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")
