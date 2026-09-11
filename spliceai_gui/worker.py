"""Background thread that drives the pipeline so the GUI never blocks.

Reuses spliceai_pipeline.cli.run_pipeline_core() directly -- no pipeline logic
is reimplemented here, just Qt signal plumbing around it.
"""

from PySide6.QtCore import QThread, Signal

from spliceai_pipeline.cli import ReferenceMismatchError, run_pipeline_core

from .spliceai_setup import display_install_command, is_frozen


def _spliceai_missing_message():
    """What to tell a user whose run needed SpliceAI. The installed app needs
    no Python or command line for it (see spliceai_setup.install_spliceai_wheel),
    so only a source checkout is shown a pip command."""
    if is_frozen():
        how = (
            'Click the red "SpliceAI not installed" button at the top of the window, then '
            '"Download and install" (needs an internet connection once) or "Install from file...".'
        )
    else:
        how = (
            f"Install it with:\n{display_install_command()}\n\n"
            "See the SpliceAI status button in the top bar for details and an automatic-install option."
        )
    return (
        "SpliceAI isn't installed, so this run couldn't complete live scoring "
        "(any variant not already covered by precomputed data).\n\n" + how
    )


SPLICEAI_MISSING_MESSAGE = _spliceai_missing_message()


class PipelineWorker(QThread):
    progress = Signal(str, object, object)
    finished_ok = Signal(list)
    failed = Signal(str)

    def __init__(self, vcf_path, build, mode, fasta_path, precomputed_dir, skip_precomputed,
                 use_snpeff=False, snpeff_dir=None, parent=None):
        super().__init__(parent)
        self.vcf_path = vcf_path
        self.build = build
        self.mode = mode
        self.fasta_path = fasta_path
        self.precomputed_dir = precomputed_dir
        self.skip_precomputed = skip_precomputed
        self.use_snpeff = use_snpeff
        self.snpeff_dir = snpeff_dir

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
                on_progress=on_progress,
            )
            self.finished_ok.emit(rows)
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
