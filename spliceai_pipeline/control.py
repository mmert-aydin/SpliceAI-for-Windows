"""Pause / end for a running pipeline -- the GUI's Pause and End buttons.

run_pipeline_core() calls RunControl.check() between variants (and while
SnpEff runs), so a pause takes effect after the current variant, and an end
stops the run promptly: at the next variant, or at once while SnpEff is
running (its Java process is killed).
"""
import threading


class RunCancelled(Exception):
    """The run was ended on request -- not an error."""


class RunControl:
    def __init__(self):
        self._running = threading.Event()
        self._running.set()
        self._cancelled = threading.Event()

    def pause(self):
        self._running.clear()

    def resume(self):
        self._running.set()

    def cancel(self):
        self._cancelled.set()
        self._running.set()  # a paused run has to wake up to stop

    @property
    def paused(self):
        return not self._running.is_set()

    @property
    def cancelled(self):
        return self._cancelled.is_set()

    def check(self):
        """Blocks while paused; raises RunCancelled once the run was ended."""
        self._running.wait()
        if self._cancelled.is_set():
            raise RunCancelled("Run ended by the user.")
