"""CPU use while a run is going, shown next to the progress bar: this
program's share of the whole processor (like Task Manager's per-app figure)
and the whole computer's. Standard library only -- os.times() for this
process, GetSystemTimes for the whole PC on Windows.
"""
import ctypes
import os
import sys
import time


def _system_times():
    """(idle, total) CPU time of the whole PC in 100 ns units, or None."""
    if sys.platform != "win32":
        return None
    idle, kernel, user = ctypes.c_ulonglong(), ctypes.c_ulonglong(), ctypes.c_ulonglong()
    if not ctypes.windll.kernel32.GetSystemTimes(ctypes.byref(idle), ctypes.byref(kernel), ctypes.byref(user)):
        return None
    return idle.value, kernel.value + user.value  # kernel time includes idle time


class CpuMeter:
    def __init__(self):
        self._cpus = os.cpu_count() or 1
        self.reset()

    def _now(self):
        t = os.times()
        return time.monotonic(), t.user + t.system, _system_times()

    def reset(self):
        self._last = self._now()

    def sample(self):
        """(app_percent, pc_percent) since the previous sample. pc_percent is
        None off Windows; both are None if less than 0.2 s have passed."""
        wall0, proc0, sys0 = self._last
        now = self._now()
        wall1, proc1, sys1 = now
        if wall1 - wall0 < 0.2:
            return None, None
        self._last = now
        app = max(0.0, min(100.0, (proc1 - proc0) / ((wall1 - wall0) * self._cpus) * 100))
        pc = None
        if sys0 and sys1 and sys1[1] > sys0[1]:
            busy = (sys1[1] - sys0[1]) - (sys1[0] - sys0[0])
            pc = max(0.0, min(100.0, busy / (sys1[1] - sys0[1]) * 100))
        return app, pc
