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


class _MemoryCounters(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + [
        (name, ctypes.c_size_t) for name in (
            "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage",
            "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]


class ProcessMeter:
    """CPU and RAM of another process (the database build): sample() gives
    (cpu_percent of the whole processor, working set in bytes) since the
    previous sample, both None off Windows or once the process is gone."""

    def __init__(self, pid):
        self.pid = pid
        self._cpus = os.cpu_count() or 1
        self._handle = None
        if sys.platform == "win32":
            # PROCESS_QUERY_LIMITED_INFORMATION | PROCESS_VM_READ
            kernel32 = ctypes.windll.kernel32
            kernel32.OpenProcess.restype = ctypes.c_void_p
            self._handle = kernel32.OpenProcess(0x1000 | 0x0010, False, pid) or None
        self._last = self._now()

    def _now(self):
        if not self._handle:
            return None
        times = [ctypes.c_ulonglong() for _ in range(4)]  # creation, exit, kernel, user (100 ns)
        if not ctypes.windll.kernel32.GetProcessTimes(ctypes.c_void_p(self._handle), *map(ctypes.byref, times)):
            return None
        return time.monotonic(), (times[2].value + times[3].value) / 1e7

    def sample(self):
        now = self._now()
        if now is None or self._last is None:
            return None, None
        (wall0, cpu0), (wall1, cpu1) = self._last, now
        if wall1 - wall0 < 0.2:
            return None, None
        self._last = now
        cpu = max(0.0, min(100.0, (cpu1 - cpu0) / ((wall1 - wall0) * self._cpus) * 100))
        counters = _MemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        ok = ctypes.windll.kernel32.K32GetProcessMemoryInfo(ctypes.c_void_p(self._handle), ctypes.byref(counters),
                                                             counters.cb)
        return cpu, (counters.WorkingSetSize if ok else None)

    def close(self):
        if self._handle:
            ctypes.windll.kernel32.CloseHandle(ctypes.c_void_p(self._handle))
            self._handle = None


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
