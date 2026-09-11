"""First-run environment check: does this machine actually have everything a
full analysis needs, and -- where it doesn't -- can the app fix it for the
user without them leaving the GUI?

Split deliberately into pure-logic checks here (each returns a CheckResult,
no Qt, no side effects, individually testable against a fake environment)
and the dialog that renders them + wires up the "Fix" buttons
(system_check_dialog.py).

What we check, and why each one is where it is on the required/optional line:

- SpliceAI importable            REQUIRED  -- live scoring is impossible without it.
- Reference FASTA set + indexed  REQUIRED  -- every run reads it; a missing/uindexed
                                             FASTA fails the pipeline immediately.
- Writable working directory     REQUIRED  -- the pipeline writes its output table
                                             and temp files next to the .exe / output path.
- Free disk space (>= 4 GB)      WARNING   -- a fresh reference download is ~1 GB
                                             compressed + ~3 GB decompressed.
- Java >= MIN_JAVA_VERSION       OPTIONAL  -- only SnpEff needs it; only checked when
                                             SnpEff is enabled in settings.
- SnpEff jar + build database    OPTIONAL  -- same: only when the user opted into SnpEff.
- MANE Select summary present    OPTIONAL  -- improves transcript choice; pipeline runs
                                             without it.

`app_root_dir()` (frozen: the folder holding the .exe) is the same base the
pipeline already uses for its default snpeff/ and mane/ folders, so the
"Fix" actions drop downloads exactly where the pipeline will later look.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from enum import Enum

from spliceai_pipeline.app_paths import app_root_dir

# Roughly what a fresh hg19/hg38 download needs on disk (compressed + decompressed
# + .fai), with headroom. Only a soft warning -- the download dialog itself will
# also fail loudly if the disk truly fills up mid-write.
MIN_FREE_BYTES = 4 * 1024**3


class Severity(str, Enum):
    REQUIRED = "required"   # blocks "Continue" until fixed or explicitly overridden
    WARNING = "warning"     # shown, but "Continue" stays enabled
    OPTIONAL = "optional"   # only surfaced when the related feature is enabled


class Status(str, Enum):
    OK = "ok"
    MISSING = "missing"
    ERROR = "error"         # check itself couldn't run (treat like MISSING for gating)


# Machine-readable id for each check, so the dialog can map a row to the right
# "Fix" handler without string-matching the title.
class CheckId(str, Enum):
    SPLICEAI = "spliceai"
    REFERENCE_FASTA = "reference_fasta"
    WORKING_DIR = "working_dir"
    DISK_SPACE = "disk_space"
    JAVA = "java"
    SNPEFF = "snpeff"
    MANE = "mane"


@dataclass
class CheckResult:
    id: CheckId
    title: str
    severity: Severity
    status: Status
    detail: str = ""
    # True when there's a concrete in-app action for this row (Download.../Install...).
    fixable: bool = False
    fix_label: str = "Fix…"
    # Free-form extras a fix handler may need (e.g. {"build": "hg19"}).
    context: dict = field(default_factory=dict)

    @property
    def blocking(self) -> bool:
        return self.severity is Severity.REQUIRED and self.status is not Status.OK


# --------------------------------------------------------------------------- #
# Individual checks. Each takes the merged settings dict (see config.load()).
# --------------------------------------------------------------------------- #

def check_spliceai(_settings: dict) -> CheckResult:
    from spliceai_gui.spliceai_setup import is_spliceai_installed

    if is_spliceai_installed():
        return CheckResult(
            CheckId.SPLICEAI, "SpliceAI model", Severity.REQUIRED, Status.OK,
            "Installed and importable.",
        )
    return CheckResult(
        CheckId.SPLICEAI, "SpliceAI model", Severity.REQUIRED, Status.MISSING,
        "Required for live variant scoring. It is licensed separately by "
        "Illumina and is not bundled with this app.",
        fixable=True, fix_label="Set up SpliceAI…",
    )


def check_reference_fasta(settings: dict) -> CheckResult:
    build = settings.get("build") or "hg19"
    path = (settings.get("fasta_path") or "").strip()
    ctx = {"build": build}

    if not path:
        return CheckResult(
            CheckId.REFERENCE_FASTA, "Reference genome (FASTA)", Severity.REQUIRED,
            Status.MISSING,
            f"No reference FASTA is set. The {build} genome (~3 GB) can be "
            "downloaded and indexed automatically.",
            fixable=True, fix_label="Download reference…", context=ctx,
        )
    if not os.path.isfile(path):
        return CheckResult(
            CheckId.REFERENCE_FASTA, "Reference genome (FASTA)", Severity.REQUIRED,
            Status.MISSING,
            f"Configured FASTA is missing:\n{path}",
            fixable=True, fix_label="Download reference…", context=ctx,
        )
    if not os.path.isfile(path + ".fai"):
        # pyfaidx will build this on first use, but doing it now surfaces a
        # permission/corruption problem before the user starts a real run.
        return CheckResult(
            CheckId.REFERENCE_FASTA, "Reference genome (FASTA)", Severity.WARNING,
            Status.MISSING,
            f"FASTA found but its .fai index is missing:\n{path}\n"
            "It will be built on first use; building it now is safer.",
            fixable=True, fix_label="Build index now", context={"fasta_path": path},
        )
    return CheckResult(
        CheckId.REFERENCE_FASTA, "Reference genome (FASTA)", Severity.REQUIRED,
        Status.OK, f"{path}",
    )


def check_working_dir(_settings: dict) -> CheckResult:
    root = app_root_dir()
    probe = os.path.join(root, ".spliceai_write_test")
    try:
        with open(probe, "w") as fh:
            fh.write("ok")
        os.remove(probe)
    except OSError as exc:
        return CheckResult(
            CheckId.WORKING_DIR, "Writable working folder", Severity.REQUIRED,
            Status.ERROR,
            f"Cannot write into:\n{root}\n({exc})\n\n"
            "Move the app somewhere you own (e.g. your Desktop or Documents), "
            "or run it as administrator.",
        )
    return CheckResult(
        CheckId.WORKING_DIR, "Writable working folder", Severity.REQUIRED, Status.OK,
        root,
    )


def check_disk_space(settings: dict) -> CheckResult:
    # Check the volume the reference would land on (its configured download dir,
    # else the app root).
    target = settings.get("last_reference_download_dir") or app_root_dir()
    # Walk up to an existing directory so shutil.disk_usage doesn't raise.
    probe = target
    while probe and not os.path.isdir(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            break
        probe = parent
    try:
        free = shutil.disk_usage(probe or app_root_dir()).free
    except OSError as exc:
        return CheckResult(
            CheckId.DISK_SPACE, "Free disk space", Severity.WARNING, Status.ERROR,
            f"Could not determine free space ({exc}).",
        )
    gb = free / 1024**3
    if free < MIN_FREE_BYTES:
        return CheckResult(
            CheckId.DISK_SPACE, "Free disk space", Severity.WARNING, Status.MISSING,
            f"Only {gb:.1f} GB free on the target drive. A fresh reference "
            "download needs about 4 GB.",
        )
    return CheckResult(
        CheckId.DISK_SPACE, "Free disk space", Severity.WARNING, Status.OK,
        f"{gb:.0f} GB free",
    )


def check_java(settings: dict) -> CheckResult:
    from spliceai_pipeline.snpeff import MIN_JAVA_VERSION, JavaNotFoundError, find_java

    severity = Severity.OPTIONAL if not settings.get("use_snpeff") else Severity.REQUIRED
    try:
        java_exe = find_java()
    except JavaNotFoundError as exc:
        return CheckResult(
            CheckId.JAVA, f"Java {MIN_JAVA_VERSION}+ (for SnpEff)", severity,
            Status.MISSING, str(exc),
            fixable=True, fix_label="How to install Java…",
        )
    return CheckResult(
        CheckId.JAVA, f"Java {MIN_JAVA_VERSION}+ (for SnpEff)", severity, Status.OK,
        str(java_exe),
    )


def check_snpeff(settings: dict) -> CheckResult:
    from spliceai_pipeline.snpeff import (
        DEFAULT_SNPEFF_DIR,
        SnpEffNotSetUpError,
        check_snpeff_setup,
    )

    build = settings.get("build") or "hg19"
    snpeff_dir = (settings.get("snpeff_dir") or "").strip() or DEFAULT_SNPEFF_DIR
    severity = Severity.OPTIONAL if not settings.get("use_snpeff") else Severity.REQUIRED
    ctx = {"build": build, "snpeff_dir": snpeff_dir}
    try:
        check_snpeff_setup(snpeff_dir, build)
    except SnpEffNotSetUpError as exc:
        return CheckResult(
            CheckId.SNPEFF, "SnpEff annotation", severity, Status.MISSING, str(exc),
            fixable=True, fix_label="Download SnpEff…", context=ctx,
        )
    except Exception as exc:  # noqa: BLE001 -- report, never crash the check
        return CheckResult(
            CheckId.SNPEFF, "SnpEff annotation", severity, Status.ERROR, str(exc),
            context=ctx,
        )
    return CheckResult(
        CheckId.SNPEFF, "SnpEff annotation", severity, Status.OK, snpeff_dir, context=ctx,
    )


def check_mane(_settings: dict) -> CheckResult:
    from spliceai_pipeline.mane import DEFAULT_MANE_DIR, find_cached_summary

    try:
        summary = find_cached_summary(DEFAULT_MANE_DIR)
    except Exception:  # noqa: BLE001
        summary = None
    if summary:
        return CheckResult(
            CheckId.MANE, "MANE Select data", Severity.OPTIONAL, Status.OK,
            str(summary),
        )
    return CheckResult(
        CheckId.MANE, "MANE Select data", Severity.OPTIONAL, Status.MISSING,
        "Not downloaded. Improves transcript selection in SnpEff annotation; "
        "the pipeline still runs without it.",
        fixable=True, fix_label="Download MANE…",
    )


ALL_CHECKS = (
    check_spliceai,
    check_reference_fasta,
    check_working_dir,
    check_disk_space,
    check_java,
    check_snpeff,
    check_mane,
)


def run_checks(settings: dict) -> list[CheckResult]:
    """Run every check against ``settings`` (a config.load() dict). Never raises;
    a check that blows up becomes a Status.ERROR row."""
    results: list[CheckResult] = []
    for fn in ALL_CHECKS:
        try:
            results.append(fn(settings))
        except Exception as exc:  # noqa: BLE001
            results.append(CheckResult(
                CheckId(fn.__name__.replace("check_", "")),
                fn.__name__.replace("check_", "").replace("_", " ").title(),
                Severity.WARNING, Status.ERROR, f"Check failed to run: {exc}",
            ))
    return results


def has_blocking_failure(results: list[CheckResult]) -> bool:
    return any(r.blocking for r in results)


def visible_results(results: list[CheckResult]) -> list[CheckResult]:
    """OPTIONAL rows that are OK are hidden to keep the list short; everything
    else (any REQUIRED/WARNING, or a failing OPTIONAL) is shown."""
    out = []
    for r in results:
        if r.severity is Severity.OPTIONAL and r.status is Status.OK:
            continue
        out.append(r)
    return out
