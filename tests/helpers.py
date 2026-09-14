"""Shared setup for the checks in tests/.

Plain scripts, no test framework: each prints PASS / FAIL / SKIP lines and
exits non-zero if anything failed. See tests/README.md.

Needs the build venv (.venv), the hg19/hg38 FASTAs with their .fai files
(SPLICEAI_REF_DIR, default ~/SpliceAI_reference_data), and SpliceAI -- the
pinned wheel is downloaded once into tests/.cache and unpacked there, never
installed into .venv (installer/build.ps1 refuses a venv that has it).
SnpEff checks use installer/thirdparty and are skipped without it.

Test data is generated from the reference genome here; no patient data is
ever needed or committed.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "tests" / ".cache"
REAL_HOME = Path(os.environ.get("USERPROFILE") or os.environ.get("HOME") or Path.home())
REF_DIR = Path(os.environ.get("SPLICEAI_REF_DIR") or REAL_HOME / "SpliceAI_reference_data")
HG19, HG38 = REF_DIR / "hg19.fa", REF_DIR / "hg38.fa"
THIRDPARTY = REPO / "installer" / "thirdparty"
SNPEFF_DIR = THIRDPARTY / "snpeff" if (THIRDPARTY / "snpeff" / "snpEff" / "snpEff.jar").exists() else None
MANE_DIR = THIRDPARTY / "mane" if (THIRDPARTY / "mane").is_dir() else None
SPLICEAI_PKG = CACHE / "spliceai-pkg"

VCF_COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

sys.path.insert(0, str(REPO))

# Check names quote what the window says, which contains ✓ and ⚠. A console on
# a non-UTF-8 codepage (cp1254 on a Turkish Windows, cp437 in some pipes) would
# otherwise kill the run with UnicodeEncodeError instead of printing a result.
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except (AttributeError, OSError):
    pass

_failures = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not ok:
        _failures.append(name)


def skip(name, why):
    print(f"SKIP {name}  [{why}]", flush=True)


def finish():
    print("\nALL PASSED" if not _failures else f"\n{len(_failures)} FAILED: {_failures}", flush=True)
    sys.exit(1 if _failures else 0)


def watchdog(seconds):
    """Ends this script with a FAIL if it's still running after `seconds` --
    e.g. an unexpected modal dialog in the invisible (offscreen) window would
    otherwise wait for a click forever."""
    import threading

    def fire():
        print(f"FAIL watchdog: still running after {seconds} s -- something is waiting (a dialog?)", flush=True)
        os._exit(3)

    timer = threading.Timer(seconds, fire)
    timer.daemon = True
    timer.start()


def isolate_home():
    """Gives this process a throw-away home folder, so the GUI's settings
    (~/.spliceai_gui) are never the user's real ones. Call it before
    importing anything from spliceai_gui (config.py reads the home folder
    once, at import)."""
    home = tempfile.mkdtemp(prefix="spliceai_test_home_")
    os.environ["USERPROFILE"] = os.environ["HOME"] = home
    return Path(home)


# Run in a separate process, so importing spliceai_gui here can't fix its
# settings path before a test has isolated the home folder.
_PREPARE_SPLICEAI = r"""
import os, sys
from spliceai_gui import spliceai_setup as s
from spliceai_gui.reference_download import download_file
target, cache = sys.argv[1], sys.argv[2]
wheel = os.path.join(cache, s.SPLICEAI_WHEEL_NAME)
if not os.path.exists(wheel):
    download_file(s.SPLICEAI_WHEEL_URL, wheel, direct=True)
s.unpack_wheel(wheel, target)
"""


def setup(need_spliceai=True):
    """Checks the reference files, makes SpliceAI importable (downloading the
    pinned wheel the first time), and puts the bundled Java on PATH."""
    missing = [str(p) for p in (HG19, HG38, Path(f"{HG19}.fai"), Path(f"{HG38}.fai")) if not p.exists()]
    if missing:
        sys.exit(f"Missing reference files {missing} -- set SPLICEAI_REF_DIR to their folder.")
    if need_spliceai:
        if not (SPLICEAI_PKG / "spliceai" / "utils.py").exists():
            CACHE.mkdir(parents=True, exist_ok=True)
            subprocess.run([sys.executable, "-c", _PREPARE_SPLICEAI, str(SPLICEAI_PKG), str(CACHE)],
                           cwd=REPO, check=True)
        sys.path.insert(0, str(SPLICEAI_PKG))
    java_root = THIRDPARTY / "java"
    if java_root.is_dir():
        for runtime in sorted(java_root.iterdir(), reverse=True):
            if (runtime / "bin" / "java.exe").exists():
                os.environ["PATH"] = str(runtime / "bin") + os.pathsep + os.environ["PATH"]
                break


def snv_lines(fasta, chrom, positions):
    """VCF data lines for SNVs at `positions`, REF read from `fasta` (a
    pyfaidx Fasta) and ALT a different base; N positions are left out."""
    lines = []
    for pos in positions:
        ref = fasta[chrom][pos - 1:pos].seq.upper()
        if ref not in ("A", "C", "G", "T"):
            continue
        alt = next(b for b in "ACGT" if b != ref)
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    return lines


def vcf_text(data_lines, header_lines=()):
    return "##fileformat=VCFv4.2\n" + "".join(h + "\n" for h in header_lines) + VCF_COLUMNS + "".join(data_lines)


def write_vcf(data_lines, header_lines=()):
    fd, path = tempfile.mkstemp(suffix=".vcf", prefix="spliceai_test_")
    with os.fdopen(fd, "w") as fh:
        fh.write(vcf_text(data_lines, header_lines))
    return path


def ankrd26_snvs(fasta_hg19, n):
    """n SNVs inside ANKRD26 (hg19), starting at the known-value variant's
    position -- in a gene, so SpliceAI really scores each one."""
    return snv_lines(fasta_hg19, "chr10", [27326999 + 37 * k for k in range(n + 5)])[:n]
