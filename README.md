# SpliceAI Variant Scoring — recovered source

This folder is the source code of `SpliceAI-VariantScoring.exe`. It was rebuilt
from the shipped program (a PyInstaller build, Python 3.13) on 2026-09-10,
because the original source was no longer available.

## How faithful is it?

It isn't a decompiler guess. Every module was compiled with Python 3.13 and
checked against the bytecode inside the shipped exe. It matches exactly:
instructions, constants, strings and docstrings, names, exception tables,
and the line where each function and class starts.

| What | Result |
|---|---|
| `spliceai_gui` (21 modules) | identical |
| `spliceai_pipeline` (13 modules) | identical |
| `run_app.py` (entry script) | identical |

**What could not be recovered.** Python throws comments away when it
compiles, so the original comments are gone.

- Where the original had comment lines, you'll find blank lines, a short
  factual comment, or this note:
  *"(Original comments weren't recoverable from the compiled program; blank
  and comment lines like these keep line numbers where they were.)"*
- These keep every statement on its original line, so tracebacks from old
  logs still point at the right place.
- Docstrings are original; they were stored in the program.

**Files that are not in the exe.**

- `spliceai_gui/__main__.py` was recreated from the description in
  `run_app.py`.
- The `.spec`, the requirements files, `.gitignore` and this README are new.

## Layout

```
run_app.py                 entry point for the frozen exe
spliceai_gui/              PySide6 desktop app (main window, dialogs, downloads)
spliceai_pipeline/         scoring pipeline (VCF -> normalize -> lookup/live score -> merge)
LICENSE-THIRD-PARTY.md     third-party licenses (shown in the app)
SpliceAI-VariantScoring.spec
requirements.txt / requirements-build.txt
tools/verify_source.py     re-runs the byte-for-byte comparison
proposed-fixes/            fixes found during testing, NOT applied (see below)
```

## Run from source (Windows, Python 3.13)

```
py -3.13 -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pip install --no-deps spliceai==1.3.1
.venv\Scripts\python -m spliceai_gui
```

**SpliceAI is installed separately on purpose.** Its license forbids
redistributing it, so the app never bundles it. The app checks for it at
startup and offers to install it.

- `--no-deps` skips pysam, which has no Windows build and isn't needed.
- The copy of the app that was tested on this PC never showed the
  "SpliceAI missing" screen. That's because spliceai had already been
  installed into its `_internal` folder.

Command line (no GUI):

```
.venv\Scripts\python -m spliceai_pipeline.cli input.vcf --build hg38 --mode raw --fasta hg38.fa
```

## Build the exe

```
.venv\Scripts\python -m pip install -r requirements-build.txt
.venv\Scripts\python -m PyInstaller SpliceAI-VariantScoring.spec --noconfirm
```

Result: `dist\SpliceAI-VariantScoring\SpliceAI-VariantScoring.exe`, a windowed
onedir build like the original.

- SnpEff and MANE data are downloaded next to the exe at runtime.
- The reference FASTA can live anywhere.

## Build the Windows installer (Setup.exe)

This is what gets handed out, e.g. on a USB drive. Setup installs the app per
user (no admin rights) into `%LOCALAPPDATA%\Programs\SpliceAI-VariantScoring`,
with SnpEff (GRCh37 + GRCh38 databases), Java, MANE and the VC++ runtime
included, plus Start Menu and Desktop shortcuts and an uninstaller.

SpliceAI itself is never in Setup.exe. During Setup it's downloaded from PyPI
(the exact file pinned in `spliceai_gui/spliceai_setup.py`, SHA-256-checked)
and unpacked into the app. If a PC is offline, Setup says so and finishes;
the app then offers the download itself. A copy of the same wheel next to
Setup.exe is used instead of downloading.

```
winget install --id JRSoftware.InnoSetup -e --scope user   # once
powershell -ExecutionPolicy Bypass -File installer\fetch-thirdparty.ps1   # fresh checkout only
powershell -ExecutionPolicy Bypass -File installer\build.ps1 -Version 1.0.0
```

Result: `installer\output\SpliceAI-VariantScoring-Setup-1.0.0.exe`.

- `installer\SpliceAI-VariantScoring.iss` is the Inno Setup script; `build.ps1`
  builds the app with PyInstaller into `installer\build` and passes it the
  SpliceAI wheel URL/hash.
- `installer\BUNDLED-VERSIONS.md` lists the exact versions in the current build.
- `installer\test\run-sandbox.ps1` tests Setup.exe in Windows Sandbox (a clean
  Windows, optionally without network): install, scoring, relaunch, uninstall.
- The installed exe also runs headless:
  `SpliceAI-VariantScoring.exe --cli input.vcf --build hg19 --mode masked --fasta hg19.fa -o out.tsv`.

## Re-checking the source against an exe

1. Extract the exe with [pyinstxtractor-ng](https://github.com/pyinstxtractor/pyinstxtractor-ng).
2. Run the comparison:

   ```
   python tools/verify_source.py <exe>_extracted/PYZ.pyz_extracted . spliceai_gui/main_window spliceai_pipeline/cli ...
   ```

Each module prints `EXACT` (identical), `NOPS` (identical except for line-only
padding instructions), or `DIFF` with the first difference.

## Changes since recovery

The recovered, byte-identical version of every file is what `tools/verify_source.py`
checks. These files have since been changed on purpose, so the tool reports
them as `DIFF`:

**"Download SnpEff…" no longer needs Java installed by hand (2026-09-10).**

- **Before:** SnpEff needs Java 21+. The app used to stop with "No Java 21+
  found … install a JDK/JRE". That happened on PCs with no Java at all, or
  with only an old Java 8.
- **Now:** when no Java 21+ is found, "Download SnpEff…" first downloads the
  official Eclipse Temurin 21 runtime from Adoptium.
  - It is checksum-verified and unpacked into a `java` folder next to the app.
  - Only this app uses it. Nothing is installed system-wide, and PATH isn't
    touched.
  - `find_java()` looks in that folder first.
  - Parts that are already installed (SnpEff, the database) aren't downloaded
    again.
  - When only Java is missing, the SnpEff status line says so.

Files touched:

- `spliceai_gui/java_download.py` (new)
- `spliceai_gui/snpeff_download.py`
- `spliceai_gui/snpeff_download_dialog.py`
- `spliceai_pipeline/snpeff.py`
- `spliceai_gui/main_window.py` (status check and SnpEff help text)

**Two more fixes on the same button.**

- `spliceai_gui/reference_download.py`: the final rename of a finished
  download is retried for up to 10 seconds. Antivirus scanning the new file
  used to make it fail with "Access denied".
- `spliceai_pipeline/snpeff.py`: every Java launch gets
  `-Duser.language=en -Duser.country=US`. On Turkish-locale Windows, SnpEff
  otherwise misreads its own `String` type name ("STRİNG").
- No console window flashes when Java starts.
- The Java and SnpEff zips are downloaded into Windows' temp folder and written
  there directly (`download_file(..., direct=True)`), with no `.part` file and
  no rename. This works around Avast's auto-sandbox: on the first launch of a
  newly built exe, Avast runs the app virtualized, and the sandbox refused
  exactly that rename (confirmed in Avast's own log). Only the unpacked files
  end up next to the app.

**SpliceAI stays installed across launches and app updates (2026-09-10).**

- **Before:** "Install automatically" (and the manual command) put SpliceAI
  inside the app's own `_internal` folder. That folder is replaced whenever
  the app is updated, re-copied or rebuilt, so the install silently
  disappeared and the app said "SpliceAI isn't installed" again.
- **Now:** SpliceAI is installed per user into `~\.spliceai_gui\python-packages`,
  next to the settings file.
  - `spliceai_setup.activate_packages_dir()` puts that folder on `sys.path` at
    startup (called first thing in `main.main()`).
  - It also registers the folder with `pkg_resources`. PyInstaller's runtime
    hook builds pkg_resources' package list before the app's code runs, and
    SpliceAI's `__init__` calls `get_distribution("spliceai")`.
- Older installs inside `_internal` are still found.
- Files: `spliceai_gui/spliceai_setup.py`, `spliceai_gui/main.py`.

**Installer and offline setup (2026-09-11).** See "Build the Windows installer" above.

- `spliceai_gui/spliceai_setup.py`, `spliceai_setup_dialog.py`: the installed app
  installs SpliceAI without Python or pip (downloads the pinned wheel, checks its
  SHA-256, unzips it); "Install from file..." takes a copy from a USB drive.
- `spliceai_gui/reference_note.py` (new), `main.py`, `config.py`: a non-blocking note
  after startup while no reference FASTA is set, so an empty field isn't mistaken
  for a broken install ("Don't show this again" is remembered).
- `run_app.py`: `--cli` runs the pipeline headless from the installed exe.

**Pipeline fixes from a code review (2026-09-11).**

- **A wrong genome build no longer fails silently.** Every REF allele is checked
  against the FASTA before scoring; if more than 10% don't match, the run stops with
  a plain explanation and nothing is scored (a wrong build gives ~75%). Pressing Run
  also warns when the VCF's `##reference` header names the other build.
- **Variants that can't be scored say why** in the source column (e.g. "skipped:
  REF doesn't match reference FASTA", "skipped: REF allele too long for SpliceAI")
  instead of a blank "none" row; "none" now only means "no gene overlaps".
  SpliceAI's own skip warnings used to be hidden. The Summary dialog counts them.
- **A chromosome the FASTA doesn't have no longer crashes the run** (it's listed
  as skipped); MT and chrM count as the same; a FASTA whose names match none of the
  VCF's (e.g. RefSeq `NC_...`) stops with a clear message. One variant that fails
  to score no longer ends a whole run.
- **Blank lines in a VCF** (e.g. pasted text) no longer crash the run.
- **Complex indels** (e.g. AT>GCA) are no longer shifted onto a different variant
  during left-alignment.
- **SnpEff's transcript/HGVS is matched to each row's gene**: a variant overlapping
  two genes no longer shows one gene's transcript on the other gene's row.
- **The run-time "SpliceAI isn't installed" error** in the installed app points to
  the in-app install button instead of a pip command.

Files: `spliceai_pipeline/normalize.py`, `vcfio.py`, `score.py`, `snpeff.py`,
`merge.py`, `cli.py`; `spliceai_gui/worker.py`, `main_window.py`, `vcf_check.py`.

**Genome build is selected from the VCF; one FASTA per build (2026-09-11).**

- Loading or pasting a VCF selects its build automatically, from the header:
  `##reference` / `##assembly` (hg19, GRCh37, b37, hs37d5, hg38, GRCh38, ...), or
  else chr1's length in the `##contig` lines. A label next to the Build dropdown
  says where it came from, that the file doesn't say, or -- if the build was
  changed by hand -- that the VCF names the other one. A hand-picked build isn't
  overridden while the same VCF is edited.
- The single "Reference FASTA" field became two rows, **hg19 reference FASTA** and
  **hg38 reference FASTA**, each with Browse... and Download...; a run uses the
  selected build's row, which is marked "(used for this run)". Browsing to a file
  whose name names the other build asks first.
- Settings: `fasta_path_hg19` / `fasta_path_hg38`; an older single `fasta_path` is
  moved to the matching one on first start (by file name, else the saved build).
- Files: `spliceai_gui/vcf_info.py` (`detect_vcf_build`), `main_window.py`,
  `config.py`, `download_dialog.py`, `reference_note.py`.

## Still open

**First-run check of the PC's prerequisites.** The reference genome, SnpEff
and so on aren't checked up front. A draft is in `proposed-fixes/`
(`preflight.py`, `system_check_dialog.py`, notes in `INTEGRATION.md`); it
isn't applied.
