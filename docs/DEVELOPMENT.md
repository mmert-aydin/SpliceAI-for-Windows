# Development notes

How to run from source, build the exe and the installer, and what changed in
each version. For using the program, see the [README](../README.md).

## History: the source was recovered from the compiled program

The source of version 1.0 was rebuilt from the shipped program (a PyInstaller
build, Python 3.13) on 2026-09-10, because the original source was no longer
available on the build PC. Everything since 1.0 was developed on top of it.

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

- Where the original had comment lines, you'll find blank lines or a short
  factual comment written during the recovery.
- Docstrings are original; they were stored in the program.
- `spliceai_gui/__main__.py` was recreated from the description in
  `run_app.py`; the `.spec` and the requirements files are new.
- The byte-for-byte comparison tool (`tools/verify_source.py`) is in the git
  history, commit `f0177ff`.

## Layout

```
run_app.py                 entry point for the frozen exe
spliceai_gui/              PySide6 desktop app (main window, dialogs, downloads)
spliceai_pipeline/         scoring pipeline (VCF -> normalize -> lookup/live score -> merge)
installer/                 Inno Setup script, build and fetch scripts, Windows Sandbox tests
tests/                     plain-script checks (tests\run-tests.ps1)
examples/                  a small demo VCF of public ClinVar variants
LICENSE-THIRD-PARTY.md     third-party licenses (shown in the app)
SpliceAI-VariantScoring.spec
requirements.txt / requirements-build.txt
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

The reference genomes aren't in Setup.exe either -- 6 GB, and they're public
reference files. When Setup is started from a drive that has a `reference-data`
folder next to Setup.exe, it offers to copy them to
`%USERPROFILE%\SpliceAI_reference_data`; the program looks there first, so
after that it never asks for them (`spliceai_gui/reference_locator.py`).

SpliceAI itself is never in Setup.exe. During Setup it's downloaded from PyPI
(the exact file pinned in `spliceai_gui/spliceai_setup.py`, SHA-256-checked)
and unpacked into the app. If a PC is offline, Setup says so and finishes;
the app then offers the download itself. A copy of the same wheel next to
Setup.exe is used instead of downloading.

```
winget install --id JRSoftware.InnoSetup -e --scope user   # once
powershell -ExecutionPolicy Bypass -File installer\fetch-thirdparty.ps1   # fresh checkout only
powershell -ExecutionPolicy Bypass -File installer\build.ps1 -Version 1.1.1
```

Result: `installer\output\SpliceAI-VariantScoring-Setup-1.1.1.exe`.

- `installer\SpliceAI-VariantScoring.iss` is the Inno Setup script; `build.ps1`
  builds the app with PyInstaller into `installer\build` and passes it the
  SpliceAI wheel URL/hash.
- `installer\BUNDLED-VERSIONS.md` lists the exact versions in the current build.
- `installer\test\run-sandbox.ps1` tests Setup.exe in Windows Sandbox (a clean
  Windows, optionally without network): install, scoring, relaunch, uninstall.
  Scenarios: `online`, `offline`, `offline-usb-wheel`, `usb-reference` (the
  genomes are copied off the drive and found by the program), `clickthrough`
  (drives the real wizard and window, with screenshots), `manual`.
- The installed exe also runs headless:
  `SpliceAI-VariantScoring.exe --cli input.vcf --build hg19 --mode masked --fasta hg19.fa -o out.tsv`.

## Changes since recovery

Each entry is dated. The 2026-09-10/11 entries make up 1.0; the 2026-09-14
entries are what 1.1.0 added; 1.1.1 is the last entry.

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

**Run controls, MANE folder, SpliceAI check, scrolling (2026-09-11).**

- **Pause / Resume and End** next to Run, enabled only during a run. The pipeline
  checks in between variants (`spliceai_pipeline/control.py`), so a pause takes
  effect after the current variant; End stops at the next one -- or at once while
  SnpEff runs (its Java process is killed; SnpEff now runs via a polled Popen
  instead of `subprocess.run`). An ended run keeps no results.
- **Closing during a run** asks "The calculation is not finished yet. Are you sure
  you want to exit?"; Yes stops the run cleanly before the window closes.
- **CPU display** during a run: this program's share of the processor, the whole
  PC's, and the elapsed time (`spliceai_gui/cpu_meter.py`, standard library only).
- **MANE Select folder** field with Browse... (checks it's a MANE summary file);
  remembered, passed to the run, and used as the download destination. The status
  says whether the file is found there.
- **SpliceAI button always shown**: green "SpliceAI found · Check" opens a check of
  the install (version, location, where it came from, all 5 model and 2 annotation
  files present); red when missing, opening the install dialog.
- **The whole window scrolls**, so small screens work; the window also starts no
  bigger than the screen, and the results table keeps a usable minimum height.

**One screen you can read, not a form to fill in (2026-09-14).** Everything
that comes with Setup stopped asking to be configured, and the results table
started showing which variants matter.

- **Order.** Results first, with Run / Pause / End under them, then the VCF,
  then Settings -- and the three are panes of a splitter, so the line between
  any two can be dragged. Where they are left is remembered (`pane_sizes`). The
  VCF box shrank to a few lines: a loaded VCF is read from its file, not from
  what is on show there.
- **Advanced settings.** Settings now holds only the two real per-run choices,
  raw/masked and precomputed data. The genome build, the reference FASTAs,
  SnpEff and MANE moved behind one button, with a line saying whether
  everything is in place. It opens by itself when something isn't. The build
  stays visible read-only -- scoring against the wrong one is wrong quietly.
- **SnpEff is on by default** (it ships with Setup), applied once to an
  existing settings file. Turning it off now says what that costs, MANE
  included, instead of leaving a gap where its rows were.
- **A Region column**, from SnpEff's own `ANN` effect and rank: `Intron 21/33`,
  `Exon 23/23`, `5' UTR`, `Upstream`, `Intergenic`. Always from the same
  transcript the Transcript column names. Also written to the results file.
- **Columns reordered** to gene, transcript, variant, region, then the max
  score, then the eight scores it came from, then its source. A column order
  saved by an older version is dropped once, otherwise the new order would
  never be seen.
- **The max score is drawn, not just printed**: a bar filled in proportion and
  a stripe down the row, coloured at SpliceAI's own 0.20 / 0.50 / 0.80 cutoffs
  -- the same three values the filter buttons use, now read from one place.
  Presentation only; sorting, filtering and export still see the number.
- **One palette** (`spliceai_gui/theme.py`) instead of a dozen colours written
  at the widgets that used them. One teal accent; green, amber and red kept for
  meaning. Windows' own fonts, nothing bundled. Scrollbars and splitter handles
  restyled to match.

Files: `theme.py`, `score_delegate.py` (both new), `main_window.py`,
`table_model.py`, `config.py`, `main.py`, `first_launch_dialog.py` (a contact
address), `spliceai_pipeline/snpeff.py`, `mane.py`, `merge.py`, `writer.py`.

**Java and MANE are found beside SnpEff (2026-09-14).** Setup puts `snpeff`,
`java` and `mane` next to each other, but `find_java()` and the MANE lookup
only ever searched next to the *app*. Running the program from anywhere else --
from source, or after the install folder is moved -- reported "No Java 21+
found" with a perfectly good Temurin 21 sitting beside the jar it had just
found, and SnpEff then switched itself off. Both now also look one level up
from the chosen SnpEff folder. Files: `spliceai_pipeline/snpeff.py`
(`sibling_java_dir`), `mane.py` (`sibling_mane_dir`), `main_window.py`.

**The reference genome finds itself (2026-09-14).** A new install asked for a
file that was, on the USB drive, right next to Setup.exe -- the first thing a
colleague saw was two empty "reference FASTA" rows and a note about browsing or
downloading 3 GB. Two changes, and between them nobody should have to answer
that question again:

- Setup, when it is run from a drive that carries a `reference-data` folder,
  offers to **copy the genomes to `~\SpliceAI_reference_data`** (about 6 GB,
  ticked by default; it checks there is room, and doesn't copy again over an
  existing install). Declining still installs everything else.
- The program **looks for them itself** every time it opens: the
  `SPLICEAI_REF_DIR` override, `~\SpliceAI_reference_data`, a `reference-data`
  folder next to the exe, then `reference-data` / `SpliceAI-USB-*` on each fixed
  and removable drive -- local disks before removable ones, so a copy on the PC
  wins over the same file on a drive that can be unplugged. A file only counts
  if its name says which build it is (a wrong build gives wrong scores, quietly)
  and its `.fa.fai` index is next to it. A path already saved in the settings is
  never overridden -- only an empty one, or one whose file has gone.

The two FASTA rows are now **folded away behind "Reference files..."**, with one
line in their place saying which builds are ready and where they were found.
They start open only when a reference still has to be found, and the
first-launch reference note is skipped entirely when the files were located.
Files: `spliceai_gui/reference_locator.py` (new), `main_window.py`,
`reference_note.py`, `main.py`, `installer/SpliceAI-VariantScoring.iss`; checks
in `tests/test_reference_locator.py`.

**Downloads work on networks that inspect HTTPS (2026-09-14).** On a PC whose
network (or antivirus) re-signs secure connections -- common in hospitals --
every download in the app failed with an error Windows itself accepts:
`CERTIFICATE_VERIFY_FAILED: Missing Authority Key Identifier`. Python checks
certificates with its own bundled list; the Windows installer had always used
Windows' own (hence SpliceAI downloading fine there while "Download
reference..." failed). All downloads -- reference FASTA, SnpEff, Java, MANE and
the app's own SpliceAI install -- now verify through Windows via `truststore`
(the approach pip takes), and a certificate failure explains what to do
(copy the files from the USB drive, or ask IT) instead of showing only the raw
error. Files: `spliceai_pipeline/net.py` (new), `mane.py`,
`spliceai_gui/reference_download.py`, `java_download.py`; `truststore==0.10.4`
in requirements.txt; checks in `tests/test_downloads.py`.

**Drag and drop anywhere on the window (2026-09-11).** Before, only the VCF box
took dropped files; elsewhere nothing happened, and a file dropped onto a path
field was pasted in as a `file:///` address. Now a file dragged from Explorer
anywhere onto the window is used by its name: a VCF (`.vcf`, `.vcf.gz`,
`.vcf.bgz`, `.bgz`) loads like "Load file..."; a FASTA (`.fa`, `.fasta`, `.fna`)
goes into its build's reference row (the row it's dropped on, else the build its
name says, else the selected one -- with the same "different build?" question as
Browse...); a MANE summary file sets the MANE folder; a folder dropped onto a
folder field sets it; anything else gets a short explanation. The VCF box shows its
dashed frame while a VCF is dragged over the window. Files:
`spliceai_gui/main_window.py` (`PathLineEdit`, `MainWindow.dropEvent`); checks in
`tests/test_drag_drop.py`.

**Files saved "with BOM" are read normally (2026-09-11).** Some Windows editors
(and Windows PowerShell) write UTF-8 with an invisible byte-order mark. Before,
such a `config.json` or gene-list file was silently ignored (all settings back to
defaults), and such a VCF stopped the run (its first header line was taken for a
broken record). Settings, gene lists and VCFs are now read as `utf-8-sig`; VCFs
also replace undecodable header characters instead of failing. Files:
`spliceai_gui/config.py`, `gene_list_storage.py`, `main_window.py`,
`spliceai_pipeline/vcfio.py`; checks in `tests/`.

**1.1.1: rebuilt for the public release (2026-09-30).** No change in what
the program does. A few docstrings and comments used the names of real VCF
files from development as examples; they now use made-up names. Python keeps
docstrings inside the compiled program, so Setup.exe was rebuilt to leave the
old names out of it too. Built on a second PC with Python 3.13.2 (1.1.0 used
3.13.1); every other pinned version and bundled file is the same, see
`installer/BUNDLED-VERSIONS.md`.

## Still open

**First-run check of the PC's prerequisites.** The reference genome, SnpEff
and so on aren't checked up front. A draft (`preflight.py`,
`system_check_dialog.py`, notes in `INTEGRATION.md`) is in the git history,
in the `proposed-fixes/` folder of commit `f0177ff`; it isn't applied.
