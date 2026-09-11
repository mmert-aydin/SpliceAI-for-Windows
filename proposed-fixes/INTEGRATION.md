# First-run system check — integration notes

## New / changed files

| File | Change |
|---|---|
| `spliceai_gui/preflight.py` | **NEW** — pure-logic checks, no Qt. `run_checks(settings) -> list[CheckResult]`. |
| `spliceai_gui/system_check_dialog.py` | **NEW** — the dialog + "Fix" button wiring. `run_system_check(app) -> bool`, `should_show_on_startup() -> bool`. |
| `spliceai_gui/main.py` | **REPLACED** — drops `show_spliceai_setup_notice_if_needed`; adds the system-check step. |
| `spliceai_gui/config.py` | **1-line** — add `"setup_complete": False` to `DEFAULTS`. |

`config.py` DEFAULTS — add the key (keep the existing 10):

```python
DEFAULTS = {
    "fasta_path": "",
    "precomputed_mode": None,
    "precomputed_dir": "",
    "build": "hg19",
    "mode": "masked",
    "skip_precomputed": False,
    "last_reference_download_dir": "",
    "column_order": [],
    "snpeff_dir": "",
    "use_snpeff": False,
    "setup_complete": False,   # <-- new: has the first-run system check passed once
}
```

## Interfaces this feature depends on (verified against the shipped 1.x build)

- `spliceai_gui.config.load() -> dict`, `save(dict)`
- `spliceai_gui.spliceai_setup.is_spliceai_installed() -> bool`
- `spliceai_gui.spliceai_setup_dialog.SpliceAISetupDialog(parent)`
- `spliceai_gui.download_dialog.ReferenceDownloadDialog(parent, default_build)` — on success `.exec()==Accepted` and `.result_path()` returns the `.fa` path. **Caller must persist it as `fasta_path`** (the dialog only saves `last_reference_download_dir`). The dialog is done in `_fix_reference`.
- `spliceai_gui.snpeff_download_dialog.SnpEffDownloadDialog(parent, build, dest_dir)`
- `spliceai_gui.mane_download_dialog.ManeDownloadDialog(parent)`
- `spliceai_pipeline.app_paths.app_root_dir() -> str` (frozen: folder holding the .exe)
- `spliceai_pipeline.snpeff`: `MIN_JAVA_VERSION` (== 21), `DEFAULT_SNPEFF_DIR`, `JavaNotFoundError`, `SnpEffNotSetUpError`, `find_java() -> str` (raises), `check_snpeff_setup(snpeff_dir, build)` (raises)
- `spliceai_pipeline.mane`: `DEFAULT_MANE_DIR`, `find_cached_summary(mane_dir) -> path | None`

> ⚠ If any of these names/signatures differ in the real source, adjust the
> imports in `preflight.py` / `system_check_dialog.py` — the logic is unchanged.

## Startup flow after this change

```
QApplication
  -> FirstLaunchDialog (licence)            Decline -> exit 0
  -> should_show_on_startup()?
       first launch (setup_complete False)  -> SystemCheckDialog
       later launch, a REQUIRED check fails  -> SystemCheckDialog   Quit -> exit 0
       otherwise                             -> skip
  -> MainWindow
```

## Checks

| id | severity | fix action |
|---|---|---|
| `spliceai` | REQUIRED | `SpliceAISetupDialog` |
| `reference_fasta` | REQUIRED (WARNING if only `.fai` missing) | `ReferenceDownloadDialog` + persist `fasta_path` |
| `working_dir` | REQUIRED | (none — advises moving the app / running as admin) |
| `disk_space` | WARNING | (none) |
| `java` | OPTIONAL (REQUIRED if `use_snpeff`) | `SnpEffDownloadDialog` — now downloads a private Java automatically when none is found (`spliceai_gui/java_download.py`) |
| `snpeff` | OPTIONAL (REQUIRED if `use_snpeff`) | `SnpEffDownloadDialog` |
| `mane` | OPTIONAL | `ManeDownloadDialog` |

OPTIONAL rows that pass are hidden. "Continue" is disabled while any REQUIRED
check fails.

## Follow-up (not in this draft)

- **Java auto-install**: DONE (2026-09-10) -- `spliceai_gui/java_download.py`
  pulls a Temurin 21 JRE zip from the Adoptium API into `<app>/java`, and
  `snpeff.find_java()` checks there first; `SnpEffDownloadDialog` triggers it.
  `preflight.py`'s java check can simply offer that dialog.
- **tests**: mirror the existing `tests/` layout — `test_preflight.py` with a
  fake settings dict + monkeypatched `find_java` / `is_spliceai_installed`.

## Build

Environment is already set up at
`…/scratchpad/build-workspace/.venv` (Python 3.13.1):
PySide6 6.11.1, tensorflow 2.20.0, keras 3.15.1, numpy 2.5.1, pyfaidx 0.9.0.4,
spliceai 1.3.1 (`--no-deps`), setuptools 80.9.0 (for `pkg_resources`),
pandas 3.0.5, pyinstaller 6.22.2.

```
cd <source project>
…/build-workspace/.venv/Scripts/pip install -e .        # or -r requirements.txt
…/build-workspace/.venv/Scripts/pyinstaller SpliceAI-VariantScoring.spec --noconfirm
```

Then smoke-test `dist/SpliceAI-VariantScoring/SpliceAI-VariantScoring.exe` from
**outside** the source tree (frozen-path bugs only show there).

> The shipped build has pandas pinned lower than 3.0.5 possibly; match the
> source's `requirements.txt` when it's available.
