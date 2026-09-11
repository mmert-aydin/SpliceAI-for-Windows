# Bundled versions — internal build note

Exact inputs of the current Setup.exe, so a rebuild can reproduce it. Update
this file whenever a pin changes. (Not user-facing.)

Built: 2026-09-11 on Windows 10 Pro 22H2 (10.0.19045), Setup version 1.0.0.

## Build tools

| Tool | Version |
|---|---|
| Python (build venv `.venv`) | 3.13.1 (python.org, `%LOCALAPPDATA%\Programs\Python\Python313`) |
| pip | 24.3.1 |
| PyInstaller | 6.22.2 (pyinstaller-hooks-contrib 2026.7) |
| Inno Setup | 6.7.3 (`winget install --id JRSoftware.InnoSetup -e --scope user`) |

## Python packages frozen into the app

Top-level pins are in `requirements.txt`. Full resolved set (`pip freeze` of `.venv`):

```
absl-py==2.5.0
altgraph==0.17.5
astunparse==1.6.3
biopython==1.88
certifi==2026.7.22
charset-normalizer==3.5.1
flatbuffers==25.12.19
gast==0.7.0
google-pasta==0.2.0
grpcio==1.83.1
h5py==3.16.0
idna==3.19
keras==3.15.1
libclang==18.1.1
Markdown==3.10.3
markdown-it-py==4.2.0
MarkupSafe==3.0.3
mdurl==0.1.2
ml_dtypes==0.6.0
namex==0.1.0
numpy==2.5.1
opt_einsum==3.4.0
optree==0.20.0
packaging==26.3
pandas==3.0.5
pefile==2024.8.26
pillow==12.3.0
protobuf==7.36.1
pyfaidx==0.9.0.4
Pygments==2.21.0
pyinstaller==6.22.2
pyinstaller-hooks-contrib==2026.7
PySide6==6.11.1
PySide6_Addons==6.11.1
PySide6_Essentials==6.11.1
python-dateutil==2.9.0.post0
pywin32-ctypes==0.2.3
requests==2.34.2
rich==15.0.0
setuptools==80.9.0
shiboken6==6.11.1
six==1.17.0
tensorboard==2.20.0
tensorboard-data-server==0.7.2
tensorflow==2.20.0
termcolor==3.3.0
typing_extensions==4.16.0
tzdata==2026.3
urllib3==2.7.0
Werkzeug==3.1.8
wheel==0.48.0
wrapt==2.4.1
```

To rebuild with exactly these: `.venv\Scripts\python -m pip install -r <this list saved as a file>`.

## SpliceAI — downloaded on the user's PC, never bundled

The single source of truth is `spliceai_gui/spliceai_setup.py`; `build.ps1`
passes these values into the installer.

| | |
|---|---|
| Package | spliceai 1.3.1, installed without dependencies (all needed ones are frozen into the app; pysam isn't used) |
| File | `spliceai-1.3.1-py2.py3-none-any.whl` (16,676,295 bytes) |
| URL | https://files.pythonhosted.org/packages/d6/2b/9dbf72fdd948cd606c21826cc3735a5beea52633dab72d95d9936a9454d4/spliceai-1.3.1-py2.py3-none-any.whl |
| SHA-256 | `63c633b8de6803d4ffb613e6ef62f1a896fec9fc13245c0f000ff7adea279776` |
| Installed to | `{app}\_internal` by Setup; `~\.spliceai_gui\python-packages` by the app's fallback |

## Third-party payload (`installer/thirdparty`, git-ignored)

| Component | Version | Pin / SHA-256 |
|---|---|---|
| SnpEff | 5.4c (2026-02-23), `snpEff_latest_core.zip` | `snpEff.jar` `5E8F75CBF908A33C6FB2E65C81E66FE31236CB21BB0195541C4703D8202C22B3` |
| SnpEff DB GRCh37.p13 | `snpEff_v5_1_GRCh37.p13.zip` (via `snpEff download`) | `snpEffectPredictor.bin` `A1AD41028A513AEB3D31E3548B327260D15D81981A9D952A1F8FF0694DE6F06C` |
| SnpEff DB GRCh38.p14 | via `snpEff download` (2026-09-10) | `snpEffectPredictor.bin` `43DB6EC1ED69C7BD3B4358FBE69F2BAE5FB392673BA5E056DA743229B08ED632` |
| Java | Eclipse Temurin JRE 21.0.12.1+1 (x64, 2026-08-18) | release `jdk-21.0.12.1+1`, checksum from Adoptium API |
| MANE | `MANE.GRCh38.v1.5.summary.txt.gz` (NCBI release_1.5) | `D10ACE2720681A3B2E0EEFD9DA4F551274A6B4141AC9BFD6A2565DFB6E9AD55C` |
| VC++ runtime | 14.44.35208.0 x64 (from the build PC's VC++ 2015-2022 redistributable) | see below |

SnpEff's `examples/`, `galaxy/` and `.claude/` folders are left out of Setup (not used).

VC++ runtime DLLs (installed next to the exe):

```
concrt140.dll             8BF82E025265D30FB25FDA7161F8E8D73B96F83FDE10F95B0FC4887BD6B09B5E
msvcp140.dll              8F141B4454FA78DB34BC1F28C571B4DA0E00CD2C43F7AD0E282F313036826AAE
msvcp140_1.dll            4D547DB5EB3AABFB1BBE25604C37FE1C253FBA3BF25FAE654D78D7CD44B7DD3D
msvcp140_2.dll            8F481C24983FE695E29DF1BFF55A44072A1C96927322E70844EF5E6A58ECF8D1
msvcp140_atomic_wait.dll  44843C7541B75E30F8FF191BF5D5785A0885F9D694F6E10E545C3346ACD5D0C3
msvcp140_codecvt_ids.dll  2A61FAA8F81E397E04E8AFCD82BE7CA95309B144F148B2BD8E6A1E93BEF04BDA
vccorlib140.dll           6FBC46F26F0ED1D6EA837890E60BB2C5B5CF7E1DF7C66D16267B8A9FD3545F9C
vcomp140.dll              E7E50211906AB1C6226A3133FDA8D8C8E35A23A8B6C6133C9606601604680B85
vcruntime140.dll          0205071C36C17F1EFBD70178C852CB7D49985C484202752B8704B7AC6B184E60
vcruntime140_1.dll        963E45EDD064545962E216C12D68071CED94DC8E11862A18F07F14EB2690A57C
vcruntime140_threads.dll  18186DD0AF0E8CB6B0111ADB585AD70DD555C1C6950ECA7BA3A4900558D6EB37
```

## Output

| | |
|---|---|
| Setup.exe | `SpliceAI-VariantScoring-Setup-1.0.0.exe` |
| Size | 1,412,618,904 bytes (1.35 GB; LZMA2 ultra64, solid) |
| SHA-256 | `5F2BA552EEBC0A526D15A9581F4CB1B340CE923EB70E40400F4B2C052CCDBD7F` |
| Includes | the pipeline fixes and the automatic build selection of 2026-09-11 (see README) |
| Installed size | app 1.30 GB + SnpEff/Java/MANE ~1.3 GB + SpliceAI 25 MB |
