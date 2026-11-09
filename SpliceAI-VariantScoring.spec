# -*- mode: python ; coding: utf-8 -*-
# Build the Windows app (onedir):
#   python -m PyInstaller SpliceAI-VariantScoring.spec --noconfirm
# Output: dist/SpliceAI-VariantScoring/SpliceAI-VariantScoring.exe
#
# Reconstructed to match the shipped build: windowed exe, no version resource,
# logo.ico icon, assets + LICENSE-THIRD-PARTY.md as data, and SpliceAI itself
# NOT bundled (license) -- the app offers to install it on first run.

from PyInstaller.utils.hooks import collect_dynamic_libs

a = Analysis(
    ["run_app.py"],
    pathex=[],
    # parasail loads parasail.dll from its own folder with ctypes, which
    # PyInstaller doesn't see by itself.
    binaries=collect_dynamic_libs("parasail"),
    datas=[
        # assets_paths.py finds these relative to its own file, so keep the
        # source-tree layout inside the bundle.
        ("spliceai_gui/assets/logo.png", "spliceai_gui/assets"),
        ("spliceai_gui/assets/logo.ico", "spliceai_gui/assets"),
        ("LICENSE-THIRD-PARTY.md", "."),
        # The curated RefSeq tables the Sequence to cDNA tab and the database
        # builder read (transcripts.default_refseq_dir finds them in the
        # bundle); the repository's pinned copies.
        ("refseq/hg19.ncbiRefSeqCurated.txt.gz", "refseq"),
        ("refseq/hg38.ncbiRefSeqCurated.txt.gz", "refseq"),
        # The WT validation sets a locally built database is verified with
        # (refdb_wizard.validation_set_path looks in the bundle).
        ("validation-sets/validation-set-hg19.h5", "validation-sets"),
        ("validation-sets/validation-set-hg38.h5", "validation-sets"),
        # The per-gene WT concordance the tab shows (wt_concordance.for_gene).
        ("spliceai_pipeline/data/wt_concordance_hg19.tsv.gz", "spliceai_pipeline/data"),
        ("spliceai_pipeline/data/wt_concordance_hg38.tsv.gz", "spliceai_pipeline/data"),
    ],
    # A user-installed spliceai imports pkg_resources at runtime; nothing in
    # this app imports it directly, so make sure it's in the bundle.
    hiddenimports=["pkg_resources"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Even if spliceai is installed in the build environment (e.g. for testing
    # from source), it must not be bundled.
    excludes=["spliceai"],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="SpliceAI-VariantScoring",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    icon="spliceai_gui/assets/logo.ico",
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="SpliceAI-VariantScoring",
)
