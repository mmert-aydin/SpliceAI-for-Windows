# -*- mode: python ; coding: utf-8 -*-
# Build the Windows app (onedir):
#   python -m PyInstaller SpliceAI-VariantScoring.spec --noconfirm
# Output: dist/SpliceAI-VariantScoring/SpliceAI-VariantScoring.exe
#
# Reconstructed to match the shipped build: windowed exe, no version resource,
# logo.ico icon, assets + LICENSE-THIRD-PARTY.md as data, and SpliceAI itself
# NOT bundled (license) -- the app offers to install it on first run.

a = Analysis(
    ["run_app.py"],
    pathex=[],
    binaries=[],
    datas=[
        # assets_paths.py finds these relative to its own file, so keep the
        # source-tree layout inside the bundle.
        ("spliceai_gui/assets/logo.png", "spliceai_gui/assets"),
        ("spliceai_gui/assets/logo.ico", "spliceai_gui/assets"),
        ("LICENSE-THIRD-PARTY.md", "."),
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
