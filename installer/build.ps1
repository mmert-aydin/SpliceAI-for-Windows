<#
.SYNOPSIS
  Builds the SpliceAI Variant Scoring Setup.exe from source.

.DESCRIPTION
  1. Makes sure the build venv (.venv, Python 3.13) exists with
     requirements.txt + requirements-build.txt, and that spliceai is NOT in it.
  2. Builds the app with PyInstaller into installer\build\dist (the repo's own
     dist\ folder is left alone).
  3. Compiles SpliceAI-VariantScoring.iss with Inno Setup 6, passing the
     SpliceAI wheel URL/hash from spliceai_gui/spliceai_setup.py so Setup and
     the app always agree on the exact file.
  Output: installer\output\SpliceAI-VariantScoring-Setup-<Version>.exe

  On a fresh checkout, run installer\fetch-thirdparty.ps1 first (SnpEff,
  Java, MANE, VC++ runtime DLLs into installer\thirdparty).

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File installer\build.ps1 -Version 1.2.0
#>
param(
    [string]$Version = "1.2.0",
    # Reuse the last PyInstaller output (only re-compile the installer).
    [switch]$SkipAppBuild
)

$ErrorActionPreference = "Stop"
$here = $PSScriptRoot
$repo = Split-Path $here -Parent
$venvPy = Join-Path $repo ".venv\Scripts\python.exe"
$third = Join-Path $here "thirdparty"
$buildDir = Join-Path $here "build"
$appDir = Join-Path $buildDir "dist\SpliceAI-VariantScoring"

function Step($message) { Write-Host "==> $message" -ForegroundColor Cyan }
function Assert-Exit($what) { if ($LASTEXITCODE -ne 0) { throw "$what failed (exit code $LASTEXITCODE)" } }

# --- prerequisites ---------------------------------------------------------
$iscc = @(
    "$env:LOCALAPPDATA\Programs\Inno Setup 6\ISCC.exe",
    "${env:ProgramFiles(x86)}\Inno Setup 6\ISCC.exe",
    "$env:ProgramFiles\Inno Setup 6\ISCC.exe"
) | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $iscc) { throw "Inno Setup 6 not found. Install it with: winget install --id JRSoftware.InnoSetup -e --scope user" }

$required = "vcruntime\msvcp140_1.dll", "java", "snpeff\snpEff\snpEff.jar",
            "snpeff\snpEff\data\GRCh37.p13\snpEffectPredictor.bin",
            "snpeff\snpEff\data\GRCh38.p14\snpEffectPredictor.bin", "mane"
foreach ($p in $required) {
    if (-not (Test-Path (Join-Path $third $p))) { throw "Missing installer\thirdparty\$p -- run installer\fetch-thirdparty.ps1 first." }
}
# Bundled from the repository itself by the .spec.
$tracked = "refseq\hg19.ncbiRefSeqCurated.txt.gz", "refseq\hg38.ncbiRefSeqCurated.txt.gz",
           "spliceai_pipeline\data\wt_concordance_hg19.tsv.gz", "spliceai_pipeline\data\wt_concordance_hg38.tsv.gz",
           "validation-sets\validation-set-hg19.h5", "validation-sets\validation-set-hg38.h5"
foreach ($p in $tracked) {
    if (-not (Test-Path (Join-Path $repo $p))) { throw "Missing $p -- it is part of the repository; check out the full tree." }
}

# --- 1. build venv ---------------------------------------------------------
# parasail (spliceai_pipeline/align.py) goes in from its pinned wheel, checked
# by SHA-256, before requirements.txt (which then finds it already there).
$ParasailWheel = "parasail-1.3.4-py2.py3-none-win_amd64.whl"
$ParasailUrl = "https://files.pythonhosted.org/packages/7d/3d/3fb505271ca157798326c3f4c357ccde2c1c6f05f05400aa19e3e55bbd01/$ParasailWheel"
$ParasailSha256 = "25B8260B922933C8E7E8CE008DDCBBFF4EF998B7D077169ED441D70AB7A78B5A"
function Install-Parasail {
    & $venvPy -c "import importlib.metadata as m, sys; sys.exit(0 if m.version('parasail') == '1.3.4' else 1)" 2>$null
    if ($LASTEXITCODE -eq 0) { return }
    $file = Join-Path ([IO.Path]::GetTempPath()) $ParasailWheel
    Invoke-WebRequest -UseBasicParsing $ParasailUrl -OutFile $file
    $actual = (Get-FileHash $file -Algorithm SHA256).Hash
    if ($actual -ne $ParasailSha256) { throw "$ParasailWheel has SHA-256 $actual, expected $ParasailSha256." }
    & $venvPy -m pip install --no-deps --no-index $file; Assert-Exit "pip install parasail"
}
if (-not (Test-Path $venvPy)) {
    Step "Creating build venv (.venv)"
    py -3.13 -m venv (Join-Path $repo ".venv"); Assert-Exit "venv"
    Install-Parasail
    & $venvPy -m pip install -r (Join-Path $repo "requirements.txt") -r (Join-Path $repo "requirements-build.txt"); Assert-Exit "pip install"
}
Install-Parasail
& $venvPy -c "import importlib.util, sys; sys.exit(1 if importlib.util.find_spec('spliceai') else 0)"
if ($LASTEXITCODE -ne 0) {
    throw ("SpliceAI is installed in .venv. Remove it with: .venv\Scripts\python -m pip uninstall -y spliceai`n" +
           "Why: SpliceAI's license forbids redistributing it, and PyInstaller would pack any copy it finds in .venv into the build. " +
           "The program doesn't need that copy, from source either: it finds SpliceAI in %USERPROFILE%\.spliceai_gui\python-packages.")
}

# --- 2. app ----------------------------------------------------------------
Push-Location $repo
try {
    if (-not $SkipAppBuild) {
        Step "Building the app with PyInstaller"
        & $venvPy -m PyInstaller SpliceAI-VariantScoring.spec --noconfirm --log-level WARN `
            --distpath (Join-Path $buildDir "dist") --workpath (Join-Path $buildDir "work")
        Assert-Exit "PyInstaller"
    }
    if (Test-Path (Join-Path $appDir "_internal\spliceai")) { throw "spliceai ended up inside the app build -- it must never be shipped." }

    $wheel = & $venvPy -c "from spliceai_gui import spliceai_setup as s; print(s.SPLICEAI_VERSION); print(s.SPLICEAI_WHEEL_NAME); print(s.SPLICEAI_WHEEL_URL); print(s.SPLICEAI_WHEEL_SHA256)"
    Assert-Exit "reading SpliceAI wheel info"
} finally { Pop-Location }

# --- 3. Setup.exe ----------------------------------------------------------
Step "Compiling Setup.exe with Inno Setup (this takes a while)"
& $iscc /Qp "/DAppVersion=$Version" "/DAppDir=$appDir" "/DThirdPartyDir=$third" `
    "/DSpliceAIVersion=$($wheel[0])" "/DSpliceAIWheelName=$($wheel[1])" `
    "/DSpliceAIWheelUrl=$($wheel[2])" "/DSpliceAIWheelSha256=$($wheel[3])" `
    (Join-Path $here "SpliceAI-VariantScoring.iss")
Assert-Exit "ISCC"

$setup = Get-Item (Join-Path $here "output\SpliceAI-VariantScoring-Setup-$Version.exe")
Step ("Done: {0}  ({1:N0} bytes, {2:N0} MB)" -f $setup.FullName, $setup.Length, ($setup.Length / 1MB))
Write-Host ("SHA-256: " + (Get-FileHash $setup -Algorithm SHA256).Hash)
