<#
.SYNOPSIS
  Runs every check in tests\ and reports which scripts failed. See tests\README.md.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File tests\run-tests.ps1
#>
$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$py = Join-Path $repo ".venv\Scripts\python.exe"
if (-not (Test-Path $py)) { throw "No build venv (.venv) -- run installer\build.ps1 once; it creates it." }

# Exit codes decide pass/fail. (With "Stop", Windows PowerShell 5.1 can abort
# on TensorFlow's harmless stderr messages when the output is redirected.)
$ErrorActionPreference = "Continue"
$failed = @()
foreach ($test in "test_pipeline.py", "test_downloads.py", "test_reference_locator.py", "test_build_selection.py", "test_drag_drop.py", "test_gui_features.py") {
    Write-Host "==> $test" -ForegroundColor Cyan
    & $py -W ignore (Join-Path $PSScriptRoot $test)
    if ($LASTEXITCODE -ne 0) { $failed += $test }
}
if ($failed) {
    Write-Host "FAILED: $($failed -join ', ')" -ForegroundColor Red
    exit 1
}
Write-Host "All checks passed." -ForegroundColor Green
