<#
.SYNOPSIS
  Tests the built Setup.exe in Windows Sandbox: a clean, throw-away Windows
  with no Python, no VC++ runtime and nothing from this PC.

.DESCRIPTION
  Scenarios:
    online             network on: install, run, relaunch, uninstall (the normal case)
    offline            network off: Setup can't download SpliceAI and must say so
                       and finish anyway; everything else must work offline
    offline-usb-wheel  network off, the SpliceAI wheel next to Setup.exe on the
                       "USB drive": Setup must use it
    usb-reference      network on, a reference-data folder next to Setup.exe on the
                       "USB drive": Setup must copy the genomes onto the PC, the
                       program must find them by itself and never show the
                       "Reference data" note, and an uninstall must keep them
    clickthrough       network on: clicks through the real Setup wizard and the
                       program like a person (SpliceAI check, a real Run, normal
                       uninstall), with a screenshot of every step (clickthrough.ps1;
                       -MemoryMB 6144 is enough)
    manual             network on, nothing automated -- click through it yourself

  Mapped into the sandbox (read-only unless noted):
    C:\USB       Setup.exe (+ the wheel for offline-usb-wheel)   -- the "flash drive"
    C:\Ref       reference FASTAs (hg19.fa/hg38.fa + .fai)
    C:\TestData  test VCFs + known-good outputs
    C:\Harness   this folder
    C:\Results   writable; installer\test\results\<scenario>

  Needs the Windows Sandbox feature (Containers-DisposableClientVM).

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File installer\test\run-sandbox.ps1 -Scenario online `
    -ExtraVcf $HOME\Desktop\my-variants.vcf -ExtraExpected <known-good my-variants.tsv>
#>
param(
    [ValidateSet("online", "offline", "offline-usb-wheel", "usb-reference", "clickthrough", "manual")]
    [string]$Scenario = "online",
    [string]$Setup,
    [string]$ReferenceDir = "$env:USERPROFILE\SpliceAI_reference_data",
    # Optional larger hg38 VCF (not committed: may be patient data) and its
    # known-good output from a source run with the same options.
    [string]$ExtraVcf,
    [string]$ExtraExpected,
    [string]$Wheel,
    [int]$TimeoutMinutes = 60,
    # Memory for the sandbox. Scoring a large VCF needs the default; the
    # click-through (one variant) is fine with 6144, which leaves more for
    # the PC itself.
    [int]$MemoryMB = 8192
)

$ErrorActionPreference = "Stop"
$installerDir = Split-Path $PSScriptRoot -Parent
if (-not $Setup) {
    $Setup = (Get-ChildItem (Join-Path $installerDir "output") -Filter "SpliceAI-VariantScoring-Setup-*.exe" |
              Sort-Object LastWriteTime -Descending | Select-Object -First 1).FullName
}
if (-not (Test-Path "$env:SystemRoot\System32\WindowsSandbox.exe")) {
    throw "Windows Sandbox isn't enabled. In an administrator PowerShell run: Enable-WindowsOptionalFeature -Online -FeatureName Containers-DisposableClientVM -All  -- then restart."
}

# Fresh folders for every run: nothing is reused, so a previous run's logs
# can never be mixed in, and nothing has to be deleted first (folders a
# sandbox has mapped can stay locked on the host for a while after it closes).
$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$work = Join-Path $env:TEMP "spliceai-sandbox\$Scenario-$stamp"
$results = Join-Path $PSScriptRoot "results\$Scenario-$stamp"
Write-Host "Results: $results"
$usb = New-Item -ItemType Directory -Force (Join-Path $work "usb")
$data = New-Item -ItemType Directory -Force (Join-Path $work "data")
$harness = New-Item -ItemType Directory -Force (Join-Path $work "harness")
New-Item -ItemType Directory -Force $results | Out-Null
# A copy, not this folder itself: results\ lives inside it, and mapping a
# folder and one of its subfolders separately can make the sandbox refuse to start.
Copy-Item (Join-Path $PSScriptRoot "sandbox-test.ps1"), (Join-Path $PSScriptRoot "clickthrough.ps1") $harness

# Hard links: no 600 MB copy (same volume).
New-Item -ItemType HardLink -Path (Join-Path $usb (Split-Path $Setup -Leaf)) -Target $Setup | Out-Null
if ($Scenario -eq "offline-usb-wheel") {
    if (-not $Wheel) { throw "-Wheel <spliceai-1.3.1-py2.py3-none-any.whl> is required for offline-usb-wheel" }
    Copy-Item $Wheel $usb
}
if ($Scenario -eq "usb-reference") {
    # The drive as a colleague gets it: reference-data beside Setup.exe. Hard
    # links where possible, so putting 6 GB "on the drive" costs nothing and
    # the files are byte-identical to the real ones.
    $refOnUsb = New-Item -ItemType Directory -Force (Join-Path $usb "reference-data")
    foreach ($name in "hg19.fa", "hg19.fa.fai", "hg38.fa", "hg38.fa.fai") {
        $src = Join-Path $ReferenceDir $name
        if (-not (Test-Path $src)) { throw "Missing $src -- the usb-reference scenario needs both genomes." }
        $dest = Join-Path $refOnUsb $name
        try { New-Item -ItemType HardLink -Path $dest -Target $src -ErrorAction Stop | Out-Null }
        catch { Copy-Item $src $dest }   # a different volume: no hard links
    }
}
Copy-Item (Join-Path $PSScriptRoot "data\*") $data
if ($ExtraVcf) { Copy-Item $ExtraVcf (Join-Path $data "extra.vcf") }
if ($ExtraExpected) { Copy-Item $ExtraExpected (Join-Path $data "extra-expected.tsv") }

function Map($hostDir, $sandboxDir, $readOnly) {
    "    <MappedFolder><HostFolder>$hostDir</HostFolder><SandboxFolder>$sandboxDir</SandboxFolder><ReadOnly>$($readOnly.ToString().ToLower())</ReadOnly></MappedFolder>"
}
$networking = if ($Scenario -like "offline*") { "Disable" } else { "Default" }
$logon = switch ($Scenario) {
    "manual" { "" }
    # Clicks through the real wizard and window like a person, with screenshots.
    "clickthrough" { "  <LogonCommand><Command>powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Minimized -File C:\Harness\clickthrough.ps1</Command></LogonCommand>" }
    default { "  <LogonCommand><Command>powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Minimized -File C:\Harness\sandbox-test.ps1 -Scenario $Scenario</Command></LogonCommand>" }
}
$wsb = @"
<Configuration>
  <Networking>$networking</Networking>
  <MemoryInMB>$MemoryMB</MemoryInMB>
  <MappedFolders>
$(Map $usb "C:\USB" $true)
$(Map $ReferenceDir "C:\Ref" $true)
$(Map $data "C:\TestData" $true)
$(Map $harness "C:\Harness" $true)
$(Map $results "C:\Results" $false)
  </MappedFolders>
$logon
</Configuration>
"@
$wsbPath = Join-Path $work "$Scenario.wsb"
[IO.File]::WriteAllText($wsbPath, $wsb, (New-Object Text.UTF8Encoding $false))

Write-Host "Starting Windows Sandbox ($Scenario, networking: $networking) ..."
# Call WindowsSandbox.exe directly: the .wsb file type isn't always associated
# with it (it wasn't on the build PC), and then opening the file does nothing.
Start-Process "$env:SystemRoot\System32\WindowsSandbox.exe" -ArgumentList "`"$wsbPath`""
if ($Scenario -eq "manual") { Write-Host "Manual session: Setup is in C:\USB inside the sandbox."; return }

$deadline = (Get-Date).AddMinutes($TimeoutMinutes)
$testLog = Join-Path $results "test-log.txt"
$stalled = $false
while (-not (Test-Path (Join-Path $results "DONE")) -and (Get-Date) -lt $deadline) {
    Start-Sleep 10
    # The longest silent step (scoring a large VCF) logs nothing for a while,
    # so only give up after a long gap.
    if ((Test-Path $testLog) -and ((Get-Date) - (Get-Item $testLog).LastWriteTime).TotalMinutes -gt 25) { $stalled = $true; break }
}
Get-Process WindowsSandbox, WindowsSandboxClient -ErrorAction SilentlyContinue | Stop-Process -Force
$summary = Join-Path $results "summary.txt"
if (Test-Path $summary) { Get-Content $summary; return }
$why = if ($stalled) { "no progress for 25 min" } elseif (Test-Path (Join-Path $results "DONE")) { "report step failed" } else { "timed out after $TimeoutMinutes min" }
Write-Warning "No summary.txt ($why); checks logged so far:"
if (Test-Path $testLog) { Select-String -Path $testLog -Pattern "PASS|FAIL|Report failed" | ForEach-Object Line }
