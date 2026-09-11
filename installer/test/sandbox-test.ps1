<#
  Runs INSIDE Windows Sandbox (started by run-sandbox.ps1). Installs Setup.exe
  from C:\USB silently, checks the install, scores test VCFs with the
  installed exe (headless --cli mode), launches the GUI twice from the
  desktop shortcut, then uninstalls -- recording every check in
  C:\Results\results.json and summary.txt, and whether the app opened any
  network connection.
#>
param([string]$Scenario = "online")

$ErrorActionPreference = "Continue"
$R = "C:\Results"
$logFile = "$R\test-log.txt"
$checks = [ordered]@{}
$AppId = "{8F3C2A51-6D4E-4B7A-9C1E-2A5B7D9E4F10}_is1"
$app = "$env:LOCALAPPDATA\Programs\SpliceAI-VariantScoring"
$exe = "$app\SpliceAI-VariantScoring.exe"
$desktopLnk = "$env:USERPROFILE\Desktop\SpliceAI Variant Scoring.lnk"
$startLnk = "$env:APPDATA\Microsoft\Windows\Start Menu\Programs\SpliceAI Variant Scoring.lnk"
$packagesDir = "$env:USERPROFILE\.spliceai_gui\python-packages"
$expectSpliceAI = $Scenario -ne "offline"

function Log($m) { Add-Content $logFile ("[{0:HH:mm:ss}] {1}" -f (Get-Date), $m) }
function Check($name, $ok, $detail = "") {
    $ok = [bool]$ok
    $checks[$name] = [ordered]@{ ok = $ok; detail = "$detail" }
    Log ("{0} {1} {2}" -f $(if ($ok) { "PASS" } else { "FAIL" }), $name, $detail)
}

Add-Type -AssemblyName UIAutomationClient, UIAutomationTypes
Add-Type @"
using System; using System.Text; using System.Collections.Generic; using System.Runtime.InteropServices;
public static class Win {
  delegate bool EnumProc(IntPtr h, IntPtr l);
  [DllImport("user32.dll")] static extern bool EnumWindows(EnumProc f, IntPtr l);
  [DllImport("user32.dll")] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] static extern bool IsWindowVisible(IntPtr h);
  [DllImport("user32.dll")] static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
  public static List<string> Titles(uint pid) {
    var r = new List<string>();
    EnumWindows((h, l) => { uint p; GetWindowThreadProcessId(h, out p);
      if (p == pid && IsWindowVisible(h)) { var sb = new StringBuilder(512); GetWindowText(h, sb, 512); if (sb.Length > 0) r.Add(sb.ToString()); }
      return true; }, IntPtr.Zero);
    return r;
  }
}
"@

# Remote TCP endpoints a process (tree) has open right now, loopback excluded.
function Remote-Connections([int[]]$pids) {
    Get-NetTCPConnection -ErrorAction SilentlyContinue |
        Where-Object { $pids -contains $_.OwningProcess -and $_.RemoteAddress -notmatch '^(127\.|::1|0\.0\.0\.0|::$)' } |
        ForEach-Object { "$($_.RemoteAddress):$($_.RemotePort)" }
}
function Tree-Pids([int]$rootPid) {
    $all = Get-CimInstance Win32_Process | Select-Object ProcessId, ParentProcessId
    $pids = @($rootPid); $added = $true
    while ($added) { $added = $false
        foreach ($p in $all) { if ($pids -contains $p.ParentProcessId -and $pids -notcontains $p.ProcessId) { $pids += $p.ProcessId; $added = $true } } }
    $pids
}

# Runs the installed exe headless; returns @{ exit; seconds; remote }.
function Run-Cli($name, [string[]]$cliArgs, [int]$timeoutSec = 3600) {
    $env:SPLICEAI_CLI_LOG = "$R\$name.log"
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $p = Start-Process $exe -ArgumentList (@("--cli") + $cliArgs) -PassThru
    $remote = @(); $nextBeat = 300
    while (-not $p.HasExited -and $sw.Elapsed.TotalSeconds -lt $timeoutSec) {
        $remote += Remote-Connections (Tree-Pids $p.Id); Start-Sleep -Milliseconds 1000
        # Heartbeat, so the host doesn't mistake a long scoring run for a stall.
        if ($sw.Elapsed.TotalSeconds -ge $nextBeat) { Log "$name still running ($([int]$sw.Elapsed.TotalMinutes) min)"; $nextBeat += 300 }
    }
    if (-not $p.HasExited) { $p | Stop-Process -Force; return @{ exit = "timeout"; seconds = [int]$sw.Elapsed.TotalSeconds; remote = $remote } }
    $p.WaitForExit()
    @{ exit = $p.ExitCode; seconds = [int]$sw.Elapsed.TotalSeconds; remote = ($remote | Sort-Object -Unique) }
}

function Click-Button($windowTitle, $buttonName) {
    $root = [Windows.Automation.AutomationElement]::RootElement
    $win = $root.FindFirst([Windows.Automation.TreeScope]::Children,
        (New-Object Windows.Automation.PropertyCondition([Windows.Automation.AutomationElement]::NameProperty, $windowTitle)))
    if (-not $win) { return $false }
    $btn = $win.FindFirst([Windows.Automation.TreeScope]::Descendants,
        (New-Object Windows.Automation.PropertyCondition([Windows.Automation.AutomationElement]::NameProperty, $buttonName)))
    if ($btn) { try { $btn.GetCurrentPattern([Windows.Automation.InvokePattern]::Pattern).Invoke(); return $true } catch {} }
    # Fallback: it's the dialog's default button.
    $shell = New-Object -ComObject WScript.Shell
    if ($shell.AppActivate($windowTitle)) { Start-Sleep 1; $shell.SendKeys("{ENTER}"); return $true }
    $false
}

# Launches the GUI from the desktop shortcut, accepts the startup notice, and
# reports the windows it showed and any network connection it opened.
function Launch-Gui($label) {
    Start-Process $desktopLnk
    $proc = $null; $titles = @(); $remote = @(); $accepted = $false; $declined = $false
    $dialogs = "Welcome to SpliceAI Variant Scoring", "SpliceAI Setup Required", "Reference data"
    for ($i = 0; $i -lt 120; $i++) {
        Start-Sleep 1
        if (-not $proc) { $proc = Get-Process SpliceAI-VariantScoring -ErrorAction SilentlyContinue | Select-Object -First 1; continue }
        $now = [Win]::Titles([uint32]$proc.Id); $titles += $now
        $remote += Remote-Connections (Tree-Pids $proc.Id)
        if (-not $accepted -and $now -contains "Welcome to SpliceAI Variant Scoring") {
            Start-Sleep 2; $accepted = Click-Button "Welcome to SpliceAI Variant Scoring" "I have read and accept the terms"
        }
        # Only when SpliceAI is missing (offline Setup): the app asks before the
        # main window opens -- record it, then continue without installing.
        if (-not $declined -and $now -contains "SpliceAI Setup Required") {
            Start-Sleep 2; $declined = Click-Button "SpliceAI Setup Required" "Continue without SpliceAI"
        }
        if ($accepted -and ($now | Where-Object { $_ -notin $dialogs })) { Start-Sleep 8
            $titles += [Win]::Titles([uint32]$proc.Id); $remote += Remote-Connections (Tree-Pids $proc.Id); break }
    }
    $titles = $titles | Sort-Object -Unique
    Log "$label windows: $($titles -join ' | ')"
    if ($proc) { Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue; Start-Sleep 3 }
    @{ started = [bool]$proc; titles = $titles; remote = ($remote | Sort-Object -Unique) }
}

# Data rows of an output TSV, skipping blank lines (Import-Csv turns a blank
# line into a row whose every field is empty).
function Read-Tsv($path) { @(Import-Csv $path -Delimiter "`t" | Where-Object { $_.chrom }) }

function Compare-Tsv($actualPath, $expectedPath) {
    $a = Read-Tsv $actualPath; $e = Read-Tsv $expectedPath
    if ($a.Count -ne $e.Count) { return "row count $($a.Count) vs expected $($e.Count)" }
    $cols = $e[0].PSObject.Properties.Name
    $diffs = 0
    for ($i = 0; $i -lt $e.Count; $i++) { foreach ($c in $cols) { if ("$($a[$i].$c)" -ne "$($e[$i].$c)") { $diffs++; if ($diffs -le 5) { Log "  diff row $i $c : '$($a[$i].$c)' vs '$($e[$i].$c)'" } } } }
    if ($diffs) { "$diffs differing cells" } else { "" }
}

Log "Scenario: $Scenario"

# --- the clean machine ------------------------------------------------------
$py = (Get-Command python, python3, py -ErrorAction SilentlyContinue | Where-Object { $_.Source -notmatch 'WindowsApps' } | ForEach-Object Source) -join ", "
Check "env.no_python" (-not $py) $py
Check "env.no_vcredist" (-not (Test-Path "HKLM:\SOFTWARE\Microsoft\VisualStudio\14.0\VC\Runtimes\x64") -and -not (Test-Path "$env:SystemRoot\System32\msvcp140_1.dll")) "msvcp140_1 in System32: $(Test-Path "$env:SystemRoot\System32\msvcp140_1.dll")"
$online = $false; try { Invoke-WebRequest -UseBasicParsing https://pypi.org/simple/ -TimeoutSec 10 | Out-Null; $online = $true } catch {}
Check "env.network_as_expected" ($online -eq ($Scenario -notlike "offline*")) "internet reachable: $online"

# --- 1. install -------------------------------------------------------------
$setupExe = (Get-ChildItem C:\USB -Filter "SpliceAI-VariantScoring-Setup-*.exe" | Select-Object -First 1).FullName
$sw = [Diagnostics.Stopwatch]::StartNew()
$p = Start-Process $setupExe -ArgumentList "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART", "/LOG=`"$R\setup.log`"" -Wait -PassThru
Check "install.exit_code" ($p.ExitCode -eq 0) "exit $($p.ExitCode) after $([int]$sw.Elapsed.TotalSeconds) s"
Check "install.exe" (Test-Path $exe)
Check "install.vc_runtime_app_local" (Test-Path "$app\msvcp140_1.dll")
Check "install.snpeff_java_mane" ((Test-Path "$app\snpeff\snpEff\data\GRCh37.p13\snpEffectPredictor.bin") -and (Test-Path "$app\snpeff\snpEff\data\GRCh38.p14\snpEffectPredictor.bin") -and (Get-ChildItem "$app\java" -Recurse -Filter java.exe) -and (Get-ChildItem "$app\mane" -Filter *.gz))
$hasSpliceAI = (Test-Path "$app\_internal\spliceai\models\spliceai5.h5") -and (Test-Path "$app\_internal\spliceai-1.3.1.dist-info\METADATA")
Check "install.spliceai_as_expected" ($hasSpliceAI -eq $expectSpliceAI) "installed: $hasSpliceAI, expected: $expectSpliceAI"
$setupLog = Get-Content "$R\setup.log" -Raw
switch ($Scenario) {
    "online" { Check "install.spliceai_source" ($setupLog -match "Downloading|Download") "downloaded from PyPI" }
    "offline" { Check "install.offline_handled" ($setupLog -match "SpliceAI download failed") "download failure logged, Setup finished" }
    "offline-usb-wheel" { Check "install.spliceai_source" ($setupLog -match "Using local SpliceAI wheel") "used wheel from C:\USB" }
}
$wsh = New-Object -ComObject WScript.Shell
Check "install.desktop_shortcut" ((Test-Path $desktopLnk) -and $wsh.CreateShortcut($desktopLnk).TargetPath -eq $exe)
Check "install.startmenu_shortcut" ((Test-Path $startLnk) -and $wsh.CreateShortcut($startLnk).TargetPath -eq $exe)
Check "install.uninstall_entry" (Test-Path "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\$AppId") (Get-ItemProperty "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\$AppId" -ErrorAction SilentlyContinue).DisplayName
$spliceaiStamp = if ($hasSpliceAI) { (Get-Item "$app\_internal\spliceai\utils.py").LastWriteTimeUtc } else { $null }

# --- 2. score -----------------------------------------------------------------
$cli = Run-Cli "ankrd26" @("C:\TestData\ankrd26.vcf", "--build", "hg19", "--mode", "masked", "--fasta", "C:\Ref\hg19.fa", "-o", "C:\Results\ankrd26.tsv") 900
if ($expectSpliceAI) {
    $row = if (Test-Path "$R\ankrd26.tsv") { Read-Tsv "$R\ankrd26.tsv" | Select-Object -First 1 }
    Check "score.ankrd26_DS_AG_0.29" ($cli.exit -eq 0 -and $row.gene -eq "ANKRD26" -and $row.DS_AG -eq "0.29") "exit $($cli.exit), $($cli.seconds) s, gene=$($row.gene) DS_AG=$($row.DS_AG)"
} else {
    Check "score.without_spliceai_fails_cleanly" ($cli.exit -eq 1) "exit $($cli.exit) (expected 1: SpliceAI not installed)"
}
Check "score.ankrd26_no_network" (-not $cli.remote) ($cli.remote -join ", ")

if ($expectSpliceAI -and (Test-Path C:\TestData\extra.vcf)) {
    $cli = Run-Cli "extra" @("C:\TestData\extra.vcf", "--build", "hg38", "--mode", "masked", "--fasta", "C:\Ref\hg38.fa", "--use-snpeff", "-o", "C:\Results\extra.tsv")
    $rows = if (Test-Path "$R\extra.tsv") { Read-Tsv "$R\extra.tsv" } else { @() }
    $withSnpEff = @($rows | Where-Object snpeff_transcript -ne "").Count
    Check "score.extra_vcf" ($cli.exit -eq 0 -and $rows.Count -gt 0) "exit $($cli.exit), $($cli.seconds) s, $($rows.Count) rows, $withSnpEff with SnpEff transcript"
    Check "score.extra_snpeff_offline" ($withSnpEff -gt 0)
    if (Test-Path C:\TestData\extra-expected.tsv) {
        $diff = Compare-Tsv "$R\extra.tsv" "C:\TestData\extra-expected.tsv"
        Check "score.extra_matches_known_good" (-not $diff) $diff
    }
    Check "score.extra_no_network" (-not $cli.remote) ($cli.remote -join ", ")
}

# --- 3. GUI, twice ------------------------------------------------------------
foreach ($n in 1, 2) {
    $g = Launch-Gui "launch$n"
    Check "gui.launch${n}_started" ($g.started -and ($g.titles | Where-Object { $_ -notin "Welcome to SpliceAI Variant Scoring", "SpliceAI Setup Required", "Reference data" })) ($g.titles -join " | ")
    $setupShown = $g.titles -contains "SpliceAI Setup Required"
    Check "gui.launch${n}_setup_dialog_as_expected" ($setupShown -ne $expectSpliceAI) "SpliceAI Setup dialog shown: $setupShown"
    Check "gui.launch${n}_no_network" (-not $g.remote) ($g.remote -join ", ")
    if ($n -eq 1) { Check "gui.reference_note_shown" ($g.titles -contains "Reference data") }
}
if ($hasSpliceAI) {
    Check "relaunch.nothing_reinstalled" (((Get-Item "$app\_internal\spliceai\utils.py").LastWriteTimeUtc -eq $spliceaiStamp) -and -not (Get-ChildItem $packagesDir -ErrorAction SilentlyContinue)) "python-packages empty, SpliceAI files untouched"
}

# --- 4. uninstall ---------------------------------------------------------------
Get-Process SpliceAI-VariantScoring -ErrorAction SilentlyContinue | Stop-Process -Force
Start-Process "$app\unins000.exe" -ArgumentList "/VERYSILENT", "/SUPPRESSMSGBOXES", "/NORESTART" -Wait
for ($i = 0; $i -lt 90 -and (Test-Path $app); $i++) { Start-Sleep 2 }
Check "uninstall.folder_removed" (-not (Test-Path $app))
Check "uninstall.shortcuts_removed" (-not (Test-Path $desktopLnk) -and -not (Test-Path $startLnk))
Check "uninstall.entry_removed" (-not (Test-Path "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\$AppId"))

# --- report -----------------------------------------------------------------------
# Summary first (plain text), JSON second; DONE is written whatever happens so
# the host never waits out its whole timeout.
Log "All checks done, writing report"
try {
    $names = @($checks.Keys)
    $failed = @($names | Where-Object { -not $checks[$_].ok })
    $lines = @("Scenario: $Scenario -- $($names.Count - $failed.Count)/$($names.Count) checks passed")
    foreach ($n in $names) { $lines += "{0}  {1,-42} {2}" -f $(if ($checks[$n].ok) { "PASS" } else { "FAIL" }), $n, $checks[$n].detail }
    [IO.File]::WriteAllLines("$R\summary.txt", [string[]]$lines)
    Log "summary.txt written"
    $plain = [ordered]@{}; foreach ($n in $names) { $plain[$n] = New-Object PSObject -Property $checks[$n] }
    [IO.File]::WriteAllText("$R\results.json", (New-Object PSObject -Property $plain | ConvertTo-Json -Depth 4))
    Log "results.json written"
} catch {
    Log "Report failed: $($_.Exception.Message)"
} finally {
    [IO.File]::WriteAllText("$R\DONE", "done")
}
