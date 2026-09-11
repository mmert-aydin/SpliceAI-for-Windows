<#
  Runs INSIDE Windows Sandbox (run-sandbox.ps1 -Scenario clickthrough), network
  on. Does what a person would: clicks through the real Setup wizard, starts
  the program, uses it (SpliceAI check, a real run with the Run button), and
  uninstalls it normally -- with a screenshot of every step in
  C:\Results\shots and PASS/FAIL lines in C:\Results\test-log.txt.

  Every click is tried three ways -- UI Automation, a real mouse click on the
  control, then the keyboard -- and each is checked by waiting for what should
  come next, since Inno Setup's controls don't all answer UI Automation.
#>
$ErrorActionPreference = "Continue"
$R = "C:\Results"
$shots = New-Item -ItemType Directory -Force "$R\shots"
$logFile = "$R\test-log.txt"
$checks = [ordered]@{}
$app = "$env:LOCALAPPDATA\Programs\SpliceAI-VariantScoring"

function Log($m) { Add-Content $logFile ("[{0:HH:mm:ss}] {1}" -f (Get-Date), $m) }
function Check($name, $ok, $detail = "") {
    $ok = [bool]$ok
    $checks[$name] = [ordered]@{ ok = $ok; detail = "$detail" }
    Log ("{0} {1} {2}" -f $(if ($ok) { "PASS" } else { "FAIL" }), $name, $detail)
}

Add-Type -AssemblyName UIAutomationClient, UIAutomationTypes, System.Windows.Forms, System.Drawing
Add-Type @"
using System; using System.Runtime.InteropServices;
public static class Mouse {
  [DllImport("user32.dll")] static extern bool SetCursorPos(int x, int y);
  [DllImport("user32.dll")] static extern void mouse_event(uint f, uint x, uint y, uint d, UIntPtr e);
  public static void Click(int x, int y) { SetCursorPos(x, y); mouse_event(2, 0, 0, 0, UIntPtr.Zero); mouse_event(4, 0, 0, 0, UIntPtr.Zero); }
}
"@
$AE = [Windows.Automation.AutomationElement]
$TS = [Windows.Automation.TreeScope]
$AnyCond = [Windows.Automation.Condition]::TrueCondition
$shell = New-Object -ComObject WScript.Shell

$shotNo = 0
function Shot($name) {
    $script:shotNo++
    $b = [System.Windows.Forms.Screen]::PrimaryScreen.Bounds
    $bmp = New-Object System.Drawing.Bitmap $b.Width, $b.Height
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.CopyFromScreen($b.Location, [System.Drawing.Point]::Empty, $b.Size)
    $bmp.Save(("{0}\{1:D2}-{2}.png" -f $shots, $script:shotNo, $name))
    $g.Dispose(); $bmp.Dispose()
}
function Windows { $AE::RootElement.FindAll($TS::Children, $AnyCond) }
# Top-level windows, and dialogs owned by one (UI Automation lists the
# program's message boxes as children of its main window, not of the desktop).
function Find-Window($titleRegex) {
    foreach ($w in Windows) {
        if ($w.Current.Name -match $titleRegex) { return $w }
        try {
            $sub = $w.FindAll($TS::Children, $AnyCond) | Where-Object {
                $_.Current.ControlType.ProgrammaticName -eq "ControlType.Window" -and $_.Current.Name -match $titleRegex
            } | Select-Object -First 1
            if ($sub) { return $sub }
        } catch {}
    }
    $null
}
function Find-In($root, $nameRegex, $type) {
    if (-not $root) { return $null }
    try {
        $root.FindAll($TS::Descendants, $AnyCond) | Where-Object {
            $_.Current.Name -match $nameRegex -and (-not $type -or $_.Current.ControlType.ProgrammaticName -eq "ControlType.$type")
        } | Select-Object -First 1
    } catch { $null }
}
function Dump($root, $label) {
    if (-not $root) { return }
    $items = $root.FindAll($TS::Descendants, $AnyCond) | ForEach-Object { "{0}:'{1}'" -f $_.Current.ControlType.ProgrammaticName.Replace("ControlType.", ""), $_.Current.Name }
    Log "  controls on $label -- $($items -join ', ')"
}
function Wait-For([scriptblock]$probe, $seconds) {
    $end = (Get-Date).AddSeconds($seconds)
    while ((Get-Date) -lt $end) { $v = & $probe; if ($v) { return $v }; Start-Sleep -Milliseconds 700 }
    $null
}
function Click-El($el) {
    $r = $el.Current.BoundingRectangle
    [Mouse]::Click([int]($r.X + $r.Width / 2), [int]($r.Y + $r.Height / 2))
}
function Try-Invoke($el) {
    foreach ($p in [Windows.Automation.InvokePattern]::Pattern, [Windows.Automation.SelectionItemPattern]::Pattern,
                   [Windows.Automation.TogglePattern]::Pattern) {
        try {
            $pat = $el.GetCurrentPattern($p)
            if ($p -eq [Windows.Automation.InvokePattern]::Pattern) { $pat.Invoke() }
            elseif ($p -eq [Windows.Automation.SelectionItemPattern]::Pattern) { $pat.Select() }
            else { $pat.Toggle() }
            return $true
        } catch {}
    }
    $false
}
function Send-Keys($windowTitle, $keys) {
    if ($shell.AppActivate($windowTitle)) { Start-Sleep -Milliseconds 500; [System.Windows.Forms.SendKeys]::SendWait($keys); return $true }
    $false
}
# Presses a control and waits until $done is true: UI Automation, then a real
# mouse click, then the keyboard. Returns how it worked, or $null.
function Act($el, $windowTitle, $keys, [scriptblock]$done, $seconds = 10) {
    if ($el -and (Try-Invoke $el) -and (Wait-For $done $seconds)) { return "automation" }
    if ($el) { Click-El $el; if (Wait-For $done $seconds) { return "mouse" } }
    if ($keys -and (Send-Keys $windowTitle $keys) -and (Wait-For $done $seconds)) { return "keyboard $keys" }
    $null
}
$WizardTitle = 'Setup - SpliceAI Variant Scoring'
function Wizard { Find-Window '^Setup - SpliceAI' }
function Wizard-Has($textRegex) { $w = Wizard; if ($w -and (Find-In $w $textRegex)) { $w } }
# Closes the program like a person would (its dialogs first, then the window);
# stops it only if it doesn't go.
function Close-App {
    foreach ($t in '^SpliceAI check$', '^Reference data$', '^SpliceAI Setup Required$') {
        $d = Find-Window $t
        if ($d) { try { $d.GetCurrentPattern([Windows.Automation.WindowPattern]::Pattern).Close() } catch {} }
    }
    $m = Find-Window '^SpliceAI Variant Scoring$'
    if ($m) { try { $m.GetCurrentPattern([Windows.Automation.WindowPattern]::Pattern).Close() } catch {} }
    if (-not (Wait-For { -not (Get-Process SpliceAI-VariantScoring -ErrorAction SilentlyContinue) } 30)) {
        Log "  the program didn't close by itself -- stopping it"
        Get-Process SpliceAI-VariantScoring -ErrorAction SilentlyContinue | Stop-Process -Force
        Start-Sleep 2
    }
}

Log "Scenario: clickthrough"

# --- 1. the real Setup wizard --------------------------------------------------------
$setupExe = (Get-ChildItem C:\USB -Filter "SpliceAI-VariantScoring-Setup-*.exe" | Select-Object -First 1).FullName
# /LOG only writes Setup's own log (checked below for the SpliceAI download); the wizard is unchanged.
Start-Process $setupExe -ArgumentList "/LOG=`"$R\setup.log`""
$w = Wait-For { Wizard-Has 'Welcome to' } 120
Shot "setup-welcome"; Check "setup.welcome_page" $w
Dump $w "the Welcome page"

$how = Act (Find-In $w '^&?Next$') $WizardTitle "{ENTER}" { Wizard-Has 'License Agreement' }
Shot "setup-terms"; Check "setup.terms_page" $how "reached by $how"
$w = Wizard
$accept = Find-In $w 'accept the agreement'
$how = Act $accept $WizardTitle "%a" {
    $a = Find-In (Wizard) 'accept the agreement'
    try { $a.GetCurrentPattern([Windows.Automation.SelectionItemPattern]::Pattern).Current.IsSelected } catch { $false }
} 5
if (-not $how) { Send-Keys $WizardTitle "%a" | Out-Null; $how = "keyboard %a (unverified)" }
Log "  accepted the terms by $how"
Shot "setup-terms-accepted"

$how = Act (Find-In $w '^&?Next$') $WizardTitle "{ENTER}" { Wizard-Has 'Additional Tasks' }
# Inno draws its task list itself, so the "Create a desktop shortcut" box
# isn't visible to UI Automation -- it's in the screenshot.
Shot "setup-shortcuts"
Check "setup.shortcut_page" $how "reached by $how"

$how = Act (Find-In (Wizard) '^&?Next$') $WizardTitle "{ENTER}" { Wizard-Has 'Ready to Install' }
Shot "setup-ready"; Check "setup.ready_page" $how "reached by $how"

$how = Act (Find-In (Wizard) '^&?Install$') $WizardTitle "{ENTER}" {
    Wizard-Has 'Downloading SpliceAI|Installing|Completing'
} 15
Log "  Install pressed by $how"
$sawDownload = $false; $polls = 0
$done = Wait-For {
    $w = Wizard
    $script:polls += 1
    if ($w -and -not $script:sawDownload -and (Find-In $w 'Downloading SpliceAI')) {
        $script:sawDownload = $true; Shot "setup-downloading-spliceai"
    }
    if ($w -and (Find-In $w 'Completing')) { $w }
    elseif ($w -and ($script:polls % 12) -eq 1) { Shot "setup-progress"; $null }
} 900
Log "  'Downloading SpliceAI' page caught on screen: $sawDownload (it can pass within a second)"
Shot "setup-finished"; Check "setup.finish_page" $done
$how = Act (Find-In $done '^&?Finish$') $WizardTitle "{ENTER}" { -not (Wizard) }
Log "  Finish pressed by $how"
Check "install.spliceai_present" (Wait-For { Test-Path "$app\_internal\spliceai\models\spliceai5.h5" } 30)
Check "setup.downloaded_spliceai_from_pypi" (Select-String -Path "$R\setup.log" -Quiet `
    -Pattern 'Downloading temporary file from https://files\.pythonhosted\.org')

# --- 2. first start ----------------------------------------------------------------------
$termsTitle = '^Welcome to SpliceAI Variant Scoring$'
$terms = Wait-For { Find-Window $termsTitle } 180
Shot "app-terms-notice"; Check "app.starts_with_terms_notice" $terms
$how = Act (Find-In $terms 'I have read and accept') "Welcome to SpliceAI Variant Scoring" "{ENTER}" { -not (Find-Window $termsTitle) }
Log "  terms accepted by $how"
$main = Wait-For { Find-Window '^SpliceAI Variant Scoring$' } 120
$note = Wait-For { Find-Window '^Reference data$' } 30
Shot "app-first-start"; Check "app.main_window_opens" $main
Check "app.reference_note_on_first_start" $note
Check "app.no_spliceai_setup_dialog" (-not (Find-Window '^SpliceAI Setup Required$'))
if ($note) { Act (Find-In $note '^OK$' "Button") "Reference data" "{ENTER}" { -not (Find-Window '^Reference data$') } | Out-Null }

# --- 3. the SpliceAI check button ------------------------------------------------------
$btn = Find-In $main 'SpliceAI found'
Check "app.spliceai_button_green" $btn ($btn.Current.Name)
$how = Act $btn "SpliceAI Variant Scoring" $null { Find-Window '^SpliceAI check$' }
$chk = Find-Window '^SpliceAI check$'
Shot "app-spliceai-check"
Check "app.spliceai_check_all_files" ($chk -and (Find-In $chk 'Everything live scoring needs is in place')) "opened by $how"
if ($chk) { Act (Find-In $chk '^OK$' "Button") "SpliceAI check" "{ENTER}" { -not (Find-Window '^SpliceAI check$') } | Out-Null }

# --- 4. a real run with the Run button -------------------------------------------------
# Close, point the settings at the reference files (as Browse... would), start again.
Close-App
$cfgPath = "$env:USERPROFILE\.spliceai_gui\config.json"
$cfg = Get-Content $cfgPath -Raw | ConvertFrom-Json
$cfg.fasta_path_hg19 = "C:\Ref\hg19.fa"; $cfg.fasta_path_hg38 = "C:\Ref\hg38.fa"; $cfg.precomputed_mode = "none"
# Without a byte-order mark: Windows PowerShell's "UTF8" writes one, and the
# program's JSON reader then ignores the whole file.
[IO.File]::WriteAllText($cfgPath, ($cfg | ConvertTo-Json -Depth 5), (New-Object Text.UTF8Encoding $false))
Start-Process "$env:USERPROFILE\Desktop\SpliceAI Variant Scoring.lnk"
$terms = Wait-For { Find-Window $termsTitle } 120
Act (Find-In $terms 'I have read and accept') "Welcome to SpliceAI Variant Scoring" "{ENTER}" { -not (Find-Window $termsTitle) } | Out-Null
$main = Wait-For { Find-Window '^SpliceAI Variant Scoring$' } 120
Start-Sleep 3

$vcfText = Get-Content C:\TestData\ankrd26-hg19.vcf -Raw
$box = $main.FindAll($TS::Descendants, $AnyCond) | Where-Object { $_.Current.ControlType.ProgrammaticName -in "ControlType.Edit", "ControlType.Document" } |
       Sort-Object { $_.Current.BoundingRectangle.Width * $_.Current.BoundingRectangle.Height } -Descending | Select-Object -First 1
$how = $null
if ($box) {
    try { $box.GetCurrentPattern([Windows.Automation.ValuePattern]::Pattern).SetValue($vcfText); $how = "automation" } catch {}
    if (-not $how) {
        Set-Clipboard -Value $vcfText; Click-El $box; Start-Sleep 1
        [System.Windows.Forms.SendKeys]::SendWait("^v"); $how = "clipboard"
    }
}
Log "  VCF pasted by $how"
$auto = Wait-For { Find-In $main 'from the VCF \(reference=hg19\)' } 15
Shot "app-vcf-pasted"; Check "app.build_selected_from_vcf (hg19)" $auto
$how = Act (Find-In $main '^Run$' "Button") "SpliceAI Variant Scoring" $null { Find-In $main 'Loading reference|Reading and normalizing|live-scored|^Done -- |^Failed' } 30
Log "  Run pressed by $how"
$result = Wait-For { Find-In $main '^Done -- ' } 300
Start-Sleep 2
Shot "app-run-result"
Check "app.run_finished" $result ($result.Current.Name)
Check "app.result_shows_0.29" (Find-In $main '^0\.29$')

Close-App

# --- 5. normal uninstall -------------------------------------------------------------------
$uninstallTitle = 'SpliceAI Variant Scoring Uninstall'
Start-Process "$app\unins000.exe"
$q = Wait-For { Windows | Where-Object { Find-In $_ 'completely remove' } | Select-Object -First 1 } 60
Shot "uninstall-confirm"; Check "uninstall.asks_first" $q
# Windows draws this box's buttons in its own language (Turkish Windows: Evet /
# Hayir) and makes No the default -- so press Yes by name, else the leftmost
# button; never plain Enter.
$yes = Find-In $q '^&?(Yes|Evet)$'
if (-not $yes -and $q) {
    $yes = $q.FindAll($TS::Descendants, $AnyCond) | Where-Object { $_.Current.ControlType.ProgrammaticName -eq "ControlType.Button" } |
           Sort-Object { $_.Current.BoundingRectangle.X } | Select-Object -First 1
}
$how = Act $yes $uninstallTitle $null { -not (Windows | Where-Object { Find-In $_ 'completely remove' }) }
Log "  confirmed by $how"
# If anything still used the program's files, the uninstaller asks to close it.
$prompt = Wait-For { Windows | Where-Object { Find-In $_ 'using files|close the applications' } | Select-Object -First 1 } 10
if ($prompt) { Shot "uninstall-close-programs"; Log "  uninstaller asked to close programs"; Send-Keys $uninstallTitle "{ENTER}" | Out-Null }
$ok = Wait-For { Windows | Where-Object { Find-In $_ 'successfully removed' } | Select-Object -First 1 } 180
Shot "uninstall-done"; Check "uninstall.success_message" $ok
if ($ok) { Act (Find-In $ok '^(OK|Tamam)$') $uninstallTitle "{ENTER}" { -not (Windows | Where-Object { Find-In $_ 'successfully removed' }) } | Out-Null }
Check "uninstall.folder_removed" (Wait-For { -not (Test-Path $app) } 60)
Check "uninstall.desktop_shortcut_removed" (-not (Test-Path "$env:USERPROFILE\Desktop\SpliceAI Variant Scoring.lnk"))

# --- report -------------------------------------------------------------------------------------
Log "All steps done, writing report"
try {
    $names = @($checks.Keys)
    $failed = @($names | Where-Object { -not $checks[$_].ok })
    $lines = @("Scenario: clickthrough -- $($names.Count - $failed.Count)/$($names.Count) checks passed")
    foreach ($nm in $names) { $lines += "{0}  {1,-48} {2}" -f $(if ($checks[$nm].ok) { "PASS" } else { "FAIL" }), $nm, $checks[$nm].detail }
    [IO.File]::WriteAllLines("$R\summary.txt", [string[]]$lines)
    Log "summary.txt written"
} catch {
    Log "Report failed: $($_.Exception.Message)"
} finally {
    [IO.File]::WriteAllText("$R\DONE", "done")
}
