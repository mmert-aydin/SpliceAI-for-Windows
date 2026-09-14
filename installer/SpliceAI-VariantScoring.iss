; Inno Setup script for SpliceAI Variant Scoring.
;
; Don't compile this directly -- run installer\build.ps1, which builds the
; app and passes the defines below to ISCC.
;
; What Setup does:
;   - installs the app (PyInstaller onedir build) + SnpEff, Java, MANE and
;     the VC++ runtime DLLs, per user, into %LOCALAPPDATA%\Programs (no admin);
;   - downloads the SpliceAI wheel from PyPI (it is never shipped in this
;     Setup -- see LICENSE-THIRD-PARTY.md), checks its SHA-256, and unpacks
;     it into the app's _internal folder. A copy of the exact same wheel
;     placed next to Setup.exe is used instead of downloading (offline PCs).
;   - copies the reference genomes to ~\SpliceAI_reference_data when a
;     reference-data folder sits next to Setup.exe (the USB layout), so the
;     program finds them by itself and never asks the user for a file.

#ifndef AppVersion
  #error Run build.ps1 instead of compiling this script directly.
#endif
#ifndef AppDir
  #error AppDir (PyInstaller output) not defined -- run build.ps1.
#endif
#ifndef ThirdPartyDir
  #error ThirdPartyDir not defined -- run build.ps1.
#endif
#ifndef SpliceAIVersion
  #error SpliceAI* defines missing -- run build.ps1.
#endif

#define AppName "SpliceAI Variant Scoring"
#define AppExe "SpliceAI-VariantScoring.exe"

[Setup]
; Keep AppId fixed forever: it is how Setup recognises an existing install
; (upgrade in place) and how the uninstaller entry is identified.
AppId={{8F3C2A51-6D4E-4B7A-9C1E-2A5B7D9E4F10}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=MMA
VersionInfoVersion={#AppVersion}
; Per user, no admin rights: {autopf} is %LOCALAPPDATA%\Programs here. The
; app also writes SnpEff/MANE/Java downloads next to its exe, so its folder
; has to be writable by the user anyway.
PrivilegesRequired=lowest
DefaultDirName={autopf}\SpliceAI-VariantScoring
; The folder is fixed (no directory page): the uninstaller removes the whole
; folder, which is only safe because it's always this app-owned one.
DisableDirPage=yes
DisableProgramGroupPage=yes
DisableWelcomePage=no
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
LicenseFile=terms.txt
SetupIconFile=..\spliceai_gui\assets\logo.ico
UninstallDisplayIcon={app}\{#AppExe}
UninstallDisplayName={#AppName}
WizardStyle=modern
; Close a running copy of the app before replacing its files.
CloseApplications=yes
RestartApplications=no
OutputDir=output
OutputBaseFilename=SpliceAI-VariantScoring-Setup-{#AppVersion}
Compression=lzma2/ultra64
SolidCompression=yes
LZMAUseSeparateProcess=yes
LZMANumBlockThreads=4
; The SpliceAI package unpacked after the download (~25 MB).
ExtraDiskSpaceRequired=30000000

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Shortcuts:"
; Only offered when Setup.exe has a reference-data folder beside it (the USB
; layout). Copying the genomes onto the PC is what lets the program find them
; by itself afterwards, with the drive unplugged -- see
; spliceai_gui/reference_locator.py, which looks here first.
Name: "referencedata"; Description: "&Copy the reference genomes to this PC (about 6 GB)"; GroupDescription: "Reference genomes:"; Check: ReferenceDataBeside

[Files]
; The app itself (PyInstaller onedir; SpliceAI is excluded from it).
Source: "{#AppDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
; VC++ 2015-2022 x64 runtime, app-local (next to the exe, which Windows always
; searches first): TensorFlow needs msvcp140/msvcp140_1, which a vanilla
; Windows doesn't have, and installing vc_redist would need admin rights.
Source: "{#ThirdPartyDir}\vcruntime\*.dll"; DestDir: "{app}"; Flags: ignoreversion
; SnpEff + its GRCh37/GRCh38 RefSeq databases, the Java runtime it runs on, and
; the MANE summary -- the same layout "Download SnpEff..."/"Download MANE..."
; would create next to the exe (spliceai_pipeline/app_paths.py).
Source: "{#ThirdPartyDir}\java\*"; DestDir: "{app}\java"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#ThirdPartyDir}\snpeff\*"; DestDir: "{app}\snpeff"; Excludes: "\snpEff\examples,\snpEff\galaxy,\snpEff\.claude,*.log"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#ThirdPartyDir}\mane\*"; DestDir: "{app}\mane"; Flags: ignoreversion recursesubdirs createallsubdirs
; The reference genomes, if they are on the drive next to Setup.exe. "external"
; = not compiled into Setup (they are far too big, and they are public
; reference files, not part of this program); "onlyifdoesntexist" so
; reinstalling doesn't copy 6 GB again; "uninsneveruninstall" because they are
; the user's data and are worth keeping for a reinstall.
Source: "{src}\reference-data\*"; DestDir: "{code:ReferenceDataDir}"; Tasks: referencedata; Flags: external skipifsourcedoesntexist onlyifdoesntexist uninsneveruninstall recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "Start {#AppName} now"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Files Setup didn't copy itself: the unpacked SpliceAI package, Python's
; __pycache__ folders, SnpEff logs, anything downloaded later from the app.
Type: filesandordirs; Name: "{app}"
; SpliceAI installed by the app's own fallback (only if Setup couldn't).
; Settings (config.json next to it) are kept for a reinstall.
Type: filesandordirs; Name: "{%USERPROFILE}\.spliceai_gui\python-packages"

[Code]
const
  WheelName = '{#SpliceAIWheelName}';
  WheelUrl = '{#SpliceAIWheelUrl}';
  WheelSha256 = '{#SpliceAIWheelSha256}';
  DistInfo = 'spliceai-{#SpliceAIVersion}.dist-info';

var
  DownloadPage: TDownloadWizardPage;
  { Verified wheel to unpack after the files are installed; '' = nothing to do. }
  WheelPath: String;

function InternalDir: String;
begin
  Result := ExpandConstant('{app}\_internal');
end;

{ Where the reference genomes are copied to. The program looks here first
  (spliceai_gui/reference_locator.py, USER_REFERENCE_DIR) and it is also where
  the program's own "Download reference..." puts them, so the two agree. }
function ReferenceDataDir(Param: String): String;
begin
  Result := ExpandConstant('{%USERPROFILE}\SpliceAI_reference_data');
end;

function ReferenceDataSrc: String;
begin
  Result := ExpandConstant('{src}\reference-data');
end;

{ True when Setup.exe was started from a drive that carries the genomes -- the
  only case where there is anything to copy. }
function ReferenceDataBeside: Boolean;
begin
  Result := FileExists(ReferenceDataSrc + '\hg19.fa') or
            FileExists(ReferenceDataSrc + '\hg38.fa');
end;

function ReferenceDataSize: Int64;
var
  Rec: TFindRec;
begin
  Result := 0;
  if FindFirst(ReferenceDataSrc + '\*', Rec) then
  try
    repeat
      if Rec.Attributes and FILE_ATTRIBUTE_DIRECTORY = 0 then
        Result := Result + (Int64(Rec.SizeHigh) shl 32) + Rec.SizeLow;
    until not FindNext(Rec);
  finally
    FindClose(Rec);
  end;
end;

{ False cancels the page change. Asks before starting a multi-GB copy that
  would fill the disk -- Setup can't undo a failed one, and the genomes go to
  the user's profile, which is often on a smaller drive than it looks. }
function ReferenceDataFits: Boolean;
var
  Needed, FreeBytes, TotalBytes: Int64;
begin
  Result := True;
  if not WizardIsTaskSelected('referencedata') then
    Exit;
  Needed := ReferenceDataSize;
  if (Needed = 0) or not GetSpaceOnDisk64(ExpandConstant('{%USERPROFILE}'), FreeBytes, TotalBytes) then
    Exit;
  { A gigabyte over, so Windows isn't left with nothing. }
  if FreeBytes >= Needed + 1073741824 then
    Exit;
  Result := SuppressibleTaskDialogMsgBox('Not enough free space for the reference genomes',
    'Copying them needs about ' + IntToStr(Needed div 1073741824) + ' GB, and this drive has ' +
    IntToStr(FreeBytes div 1073741824) + ' GB free.' + #13#10#13#10 +
    'You can install without copying them and pick the files from the drive in the program ' +
    'later (it keeps working only while the drive is plugged in), or free up space and run ' +
    'Setup again.',
    mbConfirmation, MB_YESNO, ['&Install without copying', '&Go back'], 0, IDNO) = IDYES;
  { The leading * keeps the other tasks (the desktop shortcut) as they are. }
  if Result then
    WizardSelectTasks('*!referencedata');
end;

function SpliceAIPresent: Boolean;
begin
  Result := FileExists(InternalDir + '\' + DistInfo + '\METADATA') and
            FileExists(InternalDir + '\spliceai\models\spliceai5.h5') and
            FileExists(InternalDir + '\spliceai\utils.py');
end;

function OnDownloadProgress(const Url, FileName: String; const Progress, ProgressMax: Int64): Boolean;
begin
  Result := True;
end;

procedure InitializeWizard;
begin
  DownloadPage := CreateDownloadPage('Downloading SpliceAI',
    'Setup is downloading SpliceAI {#SpliceAIVersion} from pypi.org (about 16 MB).',
    @OnDownloadProgress);
  DownloadPage.ShowBaseNameInsteadOfUrl := True;
end;

{ Returns True when a verified wheel is ready in WheelPath, False when the
  user chose to continue without it. Never returns until one of those. }
function ObtainWheel: Boolean;
var
  Local, Msg: String;
begin
  Result := False;

  { A copy next to Setup.exe (e.g. on the same USB drive) wins if it is the
    exact published file. }
  Local := ExpandConstant('{src}\') + WheelName;
  if FileExists(Local) and (CompareText(GetSHA256OfFile(Local), WheelSha256) = 0) then
  begin
    Log('Using local SpliceAI wheel: ' + Local);
    WheelPath := Local;
    Result := True;
    Exit;
  end;

  while True do
  begin
    DownloadPage.Clear;
    DownloadPage.Add(WheelUrl, WheelName, WheelSha256);
    DownloadPage.Show;
    try
      try
        DownloadPage.Download;
        WheelPath := ExpandConstant('{tmp}\') + WheelName;
        Result := True;
        Exit;
      except
        Msg := GetExceptionMessage;
        Log('SpliceAI download failed: ' + Msg);
      end;
    finally
      DownloadPage.Hide;
    end;

    case SuppressibleTaskDialogMsgBox('SpliceAI could not be downloaded',
        'Setup needs an internet connection once to download SpliceAI.' + #13#10#13#10 +
        'Details: ' + Msg + #13#10#13#10 +
        'Check that this computer is connected to the internet and try again. ' +
        'If you continue without it, everything else is installed and the program ' +
        'will offer to download SpliceAI when it starts.',
        mbError, MB_YESNO, ['&Try again', '&Continue without SpliceAI'], 0, IDNO) of
      IDYES: ;  { loop and retry }
    else
      Exit;
    end;
  end;
end;

function NextButtonClick(CurPageID: Integer): Boolean;
begin
  Result := True;
  if CurPageID = wpSelectTasks then
  begin
    Result := ReferenceDataFits;
    Exit;
  end;
  if CurPageID <> wpReady then
    Exit;
  WheelPath := '';
  { Reinstall/upgrade over an install that already has it: nothing to fetch. }
  if SpliceAIPresent then
  begin
    Log('SpliceAI already present, skipping download.');
    Exit;
  end;
  ObtainWheel;
end;

function RunHidden(const Exe, Params: String): Boolean;
var
  Code: Integer;
begin
  Result := Exec(Exe, Params, '', SW_HIDE, ewWaitUntilTerminated, Code) and (Code = 0);
  Log(Format('%s %s -> ok=%d code=%d', [Exe, Params, Ord(Result), Code]));
end;

{ A wheel is a zip file. Windows 10 1803+ ships tar.exe, which unpacks zips;
  PowerShell's Expand-Archive is the fallback. Neither shows a window. }
procedure UnpackWheel;
var
  ZipCopy: String;
begin
  WizardForm.StatusLabel.Caption := 'Installing SpliceAI...';
  if not RunHidden(ExpandConstant('{sys}\tar.exe'),
      '-xf "' + WheelPath + '" -C "' + InternalDir + '"') or not SpliceAIPresent then
  begin
    ZipCopy := ExpandConstant('{tmp}\spliceai-wheel.zip');
    if FileCopy(WheelPath, ZipCopy, False) then
      RunHidden(ExpandConstant('{sys}\WindowsPowerShell\v1.0\powershell.exe'),
        '-NoProfile -NonInteractive -ExecutionPolicy Bypass -Command "Expand-Archive -Force ' +
        '-LiteralPath ''' + ZipCopy + ''' -DestinationPath ''' + InternalDir + '''"');
  end;
  if not SpliceAIPresent then
    SuppressibleMsgBox('SpliceAI was downloaded but could not be unpacked. ' +
      'The program will offer to install it when it starts.', mbError, MB_OK, IDOK);
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if (CurStep = ssPostInstall) and (WheelPath <> '') then
    UnpackWheel;
end;
