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
;   - copies the reference score databases (spliceai-refdb\hg19, \hg38 next
;     to Setup.exe) the same way, unless a valid one is already on the PC; an
;     invalid one there is renamed <build>.invalid-<date> first, never deleted.

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
; English or Turkish, English first (the program's default too); the license page and Setup's own texts follow
; (messages-en.isl, messages-tr.isl, terms-tr.txt).
ShowLanguageDialog=yes
LanguageDetectionMethod=none
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

[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl,messages-en.isl"; LicenseFile: "terms.txt"
Name: "tr"; MessagesFile: "compiler:Languages\Turkish.isl,messages-tr.isl"; LicenseFile: "terms-tr.txt"

[Tasks]
Name: "desktopicon"; Description: "{cm:TaskDesktopIcon}"; GroupDescription: "{cm:TaskShortcuts}"
; Only offered when Setup.exe has a reference-data folder beside it (the USB
; layout). Copying the genomes onto the PC is what lets the program find them
; by itself afterwards, with the drive unplugged -- see
; spliceai_gui/reference_locator.py, which looks here first.
Name: "referencedata"; Description: "{cm:TaskGenomes}"; GroupDescription: "{cm:TaskGenomesGroup}"; Check: ReferenceDataBeside
; The reference score databases (Sequence to cDNA tab), likewise only from the
; drive. Copied only where no valid database is on the PC (see NeedRefdb).
Name: "refdbhg19"; Description: "{cm:TaskRefdb19}"; GroupDescription: "{cm:TaskRefdbGroup}"; Check: RefdbBeside('hg19')
Name: "refdbhg38"; Description: "{cm:TaskRefdb38}"; GroupDescription: "{cm:TaskRefdbGroup}"; Check: RefdbBeside('hg38')

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
; After the app's own files, so NeedRefdb can ask the installed program.
Source: "{src}\spliceai-refdb\hg19\*"; DestDir: "{code:RefdbDir|hg19}"; Tasks: refdbhg19; Check: NeedRefdb('hg19'); Flags: external skipifsourcedoesntexist ignoreversion uninsneveruninstall
Source: "{src}\spliceai-refdb\hg38\*"; DestDir: "{code:RefdbDir|hg38}"; Tasks: refdbhg38; Check: NeedRefdb('hg38'); Flags: external skipifsourcedoesntexist ignoreversion uninsneveruninstall

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "{cm:RunNow,{#AppName}}"; Flags: nowait postinstall skipifsilent

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
  { NeedRefdb's answer per build, asked once: '' not yet, 'copy', 'keep'. }
  RefdbDecision19, RefdbDecision38: String;
  { Shown on the last page: what happened to each database. }
  RefdbNotes: String;
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

function RefdbSrc(Build: String): String;
begin
  Result := ExpandConstant('{src}\spliceai-refdb\') + Build;
end;

{ Where the program looks for its databases (spliceai_pipeline/refdb.py,
  default_root) and where its own Import puts them. }
function RefdbDir(Build: String): String;
begin
  Result := ExpandConstant('{%USERPROFILE}\SpliceAI_reference_data\spliceai-refdb\') + Build;
end;

function RefdbBeside(Build: String): Boolean;
begin
  Result := FileExists(RefdbSrc(Build) + '\manifest.json');
end;

{ Bytes of the files directly in Folder. }
function FolderSize(Folder: String): Int64;
var
  Rec: TFindRec;
begin
  Result := 0;
  if FindFirst(Folder + '\*', Rec) then
  try
    repeat
      if Rec.Attributes and FILE_ATTRIBUTE_DIRECTORY = 0 then
        Result := Result + (Int64(Rec.SizeHigh) shl 32) + Rec.SizeLow;
    until not FindNext(Rec);
  finally
    FindClose(Rec);
  end;
end;

{ Called by the database [Files] entries, after the program is installed.
  The program itself says whether the PC already has a usable database
  (SpliceAI-VariantScoring.exe --check-refdb: 0 valid, 2 none, 1 refused):
  a valid one is kept, a refused one renamed so nothing is lost, then copied. }
function NeedRefdb(Build: String): Boolean;
var
  Code: Integer;
  Dest, Moved, Decision: String;
begin
  if Build = 'hg19' then Decision := RefdbDecision19 else Decision := RefdbDecision38;
  if Decision = '' then
  begin
    Decision := 'copy';
    Dest := RefdbDir(Build);
    if DirExists(Dest) then
    begin
      WizardForm.StatusLabel.Caption := FmtMessage(CustomMessage('CheckingRefdb'), [Build]);
      if not Exec(ExpandConstant('{app}\{#AppExe}'), '--check-refdb "' + Dest + '" ' + Build, '',
                  SW_HIDE, ewWaitUntilTerminated, Code) then
        Code := -1;
      Log(Format('--check-refdb %s %s -> %d', [Dest, Build, Code]));
      if Code = 0 then
      begin
        Decision := 'keep';
        RefdbNotes := RefdbNotes + FmtMessage(CustomMessage('RefdbKept'), [Build]) + #13#10;
      end
      else if Code = 1 then
      begin
        Moved := Dest + '.invalid-' + GetDateTimeString('yyyymmdd', #0, #0);
        if DirExists(Moved) then
          Moved := Moved + '-' + GetDateTimeString('hhnnss', #0, #0);
        if RenameFile(Dest, Moved) then
          RefdbNotes := RefdbNotes + FmtMessage(CustomMessage('RefdbRenamed'), [Build, ExtractFileName(Moved)]) + #13#10
        else
        begin
          Decision := 'keep';
          RefdbNotes := RefdbNotes + FmtMessage(CustomMessage('RefdbNotRenamed'), [Build]) + #13#10;
        end;
      end;
    end;
    if Build = 'hg19' then RefdbDecision19 := Decision else RefdbDecision38 := Decision;
  end;
  Result := Decision = 'copy';
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
  Needed := 0;
  if WizardIsTaskSelected('referencedata') then
    Needed := ReferenceDataSize;
  { Counted in full even if a valid database turns out to be there already. }
  if WizardIsTaskSelected('refdbhg19') then
    Needed := Needed + FolderSize(RefdbSrc('hg19'));
  if WizardIsTaskSelected('refdbhg38') then
    Needed := Needed + FolderSize(RefdbSrc('hg38'));
  if (Needed = 0) or not GetSpaceOnDisk64(ExpandConstant('{%USERPROFILE}'), FreeBytes, TotalBytes) then
    Exit;
  { A gigabyte over, so Windows isn't left with nothing. }
  if FreeBytes >= Needed + 1073741824 then
    Exit;
  Result := SuppressibleTaskDialogMsgBox(CustomMessage('NoSpaceTitle'),
    FmtMessage(CustomMessage('NoSpaceText'), [IntToStr(Needed div 1073741824), IntToStr(FreeBytes div 1073741824)]),
    mbConfirmation, MB_YESNO, [CustomMessage('NoSpaceInstall'), CustomMessage('NoSpaceBack')], 0, IDNO) = IDYES;
  { The leading * keeps the other tasks (the desktop shortcut) as they are. }
  if Result then
    WizardSelectTasks('*!referencedata,*!refdbhg19,*!refdbhg38');
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
  DownloadPage := CreateDownloadPage(CustomMessage('DownloadTitle'),
    FmtMessage(CustomMessage('DownloadText'), ['{#SpliceAIVersion}']),
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

    case SuppressibleTaskDialogMsgBox(CustomMessage('DownloadFailedTitle'),
        FmtMessage(CustomMessage('DownloadFailedText'), [Msg]),
        mbError, MB_YESNO, [CustomMessage('DownloadRetry'), CustomMessage('DownloadSkip')], 0, IDNO) of
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
    SuppressibleMsgBox(CustomMessage('UnpackFailed'), mbError, MB_OK, IDOK);
end;

{ Setup run in Turkish on a PC where the program has never run: the program starts in Turkish too (its settings
  file, spliceai_gui/config.py; an existing one is never touched). }
procedure PresetLanguage;
var
  Folder: String;
begin
  Folder := ExpandConstant('{%USERPROFILE}\.spliceai_gui');
  if (ActiveLanguage <> 'tr') or FileExists(Folder + '\config.json') then
    Exit;
  ForceDirectories(Folder);
  SaveStringToFile(Folder + '\config.json', '{"language": "tr"}', False);
end;

procedure CurStepChanged(CurStep: TSetupStep);
begin
  if (CurStep = ssPostInstall) and (WheelPath <> '') then
    UnpackWheel;
  if CurStep = ssPostInstall then
    PresetLanguage;
end;

procedure CurPageChanged(CurPageID: Integer);
begin
  if (CurPageID = wpFinished) and (RefdbNotes <> '') then
    WizardForm.FinishedLabel.Caption := WizardForm.FinishedLabel.Caption + #13#10#13#10 +
      CustomMessage('RefdbNotesTitle') + #13#10 + RefdbNotes;
end;
