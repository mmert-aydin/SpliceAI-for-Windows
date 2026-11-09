<#
.SYNOPSIS
  Fills installer\thirdparty with everything Setup ships besides the app:
  VC++ runtime DLLs, Temurin 21 JRE, SnpEff + GRCh37.p13/GRCh38.p14 RefSeq
  databases, NCBI MANE summary, UCSC's curated RefSeq transcript tables (hg19,
  hg38). Only needed on a fresh checkout (the folder
  is git-ignored). SpliceAI is NOT fetched here -- Setup downloads it on the
  user's PC.

  Versions are pinned where upstream allows it; see BUNDLED-VERSIONS.md for
  what the current Setup was built with. SnpEff is only published as
  "latest", so compare its version with BUNDLED-VERSIONS.md after fetching.
#>
param([switch]$Force)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$third = Join-Path $PSScriptRoot "thirdparty"
$tmp = Join-Path ([IO.Path]::GetTempPath()) "spliceai-thirdparty"
New-Item -ItemType Directory -Force $third, $tmp | Out-Null

$JavaRelease = "jdk-21.0.12.1+1"
$ManeUrl = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz"
$ManeSha256 = "D10ACE2720681A3B2E0EEFD9DA4F551274A6B4141AC9BFD6A2565DFB6E9AD55C"
# Exons, CDS and every curated NM_/NR_ transcript, for the reference score
# database and the Sequence-to-cDNA tab (spliceai_pipeline/transcripts.py).
# UCSC only serves the current table: if it changes, the hash check below
# stops the fetch -- compare, then update the pin here and in BUNDLED-VERSIONS.md.
$RefSeqTables = [ordered]@{
    "hg19.ncbiRefSeqCurated.txt.gz" = "DAF41B6C1A7465F71E0DE1E680B973249B9345BAB6ADAD15B48ACD620EA30062"
    "hg38.ncbiRefSeqCurated.txt.gz" = "67090326B606658BDABDFE5F8DCB97B7943AC802C7945399522BE7174EC021CA"
}
$SnpEffUrl = "https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"
$SnpEffDbs = "GRCh37.p13", "GRCh38.p14"
$VcDlls = "concrt140", "msvcp140", "msvcp140_1", "msvcp140_2", "msvcp140_atomic_wait", "msvcp140_codecvt_ids",
          "vccorlib140", "vcomp140", "vcruntime140", "vcruntime140_1", "vcruntime140_threads"

function Step($message) { Write-Host "==> $message" -ForegroundColor Cyan }
function Need($path) { $Force -or -not (Test-Path (Join-Path $third $path)) }
function Assert-Sha256($file, $expected) {
    $actual = (Get-FileHash $file -Algorithm SHA256).Hash
    if ($actual -ne $expected) { Remove-Item $file; throw "$(Split-Path $file -Leaf): SHA-256 $actual, expected $expected" }
}

# VC++ runtime: copied from this PC's installed VC++ 2015-2022 x64
# redistributable (Microsoft permits app-local deployment of these files).
if (Need "vcruntime") {
    Step "VC++ runtime DLLs"
    $dst = New-Item -ItemType Directory -Force (Join-Path $third "vcruntime")
    foreach ($n in $VcDlls) { Copy-Item "$env:SystemRoot\System32\$n.dll" $dst -Force }
}

if (Need "java") {
    Step "Temurin JRE $JavaRelease"
    $api = "https://api.adoptium.net/v3/assets/release_name/eclipse/$([uri]::EscapeDataString($JavaRelease))?architecture=x64&image_type=jre&os=windows"
    $pkg = ((Invoke-RestMethod $api).binaries | Where-Object { $_.package.name -like "*.zip" } | Select-Object -First 1).package
    $zip = Join-Path $tmp $pkg.name
    Invoke-WebRequest -UseBasicParsing $pkg.link -OutFile $zip
    if ((Get-FileHash $zip -Algorithm SHA256).Hash -ne $pkg.checksum.ToUpper()) { throw "Java checksum mismatch" }
    Expand-Archive $zip (Join-Path $third "java") -Force
}

$snpEffDir = Join-Path $third "snpeff\snpEff"
if (Need "snpeff\snpEff\snpEff.jar") {
    Step "SnpEff (latest core)"
    $zip = Join-Path $tmp "snpEff_latest_core.zip"
    Invoke-WebRequest -UseBasicParsing $SnpEffUrl -OutFile $zip
    Expand-Archive $zip (Join-Path $third "snpeff") -Force
}
$java = (Get-ChildItem (Join-Path $third "java") -Recurse -Filter java.exe | Select-Object -First 1).FullName
foreach ($db in $SnpEffDbs) {
    if (Need "snpeff\snpEff\data\$db\snpEffectPredictor.bin") {
        Step "SnpEff database $db"
        Push-Location $snpEffDir
        try { & $java "-Duser.language=en" "-Duser.country=US" -jar snpEff.jar download -v $db } finally { Pop-Location }
        if ($LASTEXITCODE -ne 0) { throw "SnpEff download $db failed" }
    }
}

if (Need "mane") {
    Step "MANE summary"
    $dst = New-Item -ItemType Directory -Force (Join-Path $third "mane")
    $out = Join-Path $dst (Split-Path $ManeUrl -Leaf)
    Invoke-WebRequest -UseBasicParsing $ManeUrl -OutFile $out
    Assert-Sha256 $out $ManeSha256
}

foreach ($name in $RefSeqTables.Keys) {
    if (Need "refseq\$name") {
        Step "UCSC curated RefSeq transcripts ($name)"
        $dst = New-Item -ItemType Directory -Force (Join-Path $third "refseq")
        $out = Join-Path $dst $name
        $build = $name.Split(".")[0]
        Invoke-WebRequest -UseBasicParsing "https://hgdownload.soe.ucsc.edu/goldenPath/$build/database/ncbiRefSeqCurated.txt.gz" -OutFile $out
        Assert-Sha256 $out $RefSeqTables[$name]
    }
}

Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
Step "installer\thirdparty is complete."
