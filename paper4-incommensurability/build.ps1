<#
.SYNOPSIS
    Builds the CLEAN (un-highlighted) copy of the Paper 4 manuscript suite.
.DESCRIPTION
    Same pipeline as build.ps1, but compiles each document with \def\CLEANCOPY{}
    prepended on the pdflatex command line. The preamble toggle
        \newif\ifhighlight \ifdefined\CLEANCOPY\highlightfalse\else\highlighttrue\fi
    therefore resolves to \highlightfalse, so every \chR/\chRB/\chE/\chA span
    renders as plain black text (no revision colors). Outputs carry a
    " (Clean Copy)" suffix so they sit alongside the highlighted copies that
    build.ps1 produces. build.ps1 is left untouched (it still builds the
    highlighted revision copy by default).
#>

$ErrorActionPreference = "Stop"

$ScriptFolder    = Split-Path -Parent $MyInvocation.MyCommand.Path
$ManuscriptDir   = Join-Path $ScriptFolder "manuscript"
$SupplementsDir  = Join-Path $ScriptFolder "supplements"
$OutputDir       = Join-Path $ScriptFolder "output"
$GlobalOutputDir = Join-Path (Split-Path $ScriptFolder -Parent) "Output_PDFs"
$FinalName       = "The Incommensurability Principle in Biological Transport.pdf"
$FinalSuppName   = "Supplemental Material - The Incommensurability Principle in Biological Transport.pdf"

# Inline TeX forcing the clean toggle, then inputting the real document.
$CleanMain = '\def\CLEANCOPY{}\input{main.tex}'
$CleanSupp = '\def\CLEANCOPY{}\input{supplemental.tex}'

# =========================================================================
# PHASE 1: Generate Data & Figures (Python)
# =========================================================================
Write-Host "`n== PHASE 1: Python compute pipeline ==" -ForegroundColor Cyan
Push-Location (Join-Path $ScriptFolder "scripts\main")
try {
    python -X utf8 compute.py
    if ($LASTEXITCODE -ne 0) { Write-Error "compute.py failed."; exit $LASTEXITCODE }
} finally { Pop-Location }

Push-Location (Join-Path $ScriptFolder "scripts\verification")
try {
    python -X utf8 generate_womersley_figure.py
    if ($LASTEXITCODE -ne 0) { Write-Error "Figure generation failed."; exit $LASTEXITCODE }
} finally { Pop-Location }

# =========================================================================
# PHASE 2: Interleaved Compilation (CLEAN toggle)
# =========================================================================
Write-Host "`n== PHASE 2: Interleaved compilation (CLEAN copy) ==" -ForegroundColor Cyan
$OriginalDir = Get-Location
try {
    Write-Host "-> PASS 1" -ForegroundColor Yellow
    Set-Location $ManuscriptDir;  pdflatex -interaction=nonstopmode -jobname=main $CleanMain | Out-Null
    Set-Location $SupplementsDir; pdflatex -interaction=nonstopmode -jobname=supplemental $CleanSupp | Out-Null

    Write-Host "-> BibTeX" -ForegroundColor Yellow
    Set-Location $ManuscriptDir;  bibtex main | Out-Null
    Set-Location $SupplementsDir; bibtex supplemental | Out-Null

    Write-Host "-> PASS 2" -ForegroundColor Yellow
    Set-Location $ManuscriptDir;  pdflatex -interaction=nonstopmode -jobname=main $CleanMain | Out-Null
    Set-Location $SupplementsDir; pdflatex -interaction=nonstopmode -jobname=supplemental $CleanSupp | Out-Null

    Write-Host "-> PASS 3" -ForegroundColor Yellow
    Set-Location $ManuscriptDir;  pdflatex -interaction=nonstopmode -jobname=main $CleanMain | Out-Null
    Set-Location $SupplementsDir; pdflatex -interaction=nonstopmode -jobname=supplemental $CleanSupp | Out-Null

    if (-not (Test-Path $OutputDir))       { New-Item -ItemType Directory -Path $OutputDir       | Out-Null }
    if (-not (Test-Path $GlobalOutputDir)) { New-Item -ItemType Directory -Path $GlobalOutputDir | Out-Null }

    # --- Move & Rename Main Manuscript (clean) ---
    Set-Location $ManuscriptDir
    $DestFile = Join-Path $OutputDir $FinalName
    if (Test-Path "main.pdf") {
        try {
            [System.IO.File]::Copy((Resolve-Path "main.pdf").Path, $DestFile, $true)
            Remove-Item -Path "main.pdf" -Force
            Copy-Item -Path $DestFile -Destination (Join-Path $GlobalOutputDir $FinalName) -Force -ErrorAction SilentlyContinue
            Write-Host " -> Clean Manuscript Built: $DestFile" -ForegroundColor Green
        } catch {
            Write-Host "`n ERROR: Failed to move main.pdf (file locked). Close PDF reader and retry." -ForegroundColor Red
        }
    } else { Write-Host " ERROR: main.pdf not found!" -ForegroundColor Red }

    # --- Move & Rename Supplemental Material (clean) ---
    Set-Location $SupplementsDir
    $DestSupp = Join-Path $OutputDir $FinalSuppName
    if (Test-Path "supplemental.pdf") {
        try {
            [System.IO.File]::Copy((Resolve-Path "supplemental.pdf").Path, $DestSupp, $true)
            Remove-Item -Path "supplemental.pdf" -Force
            Copy-Item -Path $DestSupp -Destination (Join-Path $GlobalOutputDir $FinalSuppName) -Force -ErrorAction SilentlyContinue
            Write-Host " -> Clean Supplemental Built: $DestSupp" -ForegroundColor Green
        } catch {
            Write-Host "`n ERROR: Failed to move supplemental.pdf (file locked). Close PDF reader and retry." -ForegroundColor Red
        }
    } else { Write-Host " ERROR: supplemental.pdf not found!" -ForegroundColor Red }

} finally { Set-Location $OriginalDir }

# =========================================================================
# PHASE 3: Cleanup
# =========================================================================
Write-Host "`n== PHASE 3: Cleanup ==" -ForegroundColor Cyan
Push-Location $ManuscriptDir
try { Remove-Item -Path "*.aux","*.bbl","*.log","*.out","*.toc","*.blg","*.fls","*.fdb_latexmk","*.synctex.gz" -Force -ErrorAction SilentlyContinue } finally { Pop-Location }
Push-Location $SupplementsDir
try { Remove-Item -Path "*.aux","*.bbl","*.log","*.out","*.toc","*.blg","*.fls","*.fdb_latexmk","*.synctex.gz" -Force -ErrorAction SilentlyContinue } finally { Pop-Location }

Write-Host "`n== CLEAN BUILD COMPLETE! Clean copies in output\ (suffix 'Clean') ==" -ForegroundColor Green
