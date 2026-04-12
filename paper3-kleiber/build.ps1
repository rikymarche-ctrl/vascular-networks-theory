<#
.SYNOPSIS
    Builds the Paper 3 manuscript deterministically from scratch.
.DESCRIPTION
    1. Executes compute.py to generate dynamic_variables.tex and figures.
    2. Compiles the LaTeX document (pdfLaTeX + BibTeX).
#>

$ErrorActionPreference = "Stop"

$ScriptFolder  = Split-Path -Parent $MyInvocation.MyCommand.Path
$ComputeScript = Join-Path $ScriptFolder "scripts\compute.py"
$ManuscriptDir = Join-Path $ScriptFolder "manuscript"
$OutputDir     = Join-Path $ScriptFolder "output"
$GlobalOutputDir = Join-Path (Split-Path $ScriptFolder -Parent) "Output_PDFs"
$FinalName     = "The Dynamic Origin of Kleibers Law.pdf"

# =========================================================================
# PHASE 1: Generate Data & Figures (Python)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 1: Running Python Compute Pipeline" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

Push-Location (Join-Path $ScriptFolder "scripts")
try {
    python -X utf8 compute.py
    if ($LASTEXITCODE -ne 0) { Write-Error "compute.py failed."; exit $LASTEXITCODE }
} finally { Pop-Location }

# =========================================================================
# PHASE 2: Compile Manuscript (LaTeX)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 2: Compiling LaTeX Manuscript (pdflatex + bibtex)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$OriginalDir = Get-Location
Set-Location $ManuscriptDir

try {
    Write-Host "-> Running pdflatex (Pass 1/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null

    Write-Host "-> Running bibtex..." -ForegroundColor Yellow
    bibtex main | Out-Null

    Write-Host "-> Running pdflatex (Pass 2/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null

    Write-Host "-> Running pdflatex (Pass 3/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null

    Write-Host "-> Cleaning up auxiliary files..." -ForegroundColor Yellow
    Remove-Item -Path "*.aux", "*.log", "*.out", "*.bbl", "*.blg", "*.toc", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue

    if (-not (Test-Path $OutputDir))       { New-Item -ItemType Directory -Path $OutputDir       | Out-Null }
    if (-not (Test-Path $GlobalOutputDir)) { New-Item -ItemType Directory -Path $GlobalOutputDir | Out-Null }

    $DestFile = Join-Path $OutputDir $FinalName
    if (Test-Path "main.pdf") {
        try {
            [System.IO.File]::Copy((Resolve-Path "main.pdf").Path, $DestFile, $true)
            Remove-Item -Path "main.pdf" -Force
            Copy-Item  -Path $DestFile -Destination (Join-Path $GlobalOutputDir $FinalName) -Force -ErrorAction SilentlyContinue

            Write-Host "`n============================================================" -ForegroundColor Green
            Write-Host " BUILD COMPLETE! PDF aggiornato in:" -ForegroundColor Green
            Write-Host " $DestFile" -ForegroundColor Green
            Write-Host " E specchiato in: $GlobalOutputDir" -ForegroundColor Green
            Write-Host "============================================================" -ForegroundColor Green
        } catch {
            Write-Host "`n ERRORE DI SPOSTAMENTO: $($_.Exception.Message)" -ForegroundColor Red
        }
    } else {
        Write-Host " ERRORE: main.pdf non trovato! Compilazione LaTeX fallita." -ForegroundColor Red
    }
} catch {
    Write-Error "Errore critico: $_"
} finally {
    Set-Location $OriginalDir
}
