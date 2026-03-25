<#
.SYNOPSIS
    Builds the Paper 1 manuscript deterministically from scratch.
.DESCRIPTION
    1. Sets the PYTHONPATH to include the shared/scripts folder.
    2. Executes compute.py to generate dynamic_variables.tex and figures.
    3. Compiles the LaTeX document (pdfLaTeX + BibTeX).
#>

$ErrorActionPreference = "Stop"

# 1. Setup absolute paths relative to the script location
$ScriptFolder = Split-Path -Parent $MyInvocation.MyCommand.Path
$SharedScripts = Join-Path (Split-Path $ScriptFolder -Parent) "shared\scripts"
$ComputeScript = Join-Path $ScriptFolder "scripts\compute.py"
$ManuscriptDir = Join-Path $ScriptFolder "manuscript"

# =========================================================================
# PHASE 1: Generate Data & Figures (Python)
# =========================================================================
Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " PHASE 1: Running Python Compute Pipeline" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$env:PYTHONPATH = $SharedScripts
python $ComputeScript

if ($LASTEXITCODE -ne 0) {
    Write-Error "L'esecuzione dello script Python è fallita. Compilazione interrotta."
    exit $LASTEXITCODE
}

# =========================================================================
# PHASE 2: Compile Manuscript (LaTeX)
# =========================================================================
Write-Host ""
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " PHASE 2: Compiling LaTeX Manuscript (pdflatex + bibtex)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$OriginalDir = Get-Location
Set-Location $ManuscriptDir

try {
    Write-Host "`n-> Running pdflatex (Pass 1/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Write-Host "-> Running bibtex..." -ForegroundColor Yellow
    bibtex main | Out-Null
    
    Write-Host "-> Running pdflatex (Pass 2/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Write-Host "-> Running pdflatex (Pass 3/3)..." -ForegroundColor Yellow
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Write-Host "-> Cleaning up auxiliary files (*.log, *.aux, etc.)..." -ForegroundColor Yellow
    Remove-Item -Path "*.aux", "*.log", "*.out", "*.bbl", "*.blg", "*.toc", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue
    
    $OutputDir = Join-Path $ScriptFolder "output"
    if (-not (Test-Path $OutputDir)) {
        New-Item -ItemType Directory -Path $OutputDir | Out-Null
    }
    
    Write-Host "-> Moving and renaming output to output\Beyond Murray's Law.pdf..." -ForegroundColor Yellow
    $DestFile = Join-Path $OutputDir "Beyond Murray's Law.pdf"
    if (Test-Path "main.pdf") {
        try {
            Move-Item -Path "main.pdf" -Destination $DestFile -Force -ErrorAction Stop
            
            $GlobalOutputDir = Join-Path (Split-Path $ScriptFolder -Parent) "Output_PDFs"
            if (-not (Test-Path $GlobalOutputDir)) {
                New-Item -ItemType Directory -Path $GlobalOutputDir | Out-Null
            }
            Copy-Item -Path $DestFile -Destination (Join-Path $GlobalOutputDir "Beyond Murray's Law.pdf") -Force -ErrorAction SilentlyContinue

            Write-Host "`n============================================================" -ForegroundColor Green
            Write-Host " BUILD COMPLETE! Il PDF è aggiornato in: " -ForegroundColor Green
            Write-Host " $DestFile" -ForegroundColor Green
            Write-Host " E specchiato nella cartella master: $GlobalOutputDir" -ForegroundColor Green
            Write-Host "============================================================" -ForegroundColor Green
        } catch {
            Write-Host "`n============================================================" -ForegroundColor Red
            Write-Host " ERRORE DI SPOSTAMENTO: Il file 'main.pdf' è stato compilato," -ForegroundColor Red
            Write-Host " ma non può essere rinominato in '$DestFile'." -ForegroundColor Red
            Write-Host " Probabilmente c'è un PDF reader (es. Acrobat) aperto su uno dei due file." -ForegroundColor Red
            Write-Host " Chiudi il file e riprova!" -ForegroundColor Red
            Write-Host "============================================================" -ForegroundColor Red
        }
    } else {
        Write-Host " ERRORE: main.pdf non trovato! La compilazione LaTeX è fallita." -ForegroundColor Red
    }
}
catch {
    Write-Error "Errore critico durante lo script: $_"
}
finally {
    # Torna sempre alla cartella di partenza
    Set-Location $OriginalDir
}
