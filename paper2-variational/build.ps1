<#
.SYNOPSIS
    Builds the Paper 2 manuscript and supplements deterministically from scratch.
.DESCRIPTION
    1. Sets the PYTHONPATH to include the shared/scripts folder.
    2. Executes compute.py and compute_supp.py to generate dynamic variables and figures.
    3. Compiles the LaTeX documents (main and supplemental).
#>

$ErrorActionPreference = "Stop"

$ScriptFolder = Split-Path -Parent $MyInvocation.MyCommand.Path
$SharedScripts = Join-Path (Split-Path $ScriptFolder -Parent) "shared\scripts"
$ComputeScript = Join-Path $ScriptFolder "scripts\compute.py"
$ComputeSuppScript = Join-Path $ScriptFolder "scripts\compute_supp.py"
$ManuscriptDir = Join-Path $ScriptFolder "manuscript"
$SupplementsDir = Join-Path $ScriptFolder "supplements"
$OutputDir = Join-Path $ScriptFolder "output"
$GlobalOutputDir = Join-Path (Split-Path $ScriptFolder -Parent) "Output_PDFs"

# =========================================================================
# PHASE 1: Generate Data & Figures (Python)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 1: Running Python Compute Pipelines" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$env:PYTHONPATH = $SharedScripts

Write-Host "-> Running scripts\compute.py..."
python $ComputeScript
if ($LASTEXITCODE -ne 0) { Write-Error "Execution of compute.py failed."; exit $LASTEXITCODE }

Write-Host "-> Running scripts\compute_supp.py..."
python $ComputeSuppScript
if ($LASTEXITCODE -ne 0) { Write-Error "Execution of compute_supp.py failed."; exit $LASTEXITCODE }

# =========================================================================
# PHASE 2: Compile Manuscript (LaTeX)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 2: Compiling Main Manuscript (pdflatex + bibtex)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$OriginalDir = Get-Location
try {
    Set-Location $ManuscriptDir
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    bibtex main | Out-Null
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Remove-Item -Path "*.aux", "*.log", "*.out", "*.bbl", "*.blg", "*.toc", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue
    
    if (-not (Test-Path $OutputDir)) { New-Item -ItemType Directory -Path $OutputDir | Out-Null }
    if (-not (Test-Path $GlobalOutputDir)) { New-Item -ItemType Directory -Path $GlobalOutputDir | Out-Null }
    
    $DestMain = Join-Path $OutputDir "A Unified Variational Principle for Branching Transport Networks.pdf"
    if (Test-Path "main.pdf") {
        Move-Item -Path "main.pdf" -Destination $DestMain -Force -ErrorAction Stop
        Copy-Item -Path $DestMain -Destination (Join-Path $GlobalOutputDir "A Unified Variational Principle for Branching Transport Networks.pdf") -Force -ErrorAction SilentlyContinue
        Write-Host " -> Main Manuscript Built: $DestMain" -ForegroundColor Green
    } else {
        Write-Error "main.pdf not found!"
    }
} finally {
    Set-Location $OriginalDir
}

# =========================================================================
# PHASE 3: Compile Supplements (LaTeX)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 3: Compiling Supplemental Material (pdflatex)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

try {
    Set-Location $SupplementsDir
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null
    
    Remove-Item -Path "*.aux", "*.log", "*.out", "*.bbl", "*.blg", "*.toc", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue
    
    $DestSupp = Join-Path $OutputDir "Supplemental Material - Unified Variational Principle.pdf"
    if (Test-Path "supplemental.pdf") {
        Move-Item -Path "supplemental.pdf" -Destination $DestSupp -Force -ErrorAction Stop
        Copy-Item -Path $DestSupp -Destination (Join-Path $GlobalOutputDir "Supplemental Material - Unified Variational Principle.pdf") -Force -ErrorAction SilentlyContinue
        Write-Host " -> Supplemental Material Built: $DestSupp" -ForegroundColor Green
    } else {
        Write-Error "supplemental.pdf not found!"
    }
} finally {
    Set-Location $OriginalDir
}

Write-Host "`n============================================================" -ForegroundColor Green
Write-Host " BUILD COMPLETE! All PDFs are in output\ and mirrored to Output_PDFs\" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Green
