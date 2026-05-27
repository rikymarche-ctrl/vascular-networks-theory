<#
.SYNOPSIS
    Builds the Paper 4 manuscript suite deterministically from scratch.
.DESCRIPTION
    1. Runs the Python compute pipeline to generate sensitivity variables and plots.
    2. Compiles the LaTeX document suite using an interleaved build cycle to ensure 
       bidirectional cross-document referencing is resolved.
    3. Automates a clean cleanup of all auxiliary LaTeX files (.aux, .bbl, etc.) 
       post-compilation.
#>

$ErrorActionPreference = "Stop"

$ScriptFolder  = Split-Path -Parent $MyInvocation.MyCommand.Path
$ManuscriptDir = Join-Path $ScriptFolder "manuscript"
$SupplementsDir = Join-Path $ScriptFolder "supplements"
$OutputDir     = Join-Path $ScriptFolder "output"
$GlobalOutputDir = Join-Path (Split-Path $ScriptFolder -Parent) "Output_PDFs"
$FinalName     = "The Incommensurability Principle in Biological Transport.pdf"
$FinalSuppName = "Supplemental Material - The Incommensurability Principle.pdf"

# =========================================================================
# PHASE 1: Generate Data & Figures (Python)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 1: Running Python Compute Pipeline" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

Push-Location (Join-Path $ScriptFolder "scripts\main")
try {
    python -X utf8 compute.py
    if ($LASTEXITCODE -ne 0) { Write-Error "compute.py failed."; exit $LASTEXITCODE }
} finally { Pop-Location }

Write-Host "`n-> Generating verification figures..." -ForegroundColor Yellow
Push-Location (Join-Path $ScriptFolder "scripts\verification")
try {
    python -X utf8 generate_womersley_figure.py
    if ($LASTEXITCODE -ne 0) { Write-Error "Figure generation failed."; exit $LASTEXITCODE }
} finally { Pop-Location }

# =========================================================================
# PHASE 2: Interleaved Compilation (LaTeX)
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 2: Interleaved Compilation (Main & Supplement)" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

$OriginalDir = Get-Location

try {
    # --- PASS 1: Generate initial .aux files ---
    Write-Host "-> PASS 1: Generating initial auxiliary files..." -ForegroundColor Yellow
    Set-Location $ManuscriptDir
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Set-Location $SupplementsDir
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null

    # --- BIBTEX ---
    Write-Host "-> Running BibTeX for both documents..." -ForegroundColor Yellow
    Set-Location $ManuscriptDir
    bibtex main | Out-Null
    
    Set-Location $SupplementsDir
    bibtex supplemental | Out-Null

    # --- PASS 2: Cross-references resolution ---
    Write-Host "-> PASS 2: Resolving cross-references..." -ForegroundColor Yellow
    Set-Location $ManuscriptDir
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Set-Location $SupplementsDir
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null

    # --- PASS 3: Finalizing formatting and anchors ---
    Write-Host "-> PASS 3: Finalizing formatting and anchor linking..." -ForegroundColor Yellow
    Set-Location $ManuscriptDir
    pdflatex -interaction=nonstopmode main.tex | Out-Null
    
    Set-Location $SupplementsDir
    pdflatex -interaction=nonstopmode supplemental.tex | Out-Null

    # --- Create output folders if they do not exist ---
    if (-not (Test-Path $OutputDir))       { New-Item -ItemType Directory -Path $OutputDir       | Out-Null }
    if (-not (Test-Path $GlobalOutputDir)) { New-Item -ItemType Directory -Path $GlobalOutputDir | Out-Null }

    # --- Move & Rename Main Manuscript ---
    Set-Location $ManuscriptDir
    $DestFile = Join-Path $OutputDir $FinalName
    if (Test-Path "main.pdf") {
        try {
            [System.IO.File]::Copy((Resolve-Path "main.pdf").Path, $DestFile, $true)
            Remove-Item -Path "main.pdf" -Force
            Copy-Item  -Path $DestFile -Destination (Join-Path $GlobalOutputDir $FinalName) -Force -ErrorAction SilentlyContinue
            Write-Host " -> Main Manuscript Built: $DestFile" -ForegroundColor Green
        } catch {
            Write-Host "`n ERROR: Failed to move main.pdf (file locked). Close PDF reader and retry." -ForegroundColor Red
        }
    } else {
        Write-Host " ERROR: main.pdf not found!" -ForegroundColor Red
    }

    # --- Move & Rename Supplemental Material ---
    Set-Location $SupplementsDir
    $DestSupp = Join-Path $OutputDir $FinalSuppName
    if (Test-Path "supplemental.pdf") {
        try {
            [System.IO.File]::Copy((Resolve-Path "supplemental.pdf").Path, $DestSupp, $true)
            Remove-Item -Path "supplemental.pdf" -Force
            Copy-Item -Path $DestSupp -Destination (Join-Path $GlobalOutputDir $FinalSuppName) -Force -ErrorAction SilentlyContinue
            Write-Host " -> Supplemental Material Built: $DestSupp" -ForegroundColor Green
        } catch {
            Write-Host "`n ERROR: Failed to move supplemental.pdf (file locked). Close PDF reader and retry." -ForegroundColor Red
        }
    } else {
        Write-Host " ERROR: supplemental.pdf not found!" -ForegroundColor Red
    }

} finally {
    Set-Location $OriginalDir
}

# =========================================================================
# PHASE 4: Global Cleanup of ALL Temporary Files
# =========================================================================
Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host " PHASE 4: Cleaning up ALL temporary files" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

Write-Host "-> Performing cleanup..." -ForegroundColor Yellow
Push-Location $ManuscriptDir
try {
    Remove-Item -Path "*.aux", "*.bbl", "*.log", "*.out", "*.toc", "*.blg", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue
} finally {
    Pop-Location
}

Push-Location $SupplementsDir
try {
    Remove-Item -Path "*.aux", "*.bbl", "*.log", "*.out", "*.toc", "*.blg", "*.fls", "*.fdb_latexmk", "*.synctex.gz" -Force -ErrorAction SilentlyContinue
} finally {
    Pop-Location
}

Write-Host "`n============================================================" -ForegroundColor Green
Write-Host " BUILD COMPLETE! All PDFs are in output\ and mirrored in Output_PDFs\" -ForegroundColor Green
Write-Host "============================================================" -ForegroundColor Green
