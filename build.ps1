<#
.SYNOPSIS
    Pipeline di build completa per Branching Networks Papers.

.DESCRIPTION
    1. Rigenera dynamic_variables.tex per tutti i paper (Python)
    2. Compila i PDF LaTeX
    3. Copia i PDF finali in Output_PDFs/
    4. Pulisce i file temporanei LaTeX

.PARAMETER Target
    "all"    — tutto (default)
    "compute"— solo step Python (rigenera variabili)
    "latex"  — solo step LaTeX (usa variabili esistenti)
    "paper1" / "paper2" / "paper3" — singolo paper
    "clean"  — rimuove file temporanei da tutte le cartelle

.EXAMPLE
    .\build.ps1
    .\build.ps1 -Target paper3
    .\build.ps1 -Target compute
    .\build.ps1 -Target clean
#>

param(
    [string]$Target = "all"
)

$ErrorActionPreference = "Stop"

$rootDir = $PSScriptRoot
$scriptsDir = Join-Path $rootDir "shared\scripts"
$outDir = Join-Path $rootDir "Output_PDFs"

# ─────────────────────────────────────────────────────────────
# Definizione paper
# ─────────────────────────────────────────────────────────────
$papers = @(
    @{
        Id        = "paper1"
        Dir       = "paper1-murray\manuscript"
        MainFile  = "main"
        FinalName = "Beyond Murray's Law.pdf"
        HasBibtex = $true
    },
    @{
        Id        = "paper2"
        Dir       = "paper2-variational\manuscript"
        MainFile  = "main"
        FinalName = "A Unified Variational Principle for Branching Transport Networks.pdf"
        HasBibtex = $true
    },
    @{
        Id        = "supp2"
        Dir       = "paper2-variational\supplements"
        MainFile  = "supplemental"
        FinalName = "Supplemental Material - Unified Variational Principle.pdf"
        HasBibtex = $false
    },
    @{
        Id        = "paper3"
        Dir       = "paper3-neural\manuscript"
        MainFile  = "main"
        FinalName = "Universal Energy Optimization for Dendritic Morphology.pdf"
        HasBibtex = $false
    }
)

# ─────────────────────────────────────────────────────────────
# Utility: pulizia file temporanei LaTeX
# ─────────────────────────────────────────────────────────────
function Remove-TexTempFiles($dir) {
    $exts = @("*.aux", "*.log", "*.out", "*.toc", "*.bbl", "*.blg",
        "*.fls", "*.fdb_latexmk", "*.synctex.gz", "texput.log")
    foreach ($ext in $exts) {
        Get-ChildItem -Path $dir -Filter $ext -ErrorAction SilentlyContinue |
        Remove-Item -Force
    }
}

# ─────────────────────────────────────────────────────────────
# Target: clean
# ─────────────────────────────────────────────────────────────
if ($Target -eq "clean") {
    Write-Host "Pulizia file temporanei..." -ForegroundColor Cyan
    foreach ($p in $papers) {
        $d = Join-Path $rootDir $p.Dir
        Remove-TexTempFiles $d
        Write-Host "  [OK] $($p.Dir)" -ForegroundColor DarkGray
    }
    Write-Host "Pulizia completata." -ForegroundColor Green
    exit 0
}

# ─────────────────────────────────────────────────────────────
# Step 1: Python — rigenera dynamic_variables.tex
# ─────────────────────────────────────────────────────────────
if ($Target -eq "all" -or $Target -eq "compute") {
    Write-Host "`n[STEP 1/2] Rigenerazione variabili dinamiche (Python)" -ForegroundColor Cyan

    $runAll = Join-Path $scriptsDir "run_all.py"
    if (-not (Test-Path $runAll)) {
        Write-Error "Non trovo $runAll"
        exit 1
    }

    Push-Location $scriptsDir
    try {
        $result = python run_all.py 2>&1
        Write-Host $result
        if ($LASTEXITCODE -ne 0) {
            Write-Error "run_all.py ha fallito (exit code $LASTEXITCODE)"
            exit 1
        }
        Write-Host "  [OK] dynamic_variables.tex aggiornati per tutti i paper" -ForegroundColor Green
    }
    finally {
        Pop-Location
    }

    if ($Target -eq "compute") { exit 0 }
}

# ─────────────────────────────────────────────────────────────
# Step 2: LaTeX — compila i PDF
# ─────────────────────────────────────────────────────────────
if ($Target -eq "all" -or $Target -eq "latex" -or ($Target -in ($papers | ForEach-Object { $_.Id }))) {
    Write-Host "`n[STEP 2/2] Compilazione LaTeX" -ForegroundColor Cyan

    if (-not (Test-Path $outDir)) {
        New-Item -ItemType Directory -Path $outDir | Out-Null
    }

    $errors = @()

    foreach ($p in $papers) {
        # Filtra per target specifico
        if ($Target -ne "all" -and $Target -ne "latex" -and $Target -ne $p.Id) {
            continue
        }

        $workDir = Join-Path $rootDir $p.Dir
        $main = $p.MainFile

        Write-Host "`n  >>> $($p.Id): $($p.FinalName)" -ForegroundColor Yellow

        if (-not (Test-Path $workDir)) {
            Write-Host "  [SKIP] Cartella non trovata: $workDir" -ForegroundColor DarkYellow
            continue
        }

        Push-Location $workDir
        try {
            # pdflatex scrive direttamente col nome finale (-jobname)
            # → nessuna rinomina, nessun conflitto con file aperti nel viewer
            $jobname = $p.FinalName -replace '\.pdf$', ''

            Write-Host "      pdflatex pass 1..." -ForegroundColor DarkGray
            pdflatex -interaction=nonstopmode -jobname "$jobname" "$main.tex" | Out-Null

            if ($p.HasBibtex) {
                Write-Host "      bibtex..." -ForegroundColor DarkGray
                bibtex "$jobname" | Out-Null
                Write-Host "      pdflatex pass 2..." -ForegroundColor DarkGray
                pdflatex -interaction=nonstopmode -jobname "$jobname" "$main.tex" | Out-Null
            }

            Write-Host "      pdflatex pass finale..." -ForegroundColor DarkGray
            pdflatex -interaction=nonstopmode -jobname "$jobname" "$main.tex" | Out-Null

            # Verifica errori fatali nel log
            $logFile = "$jobname.log"
            $fatalErrors = @()
            if (Test-Path $logFile) {
                $fatalErrors = @(Select-String -Path $logFile -Pattern "^!")
            }
            if ($fatalErrors.Count -gt 0) {
                Write-Host "  [ERRORE] $($p.Id): $($fatalErrors.Count) errori fatali nel log" -ForegroundColor Red
                $fatalErrors | Select-Object -First 5 | ForEach-Object {
                    Write-Host "    $($_.Line)" -ForegroundColor DarkRed
                }
                $errors += $p.Id
            }

            # Copia in Output_PDFs
            $pdf = "$jobname.pdf"
            $dest = Join-Path $outDir $p.FinalName
            if (Test-Path $pdf) {
                Copy-Item -Path $pdf -Destination $dest -Force
                Write-Host "  [OK] $($p.FinalName)" -ForegroundColor Green
            }
            else {
                Write-Host "  [ERRORE] PDF non generato per $($p.Id)" -ForegroundColor Red
                $errors += $p.Id
            }

            # Pulizia temporanei (i .aux/.log hanno il jobname, non "main")
            Remove-TexTempFiles $workDir

        }
        catch {
            Write-Host "  [ERRORE] $($p.Id): $_" -ForegroundColor Red
            $errors += $p.Id
        }
        finally {
            Pop-Location
        }
    }

    # ─── Riepilogo finale ───
    Write-Host "`n════════════════════════════════════" -ForegroundColor Cyan
    if ($errors.Count -eq 0) {
        Write-Host "BUILD COMPLETATA SENZA ERRORI" -ForegroundColor Green
        Write-Host "PDF pronti in: Output_PDFs\" -ForegroundColor Green
    }
    else {
        Write-Host "BUILD COMPLETATA CON ERRORI in: $($errors -join ', ')" -ForegroundColor Red
        Write-Host "Controlla i file .log nelle cartelle sorgente" -ForegroundColor Red
        exit 1
    }
    Write-Host "════════════════════════════════════" -ForegroundColor Cyan
}
