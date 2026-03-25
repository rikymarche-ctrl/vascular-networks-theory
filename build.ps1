<#
.SYNOPSIS
    Master Pipeline per Branching Networks Papers.
    Invia i comandi ai build.ps1 dedicati in ogni cartella paper per la massima modularità.
.PARAMETER Target
    "all"    — compila tutti i paper abilitati
    "paper1" / "paper2" — compila il singolo paper delegando al suo script
#>

param(
    [string]$Target = "all"
)

$ErrorActionPreference = "Stop"
$rootDir = $PSScriptRoot

# Ordine di compilazione e cartelle abilitate
$papers = @(
    "paper1-murray",
    "paper2-variational",
    "paper3-kleiber",
    "paper4-incommensurability"
    # "paper5-neural"   # work in progress — excluded from default build
)

Write-Host "============================================================" -ForegroundColor Magenta
Write-Host " BRANCHING NETWORKS MASTER BUILD PIPELINE" -ForegroundColor Magenta
Write-Host "============================================================" -ForegroundColor Magenta

if ($Target -eq "clean") {
    Write-Host "Clean command is distributed to local repos (to be implemented)." -ForegroundColor DarkYellow
    exit 0
}

$targetPapers = @()
if ($Target -eq "all") {
    $targetPapers = $papers
} else {
    $targetPapers = $papers | Where-Object { $_ -match $Target }
    if ($targetPapers.Length -eq 0) {
        Write-Error "Nessuna cartella corrisponde al target '$Target'."
        exit 1
    }
}

foreach ($p in $targetPapers) {
    Write-Host "`n>>> DELEGATING BUILD TO: $p" -ForegroundColor Magenta
    $pDir = Join-Path $rootDir $p
    if (Test-Path (Join-Path $pDir "build.ps1")) {
        Push-Location $pDir
        try {
            # Chiama lo script locale
            .\build.ps1
            if ($LASTEXITCODE -ne 0) {
                Write-Error "Build fallita nel modulo $p"
            }
        } finally {
            Pop-Location
        }
    } else {
        Write-Host "  [SKIP] Script build.ps1 non trovato in $p" -ForegroundColor DarkYellow
    }
}

Write-Host "`n============================================================" -ForegroundColor Magenta
Write-Host " MASTER PIPELINE COMPLETATA" -ForegroundColor Magenta
Write-Host " Tutti i PDF sono stati aggiornati globalmente in /Output_PDFs" -ForegroundColor Magenta
Write-Host "============================================================" -ForegroundColor Magenta
