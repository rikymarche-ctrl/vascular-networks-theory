# Scripts Directory

This directory contains all computational scripts for Paper 4: "The Incommensurability Principle in Biological Transport".

## Structure

### `main/`
Primary computational pipeline scripts:
- **`compute.py`**: Main computation engine that generates `dynamic_variables.tex` and data files for manuscript figures. Implements automatic dynamic pruning to remove orphaned LaTeX macros. Run by `build.ps1` during the build process.

### `verification/`
Numerical verification scripts for supplemental material:
- **`generate_womersley_figure.py`**: Generates the critical Womersley verification figure showing Wo_c = √3 (auto-run by `build.ps1`)
- **`verify_heterogeneity.py`**: Tests robustness of minimax attractor under anatomical heterogeneities (asymmetry, taper)
- **`verify_shielding.py`**: Validates topological shielding hypothesis - elastic wall perturbations have negligible effect on alpha*

### `utility/`
Utility scripts for manuscript maintenance:
- **`wrap_project_tex.py`**: Formats all LaTeX files to 80-character line width for readability and version control
- **`check_orphaned_macros.py`**: Analyzes dynamic_variables.tex for unused macro definitions

## Usage

### Building the manuscript
The build process automatically:
1. Runs `compute.py` to generate dynamic variables and data
2. Applies dynamic pruning to remove orphaned macros
3. Generates verification figures (e.g., `generate_womersley_figure.py`)
4. Compiles LaTeX documents

Simply run from project root:
```powershell
.\build.ps1
```

### Running verification scripts
Navigate to the script directory and run:
```powershell
cd scripts\verification
python generate_womersley_figure.py
python verify_heterogeneity.py
python verify_shielding.py
```

Expected output:
- **Womersley verification**: Numerical Wo_c within 0.5% of analytical √3
- **Heterogeneity**: Asymmetry/taper shifts < 0.1 in alpha*
- **Shielding**: Elastic wall corrections |Δα*| < 0.02

### Running utility scripts
Format all LaTeX files to 80-character width:
```powershell
python scripts\utility\wrap_project_tex.py
```

Check for orphaned macro definitions:
```powershell
python scripts\utility\check_orphaned_macros.py
```

## Dynamic Pruning

`compute.py` implements automatic dynamic pruning of LaTeX macros:
1. Generates all numerical variables from physics calculations
2. Scans entire manuscript suite to count macro usage
3. Exports only macros that are actually cited (usage > 0)
4. Reports number of orphaned macros removed

This ensures `dynamic_variables.tex` stays clean and contains only used definitions.

## Notes
- All scripts use UTF-8 encoding
- Output files are generated in `manuscript/figures/` or `supplements/figures/`
- Python dependencies: numpy, scipy, matplotlib
- Verification scripts use single-harmonic approximations for computational efficiency (sufficient for validation)
