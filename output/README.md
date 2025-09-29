# Output Directory

This directory contains FEM analysis output files from PLSTss.

## Directory Structure

```
output/
├── README.md              # This file
├── extract2csv.py         # Script to extract data to CSV format
│
├── test_runs/             # Test execution logs and debugging outputs
│   └── test_*.txt
│
├── debug_runs/            # Debug output files
│   └── debug_*.txt
│
├── csv_exports/           # Exported CSV files from FEM results
│   └── RES_*_*.csv
│
├── old_results/           # Archived/older analysis results
│   └── [old FEM output files]
│
└── [Active FEM Results]   # Current analysis outputs
    ├── RES_*.cml         # Result files in CML format
    ├── STS_*.txt         # Stress-strain data
    ├── DIS_*.txt         # Displacement data
    ├── ENE_*.txt         # Energy data
    ├── NOR_*.txt         # Convergence norms
    └── TMP_*.txt         # Temperature/temporary data
```

## File Naming Convention

FEM output files follow this pattern:
- `{TYPE}_{name}.{ext}`
  - TYPE: RES, STS, DIS, ENE, NOR, TMP
  - name: Input filename without extension
  - ext: cml for results, txt for data

## Current Active Results

Latest analysis results for:
- `1elem_f` - Base element test
- `1elem_phi10psi5` - φ=10°, ψ=5° Drucker-Prager
- `1elem_phi10psi5_step11` - Step 11 yielding test
- `1elem_phi15` - φ=15° associated plasticity
- `1elem_phi20` - φ=20° associated plasticity
- `1elem_phi20psi10` - φ=20°, ψ=10° non-associated
- `1elem_voce` - Voce hardening model

## Usage

To extract data from CML files to CSV:
```bash
python extract2csv.py RES_filename.cml
```

## Cleanup Policy

- Test runs and debug outputs → `test_runs/` and `debug_runs/`
- CSV exports → `csv_exports/`
- Superseded results → `old_results/`
- Keep only latest representative analysis results in root

---
Updated: 2025-09-30