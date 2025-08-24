# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PLSTss is a Fortran-based finite element analysis solver for elasto-plastic problems using the Drucker-Prager constitutive model with advanced hardening options. The codebase consists primarily of fixed-form Fortran 77 files (.f) that implement various FEM routines.

## Build and Run Commands

### Building the Project
```bash
make              # Builds the executable using Intel Fortran compiler
make clean        # Removes object files and executable
```

The executable is created at `$(HOME)/bin/plstss2393_nagasaku` by default. Object files are placed in the `obj/` directory.

### Running the Analysis

**Interactive mode:**
```bash
$(HOME)/bin/plstss2393_nagasaku
# Enter input filename without extension when prompted
```

**Batch mode (for automation):**
```bash
echo "test" > inp
plstss2393_nagasaku < inp > out
```

Input files use `.cml` format. The program reads `[name].cml` and generates output files:
- `RES_[name].cml` - Results in CML format
- `STS_[name].txt` - Stress output
- `DIS_[name].txt` - Displacement output
- `ENE_[name].txt` - Energy output
- `NOR_[name].txt` - Normal output with convergence info
- `TMP_[name].txt` - Temperature/temporary data

### Testing
```bash
python tests/test_hardening.py   # Test hardening functions
```

## Code Architecture

### Core Components

1. **Main Control Flow** (`main.f`, `analys.f`)
   - Entry point and high-level analysis control
   - Manages loading steps and convergence iterations

2. **Element Routines**
   - **Elastic stiffness**: `hexa8a.f`, `quad4a.f`, `tria3a.f`, etc.
   - **Plastic integration**: `phexa8.f`, `pquad4.f`, `ptria3.f`, etc.
   - Each element type has initialization routines (`inith8.f`, `initq4.f`, etc.)

3. **Constitutive Models** 
   - `stress.f` - Main Drucker-Prager implementation
   - `stress_dp_bn.f`, `stress_dp_rm.f` - Alternative return mapping schemes
   - `stress_vm.f` - Von Mises model
   - `plastc.f` - Consistent tangent computation
   - `hardfunc.f` - Hardening function implementations

4. **Solvers**
   - PARDISO solver integration (`parsol.f`, `pars00.f`, `pars99.f`)
   - Skyline solver (`skylin.f`, `mapsky.f`)
   - PCG solver (`pcgsol.f`)
   - Dense solvers for small systems (`gaussj.f`, `dgaussj.f`)

5. **Assembly and Matrix Operations**
   - `assemb.f` - Global matrix assembly
   - `mapcrs.f` - Compressed row storage format
   - `addres.f` - Element addressing

6. **I/O and Post-processing**
   - `redcml.f` - CML file reader
   - `output.f` - Results output (includes plastic state monitoring)
   - `postpr.f` - Post-processing computations

### Key Data Structures

- **Stress/Strain**: 3×3 matrix form internally, converted to 6-component Voigt notation for output
- **History Variables**: `ehist` (plastic strain), `beta` (back stress) stored at integration points
- **Sparse Matrix**: CRS format with `ia`, `ja`, `sk` arrays for PARDISO

### Important Variables

- `istep`, `nstep` - Current and total load steps
- `sig(3,3)`, `str(3,3)` - Stress and strain tensors
- `deltag` - Plastic multiplier
- `alpeg` - Equivalent plastic strain
- `ctens(3,3,3,3)` - Consistent tangent tensor

## Development Notes

### Fortran Conventions
- Fixed-form Fortran 77 (code must stay within columns 1-72)
- Continuation lines use character in column 6
- Comment lines start with 'c' or 'C' in column 1

### Build Requirements
- Intel Fortran Compiler (`ifort`)
- Intel MKL library (set `MKLROOT` environment variable)
- OpenMP support for parallel execution

### Key Files for Different Tasks

**Adding new element types:**
- Create elastic routine (like `hexa8a.f`)
- Create plastic routine (like `phexa8.f`)
- Create initialization routine (like `inith8.f`)
- Update Makefile to include new object files

**Modifying constitutive model:**
- Primary files: `stress.f`, `stress_dp_bn.f`, `stress_dp_rm.f`
- Hardening functions: `hardfunc.f`

**Changing solver:**
- Check `analys.f` for solver selection logic
- Variable `isolvr` controls solver choice

### Recent Changes
The codebase recently transitioned from `stress_dp.f` to `stress_dp_bn.f` naming convention. The Git history shows active development on return mapping algorithms and hardening models.