# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PLSTss is a Fortran-based finite element analysis solver for elasto-plastic problems using the Drucker-Prager constitutive model with advanced hardening options. The codebase consists primarily of fixed-form Fortran 77 files (.f) that implement various FEM routines. The project recently implemented the Block Newton method for improved convergence in plasticity problems.

## Directory Structure

```
WORK/
├── src/                 # Fortran source files
├── obj/                 # Object files (generated during build)
├── input_files/         # Sample CML input files
├── output/              # Analysis output files
├── docs/                # Technical documentation
├── tests/               # Test scripts
└── Makefile            # Build configuration
```

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
~/bin/plstss2393_nagasaku
# Enter input filename without extension when prompted
```

**Batch mode (for automation):**
```bash
echo "inputfile" | ~/bin/plstss2393_nagasaku
```

**With output redirection:**
```bash
echo "inputfile" | ~/bin/plstss2393_nagasaku > output.log 2>&1
```

### Input/Output Files

Input files use `.cml` format (CML - Custom Markup Language for FEM). The program reads `[name].cml` and generates output files:
- `RES_[name].cml` - Results in CML format
- `STS_[name].txt` - Stress output
- `DIS_[name].txt` - Displacement output
- `ENE_[name].txt` - Energy output
- `NOR_[name].txt` - Convergence and iteration info (important for debugging)
- `TMP_[name].txt` - Temporary data

### Testing

```bash
python tests/test_hardening.py   # Test hardening functions
```

## Code Architecture

### Core Components

1. **Main Control Flow** (`main.f`, `analys.f`)
   - Entry point and high-level analysis control
   - Manages loading steps and Newton-Raphson iterations
   - Controls convergence criteria (equilibrium and yield conditions)
   - Implements Block Newton method for MATYPE=5 or 6

2. **Element Routines**
   - **Elastic stiffness**: `hexa8a.f`, `quad4a.f`, `tria3a.f`, etc.
   - **Plastic integration**: `phexa8.f`, `pquad4.f`, `ptria3.f`, etc.
   - Each element type has initialization routines (`inith8.f`, `initq4.f`, etc.)
   - 8-node hexahedral elements are primary 3D elements

3. **Constitutive Models**
   - `stress.f` - Original Drucker-Prager implementation
   - `stress_dp_bn.f` - Block Newton method for Drucker-Prager
   - `stress_dp_rm.f` - Return Mapping method
   - `stress_vm.f` - Von Mises model
   - `stress_vm_bn.f` - Block Newton for Von Mises
   - `plastc.f` - Consistent tangent computation
   - `hardfunc.f` - Hardening function implementations

4. **Solvers**
   - PARDISO solver integration (`parsol.f`, `pars00.f`, `pars99.f`)
     - `mtype=2` for symmetric positive definite
     - `mtype=11` for general unsymmetric (experimental)
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
   - `neu2cml` - External converter from FEMAP neutral format to CML

### Key Data Structures

- **Stress/Strain**: 3×3 matrix form internally, converted to 6-component Voigt notation for output
- **History Variables**:
  - `ehist(20,ngaus,nelx)` - stores plastic strain and state variables
  - `histi(50,ngaus,nelx)` - iteration history for Block Newton
- **Sparse Matrix**: CRS format with `ia`, `ja`, `sk` arrays for PARDISO
- **Convergence**: Dual criteria for equilibrium (`||R_f||`) and yield (`||R_g||`)

### Important Variables

- `istep`, `nstep` - Current and total load steps
- `sig(3,3)`, `str(3,3)` - Stress and strain tensors
- `deltag` - Plastic multiplier (Δγ)
- `alpeg` - Equivalent plastic strain
- `ctens(3,3,3,3)` - Consistent tangent tensor
- `ftreg` - Yield function value
- `g_norm` - L2 norm of yield residuals (Block Newton)

## Material Properties in CML Files

Material type (`/MATER/` section, 2nd parameter):
- `1` - Linear elastic
- `2` - Elastic-plastic (plane stress)
- `3` - Anisotropic elastic
- `4` - Drucker-Prager elastic-plastic (most common)

For Drucker-Prager (type 4), key parameters:
- Line 1: E (Young's modulus), ν (Poisson's ratio)
- Line 2: σ_y (yield stress, position 10), hardening params
- Line 3: φ (friction angle, position 16), ψ (dilation angle, position 17)

## Sample Input Files

Located in `input_files/`:
- `1elem_f.cml` - Basic single element test
- `1elem_f_uniconp.cml` - Uniaxial compression
- `1elem_f_unistress_y.cml` - Uniaxial stress in y-direction
- `1elem_phi10psi5.cml` - Non-associated flow (φ=10°, ψ=5°)

## Development Notes

### Fortran Conventions
- Fixed-form Fortran 77 (code must stay within columns 1-72)
- Continuation lines use character in column 6
- Comment lines start with 'c' or 'C' in column 1
- IMPLICIT DOUBLE PRECISION (a-h,o-z) is standard

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
- Primary files: `stress_dp_bn.f` (Block Newton), `stress_dp_rm.f` (Return Mapping)
- Hardening functions: `hardfunc.f`
- Material tangent: `plastc.f`

**Debugging convergence issues:**
- Check `NOR_[name].txt` for iteration history
- Monitor `ftreg` values (yield function)
- Verify material parameters in CML file
- Consider reducing load steps or tolerance

### Recent Development (2024-2025)

1. **Block Newton Method Implementation** (September-October 2024)
   - Implemented Yamamoto et al. (2021) algorithm in `stress_dp_bn.f`
   - Simultaneous solution of equilibrium and yield conditions
   - No local iterations required (global convergence only)
   - Improved convergence for non-associated plasticity

2. **Directory Reorganization** (September 2024)
   - Source files moved to `src/` directory
   - Input files organized in `input_files/`
   - Documentation collected in `docs/`

3. **Convergence Improvements**
   - Added damping for non-associated flow cases
   - Enhanced debug output for plastic state
   - Fixed consistent tangent formulation

### Known Issues and Limitations

- Non-associated flow (φ ≠ ψ) may require careful parameter tuning
- Large friction angles (φ > 30°) can cause convergence difficulties
- PARDISO solver limited to symmetric matrices (mtype=2) for stability

### References

- Yamamoto et al. (2021) "Simultaneously iterative procedure based on block Newton method for elastoplastic problems", Int. J. Numer. Methods Eng. 122(9), 2145-2178
- de Souza Neto et al. (2008) "Computational Methods for Plasticity"