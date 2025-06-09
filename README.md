# plstss240

Edited Fortran program implementing the Drucker--Prager plasticity model with complex hardening, based on the PLSTss program family.

## Overview

This project contains a finite element analysis solver written in Fortran. It solves elasto-plastic problems using the Drucker--Prager constitutive model with advanced hardening options. The code is derived from earlier PLSTss releases and is intended for research purposes.

## Features

### Plastic State Monitoring

The program automatically monitors plastic deformation onset at the step level. When plasticity first occurs during loading, a message is written to the NOR_*.txt file:

```
First plastic deformation at step:   4
```

This feature provides clear indication of when the material transitions from elastic to plastic behavior, improving analysis workflow and result interpretation.

## Build requirements

* **Intel Fortran Compiler** (`ifort`)
* **Intel Math Kernel Library (MKL)** – set the `MKLROOT` environment variable to your installation path
## Building

The repository provides a Makefile configured for the Intel compiler toolchain. Simply run:

```bash
make
```

This compiles the source files and links against MKL, producing the executable defined by the `TARGET` variable in the Makefile (by default `$(HOME)/bin/plstss2393_nagasaku`). Object files are written to the `obj/` directory to keep the project root tidy. Use `make clean` to remove the built objects and executable.

## Running

There are two ways to run the plstss program after compilation:

### Interactive mode

Execute the produced binary. The program runs interactively and prompts for the base name of an input dataset in CML format. For example:

```bash
$ $(HOME)/bin/plstss2393_nagasaku
*) Enter the input file name without extension:  test
```

The solver reads `test.cml` and writes several result files such as `RES_test.cml` and `STS_test.txt` in the working directory.

### Non-interactive mode (batch processing)

For automated processing or to redirect output to a file, you can use input/output redirection. First, create an input file (without extension) containing the dataset name:

```bash
# Create input file containing the dataset name
$ echo "test" > inp
```

Then run the program with input/output redirection:

```bash
$ plstss2393_nagasaku < inp > out
```

Where:
- `inp` is the input file containing the dataset name (without extension)
- `out` is the output file where all program messages will be written
- `plstss2393_nagasaku` should be replaced with your actual executable name

This method is particularly useful for:
- Batch processing multiple datasets
- Running calculations in background
- Capturing all program output for later analysis
- Automated workflows and scripting

## Repository layout

To keep the different source types organized, the repository is structured as
follows:

- `scripts/` – Python utilities for post‑processing and visualisation
- `data/` – example input and output files (`.cml`, `.NEU`, `.txt`, etc.)
- Fortran source files remain in the project root so that the existing
  `Makefile` continues to work without modification.

## Source files

A detailed description of the Fortran source code can be found in
[`docs/FORTRAN_FILES.md`](docs/FORTRAN_FILES.md).  It lists the main subroutines
and typical variables used throughout the program.

All `.f` files use fixed-form Fortran. Executable code should stay within the
first 72 columns for portability, although comments may extend beyond that
limit. Use continuation lines whenever your logic exceeds the allowed width.

## Tensor representation

The program performs all constitutive calculations using `3×3` matrix forms for
stress and strain tensors. When results are saved, each matrix is converted to a
six-component array following Voigt notation. See
[`docs/tensor_representation.md`](docs/tensor_representation.md) for the mapping
used during output.

## Additional documentation

The `docs` directory contains several other references:

- [`docs/AGENT.md`](docs/AGENT.md) — tips for reliably locating files with the
  Codex agent.
- [`docs/CMLformat07.md`](docs/CMLformat07.md) — detailed structure of the CML
  input and output format.
- [`docs/plastic_monitoring.md`](docs/plastic_monitoring.md) — plastic state detection and monitoring features.
- [`docs/stress_variables.md`](docs/stress_variables.md) — variable lists for
  `stress.f` and `stress_dp.f`.
- [`docs/plstss_flow_ja.md`](docs/plstss_flow_ja.md) — PLSTss の全体フローを
  日本語で解説した資料。

## Interpretation

- **Steps 1-3**: Elastic behavior (converges in 1 iteration)
- **Step 4**: First plastic deformation (requires 4 iterations)
- **Step 5+**: Continued plastic behavior

## Advantages of Step-Level Monitoring

### Before (Integration Point Level)
- Up to 64 redundant messages (8 elements × 8 Gauss points)
- Difficult to determine onset step
- Cluttered output files

### After (Step Level)
- Single, clear message per analysis
- Precise step identification
- Clean, readable output

## Technical Details

### Modified Files

1. **`stress_vm.f`**: Removed per-integration-point output
2. **`stress_dp.f`**: Removed per-integration-point output  
3. **`output.f`**: Added step-level detection logic

### Detection Threshold

The threshold `1.d-12` is used to account for numerical precision in floating-point calculations while detecting meaningful plastic deformation.

### Energy-Based Detection

Using plastic energy (`tene_p`) for detection provides:
- Global perspective across all elements
- Reliable indication of actual plastic work
- Immunity to local numerical artifacts

## Usage

No additional user input is required. The plastic monitoring is automatically active and will output to the standard NOR_*.txt file generated during analysis.

## Development Environment

This repository has been successfully cloned and tested on the following systems:
- Original development environment
- Akita University Laboratory PC (akita)

## Remote repository

This project's remote repository is hosted at [https://github.com/kou1-n/plstss240](https://github.com/kou1-n/plstss240).

