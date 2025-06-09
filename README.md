# plstss240

Edited Fortran program implementing the Drucker--Prager plasticity model with complex hardening, based on the PLSTss program family.

## Overview

This project contains a finite element analysis solver written in Fortran. It solves elasto-plastic problems using the Drucker--Prager constitutive model with advanced hardening options. The code is derived from earlier PLSTss releases and is intended for research purposes.

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

After compilation, execute the produced binary. The program runs interactively and prompts for the base name of an input dataset in CML format. For example:

```bash
$ $(HOME)/bin/plstss2393_nagasaku
*) Enter the input file name without extension:  test
```

The solver reads `test.cml` and writes several result files such as `RES_test.cml` and `STS_test.txt` in the working directory.

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

- [`docs/stress_variables.md`](docs/stress_variables.md) — variable lists for
  `stress.f` and `stress_dp.f`.
- [`docs/plstss_flow_ja.md`](docs/plstss_flow_ja.md) — PLSTss の全体フローを
  日本語で解説した資料。


## Development Environment

This repository has been successfully cloned and tested on the following systems:
- Original development environment
- Akita University Laboratory PC (akita)

## Remote repository

This project's remote repository is hosted at [https://github.com/kou1-n/plstss240](https://github.com/kou1-n/plstss240).

