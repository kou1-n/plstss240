# Fortran Source Overview

This document provides an overview of the main Fortran source files found in this
project. It is intended to help developers navigate the code base.  The file
list is not exhaustive, but it covers the most frequently used routines.

| File | Purpose and Notes |
|------|------------------|
| `addres.f` | Routines for computing element address indices when assembling the global matrix. Uses arrays `jr` and `jc` for row/column pointers. |
| `analys.f` | High level control of the analysis. Calls initialization, assembly and solving procedures. Variables `istep`, `nstep`, `load` manage loading steps. |
| `assemb.f` | Assembles element stiffness matrices into the global matrix `SK`. Makes use of `mk`, `ir`, `jr` arrays describing the sparse structure. |
| `chkary.f` | Checks that working arrays are large enough. Adjusts dimensions `m1`–`m29` and `n1`–`n45`. Returns error codes via `ierror`. |
| `constr.f` | Applies constraints (boundary conditions). Works on displacement vector `disp` and updates the global stiffness matrix. |
| `dgaussj.f`, `gaussj.f` | Dense Gauss–Jordan solvers for small systems. Mostly used for 2×2 or 3×3 block inversions. |
| `elastc.f` | Builds elastic constitutive tensors. Variables `de` and `ce` hold fourth-order tensors used in assembling element matrices. |
| `plastc.f` | Computes consistent elasto–plastic tangent tensors. Works with stress `sig(3,3)` and strain `str(3,3)` arrays. |
| `mapcrs.f` | Creates the global stiffness matrix in Compressed Row Storage format. Arrays `ia`, `ja` and `sk` contain the CRS structure. |
| `mapsky.f` | Alternative routine that constructs the skyline (profile) format matrix. Uses `jdiag` for the profile index. |
| `prepar.f` | Counts non-zero entries and sets up arrays for the PARDISO solver. Populates `ia` and `ja` with preliminary indices. |
| `rowupr.f` | Defines the `ia` indexing array required by PARDISO for symmetric storage. |
| `parsol.f`, `pars00.f` | Wrapper routines invoking the PARDISO solver. `parsol.f` performs factorization and solves for the displacement vector. |
| `pcgsol.f` | Uses the Intel MKL RCI Conjugate Gradient solver. The solution vector is `x`, residual `r` is monitored for convergence. |
| `skylin.f` | Implements a skyline solver for symmetric systems. The main arrays are `kd` (diagonal), `au` (upper), and `al` (lower) parts. |
| `postpr.f` | Post-processing: computes energies and writes stress/strain results. Uses `sigel`, `epsel` for stress and strain per element. |
| `hexa8a.f`, `quad4a.f`, `quad8a.f`, `tria3a.f`, `tria6a.f`, `triap1.f` | Element-level routines for computing stiffness matrices and internal forces for each element type. Typical variables include local coordinates `xi`, `wei` (weights), shape functions `N`, and derivatives `dNdx`. |
| `inith8.f`, `initp1.f`, `initq4.f`, `initq8.f`, `initt3.f`, `initt4.f`, `initt6.f` | Initialization of integration points and weights for the corresponding element types. |
| `initia.f` | Global initialization before assembling matrices. Sets default tolerances and opens result files. |
| `inform.f` | Reads node and element connectivity information from the CML input file. Populates arrays `coor`, `conn`. |
| `main.f` | Program entry point. Allocates work arrays `a(NDIM)` and coordinates the overall solution process by calling other routines. Common blocks such as `/iodev/` define logical unit numbers for I/O. |
| `output.f` | Writes displacement vectors and stress data to text files. Variables `disp`, `sigout` are used for output formatting. |
| `phexa8.f`, `pquad4.f`, `pquad8.f`, `ptria3.f`, `ptria6.f`, `ptrip1.f` | Elasto–plastic computations at integration points for each element. Maintain history variables `ehist` (plastic strain) and `beta` (back stress). |
| `stress.f` | Main constitutive routine for the Drucker–Prager model. Takes strain tensor `str(3,3)` and returns stress `sig(3,3)` along with the consistent tangent `ctens(3,3,3,3)` and energy densities. Key internal variables include the plastic multiplier `deltag` and equivalent plastic strain `alpeg`. |
| `update.f` | Updates the global displacement vector after solving the system. Also accumulates incremental plastic strain. |
| `zerodt.f` | Utility to reset large arrays to zero at the beginning of each time step. |

Some auxiliary routines, such as `prepar.f` and `rowupr.f`, are only needed when
PARDISO is selected as the linear solver. The skyline-based solver uses
`mapsky.f` instead. The solver choice is controlled in the `solut` input block
and reflected by variable `isolvr` in the code.

## 詳細: hexa8a.f と phexa8.f

`hexa8a.f` は 8 節点ヘキサ要素の剛性マトリクスを計算するサブルーチンです。形状関数の空間微分 `dndx_g` と接線剛性テンソル `ctensg` を用いて、数値積分により要素剛性 `ske` を構築します。積分点ごとにヤコビアンの行列式 `det` と積分重みを掛け合わせ、
\(B^T C B\) を体積要素に対して加算していきます。

`phexa8.f` は同じく 8 節点ヘキサ要素の内部力計算を行うサブルーチンで、弾塑性構成則に基づく応力更新と履歴変数の更新を行います。積分点ループ内でひずみテンソルを求め、履歴変数 `ehist` を参照しながら `stress` (または `st_gtn`) を呼び出して応力と接線剛性を計算します。得られた応力から内部力ベクトル `finte` を形成し、要素平均の応力やエネルギー量も集計します。
