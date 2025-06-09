# Changelog

## [Unreleased]
- **Added plastic state detection**: Implemented step-level plastic deformation monitoring in NOR_*.txt output files
  - Removed redundant per-integration-point output from `stress_vm.f` and `stress_dp.f`
  - Added step-level plastic detection in `output.f` using plastic energy comparison
  - Output format: "First plastic deformation at step: X" when plasticity first occurs
  - Significantly improved readability of convergence monitoring files
- Updated `docs/stress_variables.md` to match the current `stress_dp` subroutine interface and added auto-check note.
- Added CI script `scripts/check_stress_dp_signature.py` to verify documentation against source.
