# Block Newton Method Implementation Documentation

## Overview
This document consolidates all Block Newton method related documentation for the Drucker-Prager plasticity implementation in PLSTss.

## Algorithm Summary (Yamamoto et al., 2021)

The Block Newton method solves the elastoplastic problem through simultaneous iteration of equilibrium equations and yield conditions, eliminating the need for local iterations.

### Key Features:
1. **No local iterations** - Global convergence only
2. **Simultaneous solution** of R_f (equilibrium) and R_g (yield) = 0
3. **Algebraic update** of consistency parameter Δγ
4. **Improved convergence** for non-associated plasticity

## Implementation Details

### File: `stress_dp_bn.f`

#### BOX 1: Stress Update Algorithm
- Computes trial stress (elastic predictor)
- Evaluates yield function
- Updates plastic strain and stress using Block Newton scheme
- No return mapping iterations required

#### BOX 2: Tangent Stiffness
- Computes consistent tangent: C^ep = C - N^-1 L ⊗ M
- L, M matrices derived from flow rule
- N scalar from consistency condition

### Mathematical Formulation

For Drucker-Prager yield function:
```
f = √(J_2) + η*I_1 - ξ*(σ_y + H(α))
```

Where:
- J_2 = second invariant of deviatoric stress
- I_1 = first invariant (trace)
- η, ξ = material parameters from friction angle φ
- H(α) = hardening function

### Convergence Criteria

Two residuals must converge:
1. **Equilibrium**: ||R_f|| / ||F|| < tolerance
2. **Yield condition**: ||R_g|| < tolerance

## Known Issues and Solutions

### Non-Associated Flow (φ ≠ ψ)
- May require damping for stability
- Reduced friction angle difference improves convergence
- Small load increments recommended

### Large Friction Angles (φ > 30°)
- Convergence difficulties observed
- Consider using Return Mapping (`stress_dp_rm.f`) for extreme cases

## Testing and Validation

### Test Cases in `input_files/`:
- `basic/1elem_f.cml` - Standard test
- `parameter_study/1elem_phi10psi5.cml` - Non-associated flow
- `archive/1elem_extreme_nonassoc.cml` - Extreme conditions

### Verification Procedure:
1. Check `NOR_*.txt` for iteration history
2. Monitor `ftreg` values (should approach 0)
3. Verify dual convergence (R_f and R_g)

## References

- Yamamoto, K., et al. (2021). "Simultaneously iterative procedure based on block Newton method for elastoplastic problems." International Journal for Numerical Methods in Engineering, 122(9), 2145-2178.
- de Souza Neto, E.A., et al. (2008). Computational Methods for Plasticity: Theory and Applications. Wiley.

## Development History

- **September 2024**: Initial implementation
- **October 2024**: Convergence improvements for non-associated flow
- **October 2024**: Documentation consolidation