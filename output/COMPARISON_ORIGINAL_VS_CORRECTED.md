# Comparison: Original vs Corrected Drucker-Prager Formulation

## Summary
Testing non-associated flow (φ≠ψ) with LOVE CLAUDE CODE commit (906e459) vs theoretically corrected formulation.

## Test Results

### Original Implementation (LOVE CLAUDE CODE - commit 906e459)
**Formulation (INCORRECT but works):**
- H = ξ² * (∂K/∂α)  [Line 164 of stress_dp_rm.f]
- Δα = ξ * Δγ       [Line 168 of stress_dp_rm.f]

| Test Case | φ (deg) | ψ (deg) | Return Mapping | Global Newton | Final Residual | Status |
|-----------|---------|---------|----------------|---------------|----------------|--------|
| Test 1 | 10 | 5 | 2 iterations | 17 iterations | 4.0392E-10 | **CONVERGED** |
| Test 2 | 20 | 15 | 2 iterations | 16 iterations | 2.6003E-10 | **CONVERGED** |

**Observations:**
- Convergence is LINEAR (not quadratic) - takes 16-17 iterations
- Return mapping converges excellently (2 iterations)
- Despite incorrect theory, it converges due to reduced asymmetry

### Corrected Implementation
**Formulation (THEORETICALLY CORRECT):**
- H = ξ * (∂K/∂α) / √3  [Line 163 of stress_dp_rm.f]
- Δα = Δγ / √3         [Line 167 of stress_dp_rm.f]

| Test Case | φ (deg) | ψ (deg) | Return Mapping | Global Newton | Final Residual | Status |
|-----------|---------|---------|----------------|---------------|----------------|--------|
| Test 1 | 10 | 5 | 2 iterations | 20 iterations | 1.1692E-07 | **FAILED** (Error ID: 20) |
| Test 2 | 20 | 15 | 2 iterations | 20 iterations | 9.4353E-08 | **FAILED** (Error ID: 20) |

**Observations:**
- Return mapping still converges excellently (2 iterations)
- Global Newton fails after 20 iterations
- Residual gets stuck around 1E-07 to 1E-08, cannot reach tolerance (1E-10)

## Root Cause Analysis

### Matrix Storage Limitation
The fundamental issue is in the matrix assembly routines:

**mapcrs.f (Line 90):**
```fortran
if(m_jb.le.m_ia) then    ! Only stores lower triangular part
```

**rowupr.f (Line 19):**
```fortran
c Non-Zero Elements are Stored in ROW Major UPPER Triangular Format
```

The system assumes **symmetric matrices** and only stores the upper triangular part.

### Why Original Code "Works"
The incorrect formulation H = ξ² * (∂K/∂α) accidentally reduces the asymmetry in the consistent tangent matrix:
- The extra ξ factor scales the hardening contribution
- This happens to make the matrix "more symmetric"
- The symmetric solver can handle the reduced asymmetry

### Why Correct Code Fails
The theoretically correct formulation produces a truly asymmetric tangent:
- Non-associated plasticity (φ≠ψ) inherently produces D_ep ≠ D_ep^T
- The asymmetric tangent cannot be properly stored in symmetric format
- Information is lost, preventing convergence

## Log Files
- `output/logs/LOVE_CLAUDE_phi10psi5.log` - Original code, φ=10°, ψ=5° (CONVERGED)
- `output/logs/LOVE_CLAUDE_phi20psi15.log` - Original code, φ=20°, ψ=15° (CONVERGED)
- `output/logs/corrected_phi10psi5.log` - Corrected code, φ=10°, ψ=5° (FAILED)
- `output/logs/corrected_phi20psi15.log` - Corrected code, φ=20°, ψ=15° (FAILED)

## Conclusion
There is a **fundamental architectural limitation** in the FEM code:
1. The matrix storage only supports symmetric matrices
2. Non-associated plasticity requires asymmetric tangent matrices
3. The original "incorrect" implementation works by accident - it reduces asymmetry enough for the symmetric solver to handle
4. The theoretically correct implementation cannot work without major refactoring to support full asymmetric matrix storage and solving

## Recommendations
To properly implement non-associated Drucker-Prager plasticity, one of the following is needed:
1. **Major refactoring**: Modify mapcrs.f, rowupr.f, and PARDISO settings to handle full asymmetric matrices
2. **Alternative approach**: Use iterative methods that can handle asymmetric systems
3. **Approximation**: Continue using the "incorrect" formulation that happens to work
4. **Restriction**: Only use associated flow (φ=ψ) where the tangent remains symmetric