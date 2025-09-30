# Drucker-Prager Block Newton Implementation Report

## Date: 2024-09-30
## Implementation Status: COMPLETE ✓

## Summary
Successfully implemented all fixes from `claude_dp_bn_fix_plan.md` to correct the Drucker-Prager Block Newton method in `stress_dp_bn.f`.

## Implemented Fixes

### 1. ✓ Fixed Block Newton Denominator (Lines 210-211, 257-258)
**Before:**
```fortran
N_scalar = -(2.d0*vmu + vkp*etabar_dp*eta_dp + xi_dp*xi_dp*dhard)
```

**After:**
```fortran
N_scalar = -(2.d0*vmu + vkp*etabar_dp*eta_dp + xi_dp*dhard/sqrt3)
```

**Rationale:** Matches the correct Drucker-Prager formulation from `stress_dp_rm.f` where the hardening term is ξ*(∂K/∂α)/√3, not ξ²*(∂K/∂α).

### 2. ✓ Incremental Updates in Box 1 (Lines 234-237, 265-266, 374-375)
**Before:**
```fortran
plstrg(:,:) = plstrg(:,:) + deltag*oun(:,:)*dsqrt(1.d0/2.d0)
alpeg_new = alpeg + xi_dp*deltag
```

**After:**
```fortran
plstrg(:,:) = plstrg_old(:,:) + delta_gamma_inc*oun(:,:)*dsqrt(1.d0/2.d0)
     &                + delta_gamma_inc*etabar_dp*DELTA(:,:)/3.d0
alpeg_new = alpeg_old + delta_gamma_inc/sqrt3
```

**Rationale:** Box 1 must use only the increment `delta_gamma_inc` from `histi(1)`, not the cumulative value, to properly integrate with the global Newton solver.

### 3. ✓ Added Stress Corrector `psig` (Lines 293-307)
**Added:**
```fortran
c       === Compute stress corrector psig (from st_gtn.f) ===
c       psig = -(g_val/N_scalar)*L_mat
        if(dabs(N_scalar).gt.1.d-10) then
          psig_factor = -g_val/N_scalar
          do jj=1,3
            do ii=1,3
              psig(ii,jj) = psig_factor * L_mat(ii,jj)
            enddo
          enddo
        else
          psig(:,:) = 0.d0
        endif
c       === Add stress corrector to final stress ===
        sig(:,:) = sig(:,:) + psig(:,:)
```

**Rationale:** The stress corrector is required for the Block Newton method to properly update the stress state, as demonstrated in `st_gtn.f`.

### 4. ✓ Added Geometric Softening Factors (Lines 283-291, 318, 325-328)
**Added:**
```fortran
if(delta_gamma_inc.gt.1.d-16 .and. stno.gt.1.d-16) then
  theta = 1.d0 - (dsqrt(2.d0)*vmu*delta_gamma_inc)/stno
  A_coeff = 1.d0/(vmu + vkp*etabar_dp*eta_dp + xi_dp*dhard_new/sqrt3)
  thetab = (dsqrt(2.d0)*vmu*delta_gamma_inc)/stno - vmu*A_coeff
else
  theta = 1.d0
  thetab = 0.d0
endif
```

**Tangent matrix updated:**
```fortran
bulk_coef = 1.d0 - vkp*eta_dp*etabar_dp*A_coeff
C_ijkl = bulk_coef*vkp*DELTA(ii,jj)*DELTA(kk,ll)
     &   + 2.d0*vmu*theta*( FIT(ii,jj,kk,ll)
     &   - (1.d0/3.d0)*DELTA(ii,jj)*DELTA(kk,ll) )
     &   + 2.d0*vmu*thetab*oun(ii,jj)*oun(kk,ll)
```

**Rationale:** These factors are essential for the consistent tangent to match the return mapping in the limit δγ→0, as shown in `stress_dp_rm.f`.

### 5. ✓ Synchronized `histi` Layout (Lines 345-369, 407-412)
**Storage layout:**
- `histi(1)`: δγ increment from phexa8 (INPUT - never overwritten)
- `histi(3)`: Updated αeq (equivalent plastic strain)
- `histi(4)`: Yield function value ftreg
- `histi(5)`: N_scalar
- `histi(6)`: g_val (current yield function)
- `histi(7-15)`: Current strain tensor
- `histi(16-24)`: xa sensitivity tensor (∂Δγ/∂ε)
- `histi(25-33)`: M_mat for element computation
- `histi(34)`: Cumulative Δγ

**Rationale:** Matches the layout expected by `phexa8.f` and follows the pattern from `stress_vm_bn.f`.

### 6. ✓ Added Initial State Preservation (Lines 71-72, 139)
**Added:**
```fortran
c     === Save initial state for Box 1 incremental update ===
      plstrg_old = plstrg
      alpeg_old = alpeg
```

**Usage in trial stress:**
```fortran
stry = 2.d0*vmu*(str -emean*DELTA -plstrg_old)
```

**Rationale:** Essential for incremental updates - must restore to initial state before applying new increment.

## Key Constants Added
- `sqrt3 = dsqrt(3.d0)` - Used throughout for clarity
- `plstrg_old`, `alpeg_old` - Store initial state
- `theta`, `thetab` - Geometric softening factors
- `bulk_coef` - Volumetric tangent coefficient
- `A_coeff` - Return mapping coefficient

## Files Modified
1. `src/stress_dp_bn.f` - Main implementation (replaced)
2. `src/stress_dp_bn_backup_before_fix.f` - Backup of original
3. `src/stress_dp_bn_fixed.f` - Working copy for reference

## Verification Checklist
- [x] Plastic updates use only `delta_gamma_inc` increment
- [x] N_scalar follows ξ/√3 hardening term
- [x] Stress corrector psig is computed and added
- [x] Consistent tangent includes geometric softening
- [x] histi slots properly synchronized
- [x] Elastic branch properly initialized

## Testing Status
- [ ] Compile with Intel Fortran
- [ ] Run 1elem_f.cml test case
- [ ] Verify single-iteration convergence in plastic regime
- [ ] Compare with return mapping results

## Next Steps
1. Compile the modified code
2. Run test case with MATYPE=5 (Block Newton)
3. Compare convergence behavior with return mapping
4. Verify stress-strain curves match expected results