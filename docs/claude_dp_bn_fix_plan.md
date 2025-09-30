# Claude Code Tasking: Drucker-Prager Block Newton Rehab

## References to align with
- Yamamoto, Yamada & Matsui (2021), "Simultaneously iterative procedure based on block Newton method for elastoplastic problems" (Box 1 & Box 2).
- `src/stress_dp_rm.f` for Drucker-Prager return-mapping algebra (theta/thetab, xi/sqrt(3) hardening term, volumetric split).
- `src/stress_vm_bn.f` for working Box 1/Box 2 scaffolding and `histi` bookkeeping in the von Mises block Newton flow.
- `src/st_gtn.f` (lines ~29-391) for an existing multi-variable Block Newton implementation: note how it stores iteration history (`histi` slots 1-8, 9-35), computes stress correctors (`psig`) and pushes them back to the element residual.

## Current repository diagnosis
- `src/stress_dp_bn.f:13-36` allocates no buffer for the stress corrector `psig`. `st_gtn.f:372-396` shows that the block Newton assembly expects a pseudostress correction to be formed and added back into the updated Cauchy stress.
- `src/stress_dp_bn.f:142-170` builds `N_scalar` with `xi_dp*xi_dp*dhard_new`; the correct Drucker-Prager denominator, matching both the return mapping implementation and `st_gtn` use of mixed hardening, is `xi_dp*dhard_new/sqrt(3)`.
- `src/stress_dp_bn.f:111-140` updates the full plastic strain using `dmax1(deltagi,0.d0)` (cumulative value) instead of the incremental `delta_gamma_inc` handed in via `histi(1)`; Box 1 must apply only the latest increment, mirroring `st_gtn.f:61-87` where the prior state is restored before applying the new step.
- The current predictor path (`itr == 1`) does not seed `g_val`, `M_mat`, `xa`, or `psig`, so `phexa8.f` receives stale zeros and cannot form the Box 2 step 3.2 update. In contrast, `st_gtn.f:307-360` fully refreshes `histi(5:8)` with residuals and `histi(9:...)` with strain/gradient tensors each call.
- Tangent formation in `src/stress_dp_bn.f:180-221` drops the geometric softening factors (`theta`, `thetab`) present in `stress_dp_rm.f:327-355` and the pseudo-stress contribution seen in `st_gtn.f:372-391`, so the current tensor is inconsistent with both the return-mapping baseline and the block Newton requirements.

## Implementation directives for Claude Code
1. **Rebuild Box 1 using incremental updates**
   - Restore trial state, flow direction, and hardening exactly as in `stress_dp_rm.f`.
   - Apply only `delta_gamma_inc = histi(1)` to update `plstrg`, `alpeg`, and the volumetric pressure (`p_corr`). Store the cumulative value in `histi(34)` after the update, keeping `delta_gamma_inc` available for the global solve.
   - Reintroduce `theta` and `thetab` factors from `stress_dp_rm.f` when assembling the deviatoric part of `ctens` so the tangent matches the return mapping in the limit `delta_gamma_inc -> 0`.

2. **Fix the block Newton denominator**
   - Replace `xi_dp*xi_dp*dhard_new` with `xi_dp*dhard_new/sqrt(3)` when computing `N_scalar`.
   - Guard `N_scalar` against singular values using the pattern in `stress_vm_bn.f:329-343` (fallback to `-2*vmu`).

3. **Produce and return the stress corrector**
   - Allocate `psig(3,3)` in `stress_dp_bn.f`, compute it as `-(g_val/N_scalar)*L_mat` after the updated state is known (see `st_gtn.f:372-391` for structure), and add it to the corrected Cauchy stress.
   - Either add `psig` as an explicit output argument (preferred for clarity) or store it in reserved `histi` slots; update the element routines accordingly so the residual vector incorporates the correction.

4. **Synchronise `histi` layout with the element expectations**
   - On every call, write the predictor data so that `phexa8.f` can recover it: `[1]=delta_gamma_inc`, `[5]=N_scalar`, `[6]=g_val`, `[7-15]=current strain`, `[16-24]=xa`, `[25-33]=M_mat` (mirroring the existing von Mises implementation).
   - When `itr == 1`, fill these slots using the trial state (do not skip), just as `st_gtn` seeds its history before any Newton correction.

5. **Mirror diagnostic guardrails**
   - Follow `st_gtn` pattern of resetting `psig`, `xa`, and strain slots to zero in elastic steps so stale data never leaks into the next Newton iteration.
   - Keep comments concise, but document any non-obvious scaling (for example why `sqrt(0.5)` appears in the flow direction) since both `st_gtn` and `stress_dp_rm` rely on those conventions.

## Acceptance checklist
- Plastic updates consume only the new `delta_gamma_inc` and leave cumulative values stored for later use.
- `N_scalar` follows the `xi/sqrt(3)` hardening term and stays numerically stable in both elastic and plastic branches.
- Stress corrector `psig` is formed and either returned or stored so the element residual can add it, matching the expectation demonstrated in `st_gtn.f`.
- Consistent tangent reduces to the return-mapping tensor when `delta_gamma_inc -> 0`, and exhibits the expected asymmetry for non-associated flow once the block Newton terms are included.
- `histi` slots align with `phexa8.f` (and the precedent in `st_gtn.f`) on every iteration.

## Verification recommendations
- Rebuild and run the one-element Drucker-Prager block Newton test case (`run_analysis` with `MATYPE=5`). Confirm that the global Newton converges in one iteration once the plastic regime is reached, and that `delta_gamma_inc` reported by `phexa8.f` matches the internal update.
- Instrument temporary prints guarded by `nel_current` and `ig_current` to verify `N_scalar`, `g_val`, and `psig` during the first few plastic steps, similar to the diagnostics in `st_gtn.f`. Remove the prints after validation.
