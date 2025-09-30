#!/bin/bash
# Fix A_coeff initialization in stress_dp_bn.f

# Add A_coeff initialization after thetab = 0.d0 in the else clause (line 354)
sed -i '354a\          A_coeff = 1.d0/(vmu + vkp*etabar_dp*eta_dp     \&                   + xi_dp*dhard_new/sqrt3)' stress_dp_bn.f

# Also initialize A_coeff in the elastic case (around line 484)
sed -i '/g_val = 0.d0/a\c       Initialize tangent factors for elastic case        theta = 1.d0        thetab = 0.d0        A_coeff = 0.d0        bulk_coef = 1.d0' stress_dp_bn.f

echo Fixed A_coeff initialization
