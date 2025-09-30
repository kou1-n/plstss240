#!/usr/bin/env python3
# Fix A_coeff initialization in stress_dp_bn.f

with open('stress_dp_bn_tmp.f', 'r') as f:
    lines = f.readlines()

# Find and fix the plastic case (around line 354)
for i in range(len(lines)):
    if 'thetab = 0.d0' in lines[i] and i < 400:
        # Add A_coeff initialization after this line
        if i+1 < len(lines) and 'A_coeff' not in lines[i+1]:
            lines.insert(i+1, '          A_coeff = 1.d0/(vmu + vkp*etabar_dp*eta_dp\n')
            lines.insert(i+2, '     &                   + xi_dp*dhard_new/sqrt3)\n')
        break

# Find and fix the elastic case (search for g_val = 0.d0 in elastic section)
for i in range(450, min(500, len(lines))):
    if 'g_val = 0.d0' in lines[i]:
        # Add initialization after this line
        if i+1 < len(lines) and 'theta' not in lines[i+1]:
            lines.insert(i+1, 'c       Initialize tangent factors for elastic case\n')
            lines.insert(i+2, '        theta = 1.d0\n')
            lines.insert(i+3, '        thetab = 0.d0\n')
            lines.insert(i+4, '        A_coeff = 0.d0\n')
            lines.insert(i+5, '        bulk_coef = 1.d0\n')
        break

with open('stress_dp_bn.f', 'w') as f:
    f.writelines(lines)

print('Fixed A_coeff initialization')
