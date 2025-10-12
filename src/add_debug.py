# Add debug output to stress_dp_bn.f
with open('stress_dp_bn.f', 'r') as f:
    lines = f.readlines()

# Find line 428 where histi(6) = g_val
for i, line in enumerate(lines):
    if 'histi(6) = g_val' in line:
        # Add debug output before this line
        debug_lines = [
            '        if(itr.le.3.and.istep.eq.13.and.idepg.eq.1) then\n',
            '          write(*,*) "DBG: itr=",itr," g_val=",g_val,\n',
            '     &               " hard_new=",hard_new\n',
            '        endif\n'
        ]
        lines[i:i] = debug_lines
        break

# Also add debug after receiving delta_gamma_inc
for i, line in enumerate(lines):
    if 'delta_gamma_inc = histi(1)' in line and 'phexa8' in lines[i-1]:
        debug_lines = [
            '        if(itr.le.3.and.istep.eq.13.and.idepg.eq.1) then\n',
            '          write(*,*) "DBG: Received delta_gamma_inc=",\n',
            '     &               delta_gamma_inc," from phexa8"\n',
            '        endif\n'
        ]
        lines[i+1:i+1] = debug_lines
        break

with open('stress_dp_bn.f', 'w') as f:
    f.writelines(lines)
    
print('Debug output added')
