# Plastic State Monitoring

## Overview

The plstss program provides step-level monitoring of plastic deformation onset through the NOR_*.txt output files. This feature allows users to track when plasticity first occurs during the loading process.

## Implementation

### Step-Level Detection

The plastic state detection is implemented in `output.f` using the following mechanism:

1. **Plastic Energy Tracking**: Monitors `tene_p` (total plastic energy) at each step
2. **Comparison Logic**: Compares current step plastic energy with previous step
3. **Static Variable**: Uses `tene_p_prev` to store previous step plastic energy
4. **Detection Condition**: 
   ```fortran
   if(tene_p.gt.1.d-12 .and. tene_p_prev.le.1.d-12) then
     WRITE(lwd,'(a,i3)') 'First plastic deformation at step: ',in
   endif
   ```

### Output Format

When plasticity first occurs, the following message is written to the NOR_*.txt file:

```
First plastic deformation at step:   4
```

This message appears immediately after the step information line in the convergence monitoring output.

## Example Output

```
%Step Iter     dfact          unorm          tnorm
    0    0  0.0000000E+00  0.0000000E+00  0.0000000E+00
    1    1  0.1000000E+00  0.8957688E-03  0.0000000E+00
    2    1  0.2000000E+00  0.1791538E-02  0.0000000E+00
    3    1  0.3000000E+00  0.2687306E-02  0.0000000E+00
    4    4  0.4000000E+00  0.3626089E-02  0.0000000E+00
First plastic deformation at step:   4
    5    4  0.5000000E+00  0.4710834E-02  0.0000000E+00
```

## Interpretation

- **Steps 1-3**: Elastic behavior (converges in 1 iteration)
- **Step 4**: First plastic deformation (requires 4 iterations)
- **Step 5+**: Continued plastic behavior

## Advantages of Step-Level Monitoring

### Before (Integration Point Level)
- Up to 64 redundant messages (8 elements × 8 Gauss points)
- Difficult to determine onset step
- Cluttered output files

### After (Step Level)
- Single, clear message per analysis
- Precise step identification
- Clean, readable output

## Technical Details

### Modified Files

1. **`stress_vm.f`**: Removed per-integration-point output
2. **`stress_dp.f`**: Removed per-integration-point output  
3. **`output.f`**: Added step-level detection logic

### Detection Threshold

The threshold `1.d-12` is used to account for numerical precision in floating-point calculations while detecting meaningful plastic deformation.

### Energy-Based Detection

Using plastic energy (`tene_p`) for detection provides:
- Global perspective across all elements
- Reliable indication of actual plastic work
- Immunity to local numerical artifacts

## Usage

No additional user input is required. The plastic monitoring is automatically active and will output to the standard NOR_*.txt file generated during analysis.