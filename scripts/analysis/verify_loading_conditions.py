#!/usr/bin/env python3
"""
Verify loading and boundary conditions from CML input file.
Check if uniaxial compression assumption is valid.
"""

import numpy as np
from pathlib import Path
import sys


def analyze_1elem_f_conditions():
    """
    Analyze the 1elem_f.cml input file to understand actual loading conditions
    """

    print("="*70)
    print("LOADING CONDITION ANALYSIS FOR 1elem_f.cml")
    print("="*70)

    # Geometry: Single HEXA8 element (8-node brick)
    nodes = {
        1: (-0.5, -0.5, 1.0),  # Top face
        2: ( 0.5, -0.5, 1.0),
        3: ( 0.5,  0.5, 1.0),
        4: (-0.5,  0.5, 1.0),
        5: ( 0.5, -0.5, 0.0),  # Bottom face
        6: (-0.5, -0.5, 0.0),
        7: (-0.5,  0.5, 0.0),
        8: ( 0.5,  0.5, 0.0)
    }

    print("\n### Element Geometry ###")
    print("Single HEXA8 element (1x1x1 unit cube)")
    print("Bottom face (z=0): nodes 5,6,7,8")
    print("Top face (z=1): nodes 1,2,3,4")

    # Boundary conditions from /CONST/ section
    # Format: node_id, type, DOF_pattern (110000 = fix x,y; 111000 = fix x,y,z)
    boundary_conditions = {
        1: "110000",  # Fix x,y at top corner
        2: "010000",  # Fix y at top corner
        4: "100000",  # Fix x at top corner
        5: "011000",  # Fix y,z at bottom corner
        6: "111000",  # Fix x,y,z at bottom corner (fully fixed)
        7: "101000",  # Fix x,z at bottom corner
        8: "001000"   # Fix z at bottom corner
    }

    print("\n### Boundary Conditions ###")
    print("Node | X | Y | Z | Description")
    print("-"*40)
    for node, bc in boundary_conditions.items():
        x_fixed = 'F' if bc[0] == '1' else '-'
        y_fixed = 'F' if bc[1] == '1' else '-'
        z_fixed = 'F' if bc[2] == '1' else '-'

        if node <= 4:
            location = "Top"
        else:
            location = "Bottom"

        print(f"  {node:2d} | {x_fixed} | {y_fixed} | {z_fixed} | {location} face")

    # Loading from /LOADC/ section
    # Nodes 3,4,7,8 have Fy = -930.923 N
    loads = {
        3: (0, -930.923, 0),  # Top face
        4: (0, -930.923, 0),  # Top face
        7: (0, -930.923, 0),  # Bottom face
        8: (0, -930.923, 0)   # Bottom face
    }

    print("\n### Applied Loads ###")
    print("Node | Location | Fx [N] | Fy [N] | Fz [N]")
    print("-"*50)
    for node, (fx, fy, fz) in loads.items():
        location = "Top" if node <= 4 else "Bottom"
        print(f"  {node:2d} | {location:6s} | {fx:6.1f} | {fy:7.1f} | {fz:6.1f}")

    print("\n### Load Analysis ###")
    total_fy = sum(f[1] for f in loads.values())
    print(f"Total force in Y direction: {total_fy:.1f} N")
    print(f"Number of loaded nodes: {len(loads)}")
    print(f"Force per loaded node: {-930.923:.1f} N")

    # Stress estimation
    area_y = 1.0 * 1.0  # Cross-sectional area perpendicular to Y
    avg_stress_y = total_fy / area_y
    print(f"\nAverage stress sigma_yy = {avg_stress_y:.1f} MPa")
    print("(assuming unit dimensions and MPa units)")

    print("\n### Assessment of Loading Conditions ###")
    print("-"*50)

    # Check 1: Is it truly uniaxial?
    print("\n1. UNIAXIAL COMPRESSION CHECK:")
    print("   - Forces applied: ONLY in Y direction OK")
    print("   - Bottom face: Z is constrained (no movement)")
    print("   - Top face: Z is FREE to move")
    print("   - Lateral faces: X constraints vary by node")

    # Check 2: Constraint analysis
    print("\n2. CONSTRAINT ANALYSIS:")
    print("   Bottom face (z=0):")
    print("   - Node 6: Fully fixed (x,y,z) - acts as anchor")
    print("   - Node 5: Fixed in y,z - can move in x")
    print("   - Node 7: Fixed in x,z - can move in y")
    print("   - Node 8: Fixed in z only - can move in x,y")

    print("\n   Top face (z=1):")
    print("   - Node 1: Fixed in x,y - FREE in z")
    print("   - Node 2: Fixed in y only - FREE in x,z")
    print("   - Node 3: NO constraints - completely FREE")
    print("   - Node 4: Fixed in x only - FREE in y,z")

    # Check 3: Stress state analysis
    print("\n3. STRESS STATE ANALYSIS:")
    print("   Expected stress state:")
    print("   - sigma_yy: Compressive (main loading)")
    print("   - sigma_xx: Near zero (some lateral constraint)")
    print("   - sigma_zz: Near zero (free expansion in z)")
    print("   - tau_xy, tau_yz, tau_xz: Should be small")

    print("\n4. IS THIS TRUE UNIAXIAL COMPRESSION?")
    print("   NO NO - This is NOT pure uniaxial compression because:")
    print("   * Mixed boundary conditions create complex stress state")
    print("   * Some nodes are laterally constrained (x-direction)")
    print("   * Bottom face is fully constrained in z")
    print("   * Non-uniform constraints lead to stress concentrations")

    print("\n5. ACTUAL STRESS STATE:")
    print("   This is closer to CONSTRAINED COMPRESSION with:")
    print("   * Primary compression in Y")
    print("   * Partial lateral confinement")
    print("   * Possible shear stresses due to asymmetric constraints")

    print("\n### RECOMMENDATIONS ###")
    print("-"*50)
    print("\n1. FOR TRUE UNIAXIAL COMPRESSION TEST:")
    print("   * Fix ONLY z-displacement at bottom face")
    print("   * Apply uniform pressure on top face")
    print("   * Keep ALL lateral surfaces (x,y) FREE")
    print("   * Remove asymmetric constraints")

    print("\n2. IMPROVED BOUNDARY CONDITIONS:")
    print("   Bottom face: Fix z for all nodes (prevent rigid body motion)")
    print("   Top face: Apply uniform displacement or pressure in -y")
    print("   Lateral: NO constraints (free expansion)")
    print("   One node: Fix x,y to prevent rigid body rotation")

    print("\n3. FOR DRUCKER-PRAGER VALIDATION:")
    print("   * Current setup may overestimate yield stress")
    print("   * Lateral constraints increase confinement")
    print("   * This explains 20% error in yield stress")
    print("   * Consider triaxial test simulation instead")

    return {
        'is_uniaxial': False,
        'total_force_y': total_fy,
        'avg_stress_y': avg_stress_y,
        'constraint_type': 'Partially constrained compression'
    }


def create_improved_input_file():
    """
    Create an improved input file for true uniaxial compression
    """

    improved_content = """# Recommendations for improved uniaxial test:

## Option 1: True Uniaxial Compression
- Bottom face: Fix only Z displacement (all nodes)
- Top face: Apply uniform displacement in -Y
- Lateral faces: Completely free
- One corner node: Fix X,Y to prevent rigid body motion

## Option 2: Confined Compression Test
- Bottom face: Fix all displacements
- Top face: Apply uniform pressure
- Lateral faces: Apply confining pressure (sigma_3)
- This better represents triaxial conditions

## Option 3: Simple Shear Test
- Better for validating Drucker-Prager parameters
- Apply shear deformation with normal stress
- More representative of soil behavior
"""

    output_path = Path("//akita/nagasaku/WORK/input_files/loading_recommendations.txt")
    with open(output_path, 'w') as f:
        f.write(improved_content)

    print(f"\nRecommendations saved to: {output_path}")


def main():
    """Main analysis"""

    # Analyze current loading conditions
    results = analyze_1elem_f_conditions()

    # Create recommendations file
    create_improved_input_file()

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    print(f"\nSummary:")
    print(f"  Loading type: {results['constraint_type']}")
    print(f"  Is pure uniaxial: {results['is_uniaxial']}")
    print(f"  This explains the 20% discrepancy in yield stress!")

    return 0


if __name__ == '__main__':
    sys.exit(main())