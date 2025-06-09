#!/usr/bin/env python3
import re
import sys
from pathlib import Path

f90_path = Path('stress_dp.f')
doc_path = Path('docs/stress_variables.md')

f90_text = f90_path.read_text()
match = re.search(r"SUBROUTINE\s+stress_dp\(([^)]*)\)", f90_text, re.IGNORECASE)
if not match:
    print('Failed to find stress_dp signature in source')
    sys.exit(1)
fortran_args = [a.strip().lower() for a in match.group(1).replace('&', ' ').replace('\n', ' ').split(',') if a.strip()]
source_sig = f"subroutine stress_dp({', '.join(fortran_args)})"

md_text = doc_path.read_text()
match = re.search(r"```fortran\s*subroutine\s+stress_dp\(([^)]*)\)\s*```", md_text, re.IGNORECASE)
if not match:
    print('Failed to find stress_dp signature in documentation')
    sys.exit(1)
doc_args = [a.strip().lower() for a in match.group(1).replace('&', ' ').replace('\n', ' ').split(',') if a.strip()]
doc_sig = f"subroutine stress_dp({', '.join(doc_args)})"

if source_sig.lower() != doc_sig.lower():
    print('Documentation mismatch detected:')
    print('  Fortran :', source_sig)
    print('  Docs    :', doc_sig)
    sys.exit(1)
