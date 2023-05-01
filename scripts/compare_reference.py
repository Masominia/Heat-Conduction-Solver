"""Solve one cross-section again with SciPy and print the differences.

python scripts/compare_reference.py cases/array_ramp.heat
"""

import sys

from heatphase import engine
from heatphase.compare import compare

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    cmp = compare(sys.argv[1], engine.run(sys.argv[1]))
    for name, diff in cmp["differences"].items():
        print(f"{name:<12}{diff:>12.1e}")
    sys.exit(1 if cmp["outside"] else 0)
