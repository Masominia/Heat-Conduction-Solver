"""Draw the temperature field of one cross-section.

python scripts/plot_section.py cases/trench_undercut.heat field.png
"""

import sys

from heatphase import engine, plotting, reference

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    result = engine.run(sys.argv[1])
    plotting.case_figure(result, reference.read_section(sys.argv[1]), result["title"], sys.argv[2])
    print(f"wrote {sys.argv[2]}")
