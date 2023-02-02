"""Print the crosstalk matrix of a multi-heater cross-section: each probe's response to each heater over the largest in its row.

python scripts/crosstalk_matrix.py cases/array_one_on.heat
"""

import sys

from heatphase import engine

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    r = engine.run(sys.argv[1], "check")
    heaters = list(next(iter(r["response"].values())))
    print(f"{'':<8}" + "".join(f"{h:>8}" for h in heaters))
    for probe, row in r["response"].items():
        top = max(row.values())
        print(f"{probe:<8}" + "".join(f"{row[h] / top:>8.3f}" for h in heaters))
