"""Solve one cross-section file and print what the solver found.

python scripts/solve.py cases/trench_both.heat
"""

import sys

from heatphase import engine

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    r = engine.run(sys.argv[1], "check")
    m = r["metrics"]
    print(r["title"])
    print(
        f"  {int(m['nx'])} x {int(m['ny'])} cells, peak {m['peak']:.2f} K at ({m['peak_x']:.2f}, {m['peak_y']:.2f}) um"
    )
    print(f"  heat in {m['injected']:.3f} W/m, to the sink {m['to_sink']:.3f}, through the top {m['to_surface']:.3f}")
    for probe, row in r["response"].items():
        for heater, value in row.items():
            print(
                f"  {probe:<8} from {heater:<6} {value:10.4f} K per mW/um   power for pi {r['pi_power'][probe][heater]:10.2f} mW"
            )
    for probe, step in r["step"].items():
        print(f"  {probe:<8} step response: 63% after {step['tau']:.2f} us, 10 to 90% in {step['rise']:.2f} us")
    for c in r["checks"]:
        print(f"  {c['status']}  {c['check']}  ({c['detail']})")
