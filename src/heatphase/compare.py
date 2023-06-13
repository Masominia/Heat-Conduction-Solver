"""Compares what the solver reports with the independent SciPy solution."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import reference

TOLERANCES = {
    "response": 1e-8,
    "temperature": 1e-8,
    "peak": 1e-8,
    "balance": 1e-9,
    "outflow": 1e-8,
    "time": 1e-7,
    "field": 2e-6,
}


def _relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


def compare(path: Path | str, result: dict) -> dict:
    """Solves one cross-section again. Returns the quantities out of tolerance and the largest differences.

    The reference field is one solve with all heaters on, so it also tests that the solver's sum of
    single-heater fields is the field of the sum.
    """
    ref = reference.evaluate(reference.read_section(path))
    diff = {"response": 0.0, "temperature": 0.0, "time": 0.0}
    for probe, row in ref["response"].items():
        for heater, value in row.items():
            diff["response"] = max(diff["response"], _relative(result["response"][probe][heater], value))
    for probe, value in ref["temperature"].items():
        if abs(value) > 0:
            diff["temperature"] = max(diff["temperature"], _relative(result["temperature"][probe], value))
    m = result["metrics"]
    diff["peak"] = _relative(m["peak"], ref["peak"])
    diff["balance"] = abs((m["to_sink"] + m["to_surface"]) / m["injected"] - 1)
    diff["outflow"] = (
        max(abs(m["to_sink"] - ref["to_sink"]), abs(m["to_surface"] - ref["to_surface"])) / ref["injected"]
    )
    for probe, value in ref.get("tau", {}).items():
        for kind, want in (("tau", value), ("rise", ref["rise"][probe])):
            got = result["step"][probe][kind]
            if want == want or got == got:  # NaN on both sides means the run was too short for both
                diff["time"] = max(diff["time"], _relative(got, want) if want == want and got == got else 1.0)
    diff["field"] = float(np.max(np.abs(result["field"] - ref["field_all_on"])) / ref["peak"])
    outside = [f"{name} differs by {value:.2e}" for name, value in diff.items() if value > TOLERANCES[name]]
    return {"outside": outside, "differences": {k: float(v) for k, v in diff.items()}}
