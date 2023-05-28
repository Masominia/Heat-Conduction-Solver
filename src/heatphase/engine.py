"""Python side of the Rust solver: runs the binary on a cross-section file and parses its output."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BINARY = ROOT / "engine" / "target" / "release" / "heatphase"


def available() -> bool:
    return BINARY.exists()


def run(section: Path | str, command: str = "run") -> dict:
    """Solve a cross-section. Returns metrics, probes, responses, checks and, for `run`, the fields.

    Raises RuntimeError with the solver's message if the file is invalid.
    """
    proc = subprocess.run([str(BINARY), command, str(section)], capture_output=True, text=True)
    if proc.returncode == 2:
        raise RuntimeError(proc.stderr.strip())
    res: dict = {
        "metrics": {},
        "temperature": {},
        "response": {},
        "pi_power": {},
        "step": {},
        "checks": [],
        "title": "",
    }
    conductivity, field, trace = [], [], []
    for line in proc.stdout.splitlines():
        key, _, rest = line.partition(" ")
        if key == "field":
            field.append(rest)
        elif key == "conductivity":
            conductivity.append(rest)
        elif key == "trace":
            trace.append(rest.split())
        elif key == "probe":
            name, value = rest.split()
            res["temperature"][name] = float(value)
        elif key == "response":
            probe, heater, value, pi = rest.split()
            res["response"].setdefault(probe, {})[heater] = float(value)
            res["pi_power"].setdefault(probe, {})[heater] = float(pi)
        elif key == "step":
            probe, tau, rise = rest.split()
            res["step"][probe] = {"tau": float(tau), "rise": float(rise)}
        elif key == "check":
            status, _, tail = rest.partition(" ")
            name, _, detail = tail.partition(" :: ")
            res["checks"].append({"check": name, "status": status, "detail": detail})
        elif key == "title":
            res["title"] = rest
        else:
            res["metrics"][key] = float(rest)
    if field:
        res["field"] = np.array([r.split() for r in field], dtype=float)  # rows bottom to top
        res["conductivity"] = np.array([r.split() for r in conductivity], dtype=float)
    if trace:
        table = np.array(trace, dtype=float)
        res["time"], res["traces"] = table[:, 0], dict(zip(res["temperature"], table[:, 1:].T, strict=True))
    return res
