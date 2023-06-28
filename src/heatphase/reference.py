"""An independent solution of every cross-section with SciPy's sparse solvers, and closed forms.

It shares no code with the Rust solver: its own file reader, mesh painting, a vectorised assembly of the
conduction matrix as a sparse matrix, SuperLU instead of banded Cholesky, and its own time stepping.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

MATERIALS = {
    "si": (148.0, 1.63e6),
    "silicon": (148.0, 1.63e6),
    "sio2": (1.38, 1.61e6),
    "silica": (1.38, 1.61e6),
    "oxide": (1.38, 1.61e6),
    "sin": (20.0, 2.2e6),
    "nitride": (20.0, 2.2e6),
    "tin": (28.0, 3.2e6),
    "al": (237.0, 2.42e6),
    "aluminium": (237.0, 2.42e6),
    "air": (0.026, 1.2e3),
}


def read_section(path: Path | str) -> dict:
    """Reads a cross-section file into a dictionary, ignoring its expectations."""
    s = {
        "title": "",
        "dx": 0.25,
        "dy": 0.25,
        "width": 60.0,
        "ymin": -20.0,
        "ymax": 5.0,
        "background": "sio2",
        "materials": dict(MATERIALS),
        "rects": [],
        "heaters": [],
        "probes": [],
        "h": 0.0,
        "wavelength": 1.55,
        "dndt": 1.86e-4,
        "transient": None,
    }
    for raw in Path(path).read_text().splitlines():
        w = raw.split("#")[0].split()
        if not w or w[0] in ("expect", "solver"):
            continue
        kv = dict(x.split("=") for x in w[1:] if "=" in x)
        if w[0] == "title":
            s["title"] = " ".join(w[1:])
        elif w[0] in ("grid", "window", "optics"):
            s.update({k: float(v) for k, v in kv.items()})
        elif w[0] == "material":
            s["materials"][w[1]] = (float(kv["k"]), float(kv["c"]))
        elif w[0] == "background":
            s["background"] = w[1]
        elif w[0] == "layer":
            s["rects"].append((w[1], -1e9, 1e9, float(w[2]), float(w[3])))
        elif w[0] == "rect":
            s["rects"].append((w[1], *(float(v) for v in w[2:6])))
        elif w[0] == "heater":
            box = tuple(float(kv[k]) for k in ("xmin", "xmax", "ymin", "ymax"))
            s["heaters"].append({"name": w[1], "box": box, "power": float(kv.get("power", 0.1))})
        elif w[0] == "probe":
            box = (
                (float(kv["x"]),) * 2 + (float(kv["y"]),) * 2
                if "x" in kv
                else tuple(float(kv[k]) for k in ("xmin", "xmax", "ymin", "ymax"))
            )
            s["probes"].append({"name": w[1], "box": box})
        elif w[0] == "top":
            s["h"] = float(kv["h"]) if w[1] == "convective" else 0.0
        elif w[0] == "transient":
            s["transient"] = {"heater": kv["heater"], "dt": float(kv["dt"]), "steps": int(float(kv["steps"]))}
    return s


def mesh(s: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cell-centre coordinates and the conductivity and capacity of every cell (rows bottom to top)."""
    nx, ny = round(s["width"] / s["dx"]), round((s["ymax"] - s["ymin"]) / s["dy"])
    x = -0.5 * s["width"] + (np.arange(nx) + 0.5) * s["dx"]
    y = s["ymin"] + (np.arange(ny) + 0.5) * s["dy"]
    k = np.full((ny, nx), s["materials"][s["background"]][0])
    c = np.full((ny, nx), s["materials"][s["background"]][1])
    for name, xa, xb, ya, yb in s["rects"]:
        inside = np.outer((y >= ya) & (y <= yb), (x >= xa) & (x <= xb))
        k[inside], c[inside] = s["materials"][name]
    return x, y, k, c


def coverage(s: dict, x: np.ndarray, y: np.ndarray, box: tuple) -> np.ndarray:
    """Area of each cell covered by a rectangle, um^2."""
    xa, xb, ya, yb = box
    fx = np.clip(np.minimum(x + 0.5 * s["dx"], xb) - np.maximum(x - 0.5 * s["dx"], xa), 0, None)
    fy = np.clip(np.minimum(y + 0.5 * s["dy"], yb) - np.maximum(y - 0.5 * s["dy"], ya), 0, None)
    return np.outer(fy, fx)


def edge_conductances(s: dict, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Conductance per unit length from the bottom row to the sink and from the top row to the surroundings."""
    bottom = k[0] * s["dx"] / (0.5 * s["dy"])
    top = (
        1.0 / (1.0 / (s["h"] * s["dx"] * 1e-6) + 0.5 * s["dy"] / (k[-1] * s["dx"]))
        if s["h"] > 0
        else np.zeros(k.shape[1])
    )
    return bottom, top


def conduction_matrix(s: dict, k: np.ndarray) -> sp.csc_matrix:
    ny, nx = k.shape
    number = np.arange(nx * ny).reshape(ny, nx)
    gx = 2 * k[:, :-1] * k[:, 1:] / (k[:, :-1] + k[:, 1:]) * s["dy"] / s["dx"]
    gy = 2 * k[:-1, :] * k[1:, :] / (k[:-1, :] + k[1:, :]) * s["dx"] / s["dy"]
    diagonal = np.zeros((ny, nx))
    diagonal[:, :-1] += gx
    diagonal[:, 1:] += gx
    diagonal[:-1, :] += gy
    diagonal[1:, :] += gy
    bottom, top = edge_conductances(s, k)
    diagonal[0] += bottom
    diagonal[-1] += top
    rows = np.concatenate(
        [number[:, :-1].ravel(), number[:, 1:].ravel(), number[:-1, :].ravel(), number[1:, :].ravel(), number.ravel()]
    )
    cols = np.concatenate(
        [number[:, 1:].ravel(), number[:, :-1].ravel(), number[1:, :].ravel(), number[:-1, :].ravel(), number.ravel()]
    )
    vals = np.concatenate([-gx.ravel(), -gx.ravel(), -gy.ravel(), -gy.ravel(), diagonal.ravel()])
    return sp.csc_matrix((vals, (rows, cols)), shape=(nx * ny, nx * ny))


def _probe(s, x, y, box, field):
    area = coverage(s, x, y, box)
    if area.sum() > 0:
        return float((field * area).sum() / area.sum())
    ix, iy = int((box[0] + 0.5 * s["width"]) / s["dx"]), int((box[2] - s["ymin"]) / s["dy"])
    return float(field[iy, ix])


def evaluate(s: dict) -> dict:
    """Every quantity the solver reports for a cross-section, computed here from scratch."""
    x, y, k, c = mesh(s)
    a = conduction_matrix(s, k)
    lu = splu(a)
    units = []
    for h in s["heaters"]:
        area = coverage(s, x, y, h["box"])
        box_area = (h["box"][1] - h["box"][0]) * (h["box"][3] - h["box"][2])
        units.append((lu.solve((1e3 * area / box_area).ravel()).reshape(k.shape), 1e3 * area.sum() / box_area))
    field = sum(h["power"] * u for h, (u, _) in zip(s["heaters"], units, strict=True))
    injected = sum(h["power"] * q for h, (_, q) in zip(s["heaters"], units, strict=True))
    response = {
        p["name"]: {h["name"]: _probe(s, x, y, p["box"], u) for h, (u, _) in zip(s["heaters"], units, strict=True)}
        for p in s["probes"]
    }
    bottom, top = edge_conductances(s, k)
    out = {
        "nx": k.shape[1],
        "ny": k.shape[0],
        "response": response,
        "pi_power": {
            p: {h: s["wavelength"] / (2 * s["dndt"] * r) for h, r in row.items()} for p, row in response.items()
        },
        "temperature": {p["name"]: _probe(s, x, y, p["box"], field) for p in s["probes"]},
        "peak": float(field.max()),
        "injected": float(injected),
        "to_sink": float((bottom * field[0]).sum()),
        "to_surface": float((top * field[-1]).sum()),
        "field": field,
    }
    # one solve with every heater on: the field of the sum, not the sum of the fields
    total = sum(
        h["power"] * 1e3 * coverage(s, x, y, h["box"]) / ((h["box"][1] - h["box"][0]) * (h["box"][3] - h["box"][2]))
        for h in s["heaters"]
    )
    out["field_all_on"] = lu.solve(total.ravel()).reshape(k.shape)
    if s["transient"]:
        t = s["transient"]
        index = [h["name"] for h in s["heaters"]].index(t["heater"])
        heater = s["heaters"][index]
        source = heater["power"] * (a @ units[index][0].ravel())
        capacity = (c * s["dx"] * s["dy"] * 1e-12).ravel() / (t["dt"] * 1e-6)
        stepper = splu((a + sp.diags(capacity)).tocsc())
        temperature = np.zeros(a.shape[0])
        time, traces = [], {p["name"]: [] for p in s["probes"]}
        for step in range(1, t["steps"] + 1):
            temperature = stepper.solve(capacity * temperature + source)
            time.append(step * t["dt"])
            for p in s["probes"]:
                traces[p["name"]].append(_probe(s, x, y, p["box"], temperature.reshape(k.shape)))
        out["tau"], out["rise"] = {}, {}
        for p in s["probes"]:
            steady = heater["power"] * response[p["name"]][heater["name"]]
            levels = [crossing(time, traces[p["name"]], f * steady) for f in (1 - np.exp(-1), 0.1, 0.9)]
            out["tau"][p["name"]], out["rise"][p["name"]] = levels[0], levels[2] - levels[1]
        out["time"], out["traces"] = np.array(time), {n: np.array(v) for n, v in traces.items()}
    return out


def crossing(time, trace, level: float) -> float:
    """First time a rising trace that starts at zero reaches a level, by linear interpolation; NaN if never."""
    t0, v0 = 0.0, 0.0
    for t, v in zip(time, trace, strict=True):
        if v >= level:
            return t0 + (level - v0) / (v - v0) * (t - t0) if v > v0 else t
        t0, v0 = t, v
    return float("nan")


# closed forms


def stack_temperature(layers: list[tuple[float, float]], power: float, width: float) -> float:
    """Top temperature of a layered stack heated uniformly over its whole width: flux times the sum of t / k."""
    return 1e3 * power / width * sum(t / k for t, k in layers)


def line_source_surface(power: float, conductivity: float, thickness: float, x):
    """Surface temperature at a distance x from a line heater on a slab over a sink: P'/(pi k) ln coth(pi|x|/(4d))."""
    u = np.pi * np.abs(x) / (4 * thickness)
    return 1e3 * power / (np.pi * conductivity) * np.log(1 / np.tanh(u))


def line_source_images(power: float, conductivity: float, thickness: float, x: float, terms: int = 20000) -> float:
    """The same temperature from image sources of alternating sign at heights 2 n d, as a check on the closed form.

    The series alternates and converges slowly; the mean of two consecutive partial sums is used.
    """
    n = np.arange(1, terms + 2)
    partial = -np.log(abs(x)) - np.cumsum((-1.0) ** n * np.log(x * x + (2 * n * thickness) ** 2))
    return float(1e3 * power / (np.pi * conductivity) * 0.5 * (partial[-1] + partial[-2]))


def crosstalk_decay_length(thickness: float) -> float:
    return 2 * thickness / np.pi


def slab_time_constant(thickness: float, conductivity: float, capacity: float) -> float:
    """Slowest time constant of a slab with an adiabatic top over a sink, microseconds."""
    return 4 * thickness**2 / (np.pi**2 * conductivity / capacity * 1e6)


def slab_step_response(thickness: float, conductivity: float, capacity: float, t) -> np.ndarray:
    """Top temperature of a slab heated uniformly at its top surface from t = 0, as a fraction of its steady value."""
    t = np.asarray(t, dtype=float)
    alpha = conductivity / capacity * 1e6
    out = np.ones_like(t)
    for m in range(400):
        q = (2 * m + 1) * np.pi / (2 * thickness)
        out -= 8 / ((2 * m + 1) ** 2 * np.pi**2) * np.exp(-alpha * q * q * t)
    return out
