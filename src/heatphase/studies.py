"""Results drawn from several cross-sections of the library at once."""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import reference


def _single(result: dict, probe: str) -> float:
    return next(iter(result["response"][probe].values()))


def crosstalk(results: dict[str, dict], prefix: str) -> pd.DataFrame:
    """Response of the heated waveguide and of its neighbour against pitch, with the fitted decay length."""
    rows = []
    for name, r in results.items():
        if name.startswith(prefix):
            own, other = _single(r, "wg1"), _single(r, "wg2")
            rows.append(
                {
                    "pitch, um": int(name.removeprefix(prefix)),
                    "own waveguide, K/(mW/um)": own,
                    "neighbour, K/(mW/um)": other,
                    "crosstalk": other / own,
                }
            )
    if len(rows) < 3:
        return pd.DataFrame()
    table = pd.DataFrame(rows).sort_values("pitch, um").reset_index(drop=True)
    mid = table[(table["pitch, um"] >= 6) & (table["pitch, um"] <= 20)]
    slope = np.polyfit(mid["pitch, um"], np.log(mid["crosstalk"]), 1)[0]
    table.attrs["decay_um"] = float(-1 / slope)
    return table


def oxide_tradeoff(results: dict[str, dict]) -> pd.DataFrame:
    """Power for pi and time constant against buried-oxide thickness."""
    pairs = {
        1: ("shifter_box1", "transient_box1"),
        2: ("shifter_standard", "transient_standard"),
        3: ("shifter_box3", "transient_box3"),
    }
    rows = []
    for box, (steady, transient) in pairs.items():
        if steady in results and transient in results:
            p_pi = next(iter(results[steady]["pi_power"]["wg"].values()))
            tau = results[transient]["step"]["wg"]["tau"]
            rows.append(
                {"buried oxide, um": box, "power for pi, mW": p_pi, "time constant, us": tau, "product, nJ": p_pi * tau}
            )
    return pd.DataFrame(rows)


def trenches(results: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for name, r in results.items():
        if name.startswith("trench_"):
            own, other = _single(r, "wg1"), _single(r, "wg2")
            rows.append(
                {
                    "variant": name.removeprefix("trench_"),
                    "power for pi, mW": next(iter(r["pi_power"]["wg1"].values())),
                    "crosstalk at 10 um": other / own,
                    "time constant, us": r["step"]["wg1"]["tau"] if r["step"] else float("nan"),
                }
            )
    return (
        pd.DataFrame(rows).sort_values("power for pi, mW", ascending=False).reset_index(drop=True)
        if rows
        else pd.DataFrame()
    )


def sensitivity(results: dict[str, dict], names: dict[str, str], baseline: str) -> pd.DataFrame:
    """Response of the waveguide for a set of variants, relative to a baseline."""
    if baseline not in results:
        return pd.DataFrame()
    base = _single(results[baseline], "wg")
    rows = [
        {
            "variant": label,
            "response, K/(mW/um)": _single(results[n], "wg"),
            "change": _single(results[n], "wg") / base - 1,
        }
        for n, label in names.items()
        if n in results
    ]
    return pd.DataFrame(rows)


def line_heaters(results: dict[str, dict], sections: dict[str, dict]) -> pd.DataFrame:
    """Solver against the closed form at every probe of the line-heater family."""
    rows = []
    for name, r in results.items():
        if not name.startswith("line_"):
            continue
        s = sections[name]
        k, d = s["materials"][s["background"]][0], -s["ymin"]
        for p in s["probes"]:
            x = p["box"][0]
            exact = 10 * float(reference.line_source_surface(0.1, k, d, x))
            rows.append(
                {
                    "section": name,
                    "distance over thickness": x / d,
                    "solver over closed form": r["response"][p["name"]]["h"] / exact,
                }
            )
    return pd.DataFrame(rows)


def surface(results: dict[str, dict], sections: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for name, r in results.items():
        if name.startswith("top_"):
            m = r["metrics"]
            rows.append(
                {
                    "film coefficient, W/(m2 K)": sections[name]["h"],
                    "share of heat leaving through the top": m["to_surface"] / m["injected"],
                    "power for pi, mW": next(iter(r["pi_power"]["wg"].values())),
                }
            )
    return (
        pd.DataFrame(rows).sort_values("film coefficient, W/(m2 K)").reset_index(drop=True) if rows else pd.DataFrame()
    )
