"""Solves the whole cross-section library, checks it and writes the report."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from . import cases, engine, plotting, reference, studies
from .compare import compare
from .engine import ROOT
from .report import md_table


def _check(area, check, ok, detail):
    return {"area": area, "check": check, "status": "PASS" if ok else "FAIL", "detail": detail}


def run(cfg_path=None, out=None, fast: bool = False) -> dict:
    cfg = yaml.safe_load(open(cfg_path or ROOT / "configs" / "default.yaml"))
    out = Path(out) if out else ROOT / "reports"
    (out / "figures" / "cases").mkdir(parents=True, exist_ok=True)
    if not engine.available():
        raise SystemExit("Rust solver not built: run `make build` first")
    paths = cases.discover(ROOT / cfg["cases_dir"])
    if fast:
        keep = tuple(cfg["fast_families"])
        paths = [
            p
            for i, p in enumerate(paths)
            if i % cfg["fast_every"] == 0 or p.stem in cfg["gallery"] or p.stem.startswith(keep)
        ]

    checks, rows, worst, results, sections = [], [], {}, {}, {}
    for path in paths:
        record, result = cases.run_case(path)
        family = cases.family_of(path.stem, cfg["families"])
        cmp = compare(path, result)
        record.update(
            temperature=result["temperature"],
            response=result["response"],
            pi_power=result["pi_power"],
            step=result["step"],
            against_reference=cmp["differences"],
        )
        for name, diff in cmp["differences"].items():
            worst[name] = max(worst.get(name, 0.0), diff)
        checks.append(
            _check(
                "independent solution",
                f"{path.stem}: responses, temperatures, heat flows, step response and the field solved again with SciPy",
                not cmp["outside"],
                "; ".join(cmp["outside"]) or f"field within {cmp['differences']['field']:.1e} of its peak",
            )
        )
        cases.write_record(record, out / "cases")
        section = reference.read_section(path)
        plotting.case_figure(result, section, record["title"], out / "figures" / "cases" / f"{path.stem}.png")
        for c in record["checks"]:
            checks.append(
                {"area": family, "check": f"{path.stem}: {c['check']}", "status": c["status"], "detail": c["detail"]}
            )
        results[path.stem], sections[path.stem] = result, section
        m = result["metrics"]
        first = next(iter(result["pi_power"].values()), {})
        rows.append(
            {
                "section": path.stem,
                "family": family,
                "cells": int(m["nx"] * m["ny"]),
                "heaters": len(section["heaters"]),
                "peak, K": m["peak"],
                "power for pi, mW": next(iter(first.values()))
                if len(section["heaters"]) == 1 and first
                else float("nan"),
                "checks": len(record["checks"]),
                "passed": sum(c["status"] == "PASS" for c in record["checks"]),
            }
        )

    table = pd.DataFrame(rows)
    fam = table.groupby("family", sort=False).agg(
        sections=("section", "size"), checks=("checks", "sum"), passed=("passed", "sum")
    )
    gallery = {k: (results[k], sections[k], results[k]["title"]) for k in cfg["gallery"] if k in results}
    if gallery:
        plotting.overview(gallery, out / "figures" / "overview.png")

    lines, study = [], {}
    for prefix, oxide in (("pitch_box2_p", 5.0), ("pitch_box3_p", 6.0)):
        t = studies.crosstalk(results, prefix)
        if len(t):
            values, slab = t["crosstalk"].to_numpy(), reference.crosstalk_decay_length(oxide)
            ok = bool((values[1:] < values[:-1]).all()) and 0.7 * slab < t.attrs["decay_um"] < 1.6 * slab
            checks.append(
                _check(
                    "crosstalk",
                    f"crosstalk falls with pitch at a rate set by the oxide thickness ({oxide:g} um of oxide)",
                    ok,
                    f"fitted decay length {t.attrs['decay_um']:.2f} um; a slab over a perfect sink gives {slab:.2f} um",
                )
            )
            lines += [
                f"### Crosstalk against pitch, {oxide - 3:g} um buried oxide",
                "",
                f"Fitted decay length between 6 and 20 um: {t.attrs['decay_um']:.2f} um (2 d / pi for {oxide:g} um of oxide over a perfect sink: {slab:.2f} um).",
                "",
                md_table(t, 4, index=False),
                "",
            ]
            study[prefix.rstrip("_p")] = {"decay_um": t.attrs["decay_um"], "rows": t.to_dict("records")}
    trade = studies.oxide_tradeoff(results)
    if len(trade) >= 3:
        spread = float(trade["product, nJ"].max() / trade["product, nJ"].min() - 1)
        checks.append(
            _check(
                "efficiency against speed",
                "power for pi times time constant changes little with the buried oxide",
                spread < 0.15,
                f"products {', '.join(f'{v:.0f}' for v in trade['product, nJ'])} nJ, spread {100 * spread:.0f}%",
            )
        )
        lines += ["### Efficiency against speed", "", md_table(trade, 3, index=False), ""]
        study["oxide_tradeoff"] = trade.to_dict("records")
    tr = studies.trenches(results)
    if len(tr):
        lines += ["### Trenches and undercut (neighbour at 10 um)", "", md_table(tr, 4, index=False), ""]
        study["trenches"] = json.loads(tr.to_json(orient="records"))
        named = tr.set_index("variant")
        if {"none", "box", "both"} <= set(named.index):
            ok = (
                named.loc["box", "crosstalk at 10 um"] < 0.5 * named.loc["none", "crosstalk at 10 um"]
                and named.loc["both", "power for pi, mW"] < named.loc["none", "power for pi, mW"]
            )
            checks.append(
                _check(
                    "trenches",
                    "a trench to the substrate more than halves the crosstalk; trenches on both sides also save power",
                    bool(ok),
                    f"crosstalk {named.loc['none', 'crosstalk at 10 um']:.4f} to {named.loc['box', 'crosstalk at 10 um']:.4f}",
                )
            )
    mesh = studies.sensitivity(results, {"mesh_d500": "500 nm cells", "mesh_d250": "250 nm cells"}, "mesh_d125")
    if len(mesh) == 2:
        change = float(mesh["change"].abs().iloc[1])
        checks.append(
            _check(
                "mesh",
                "halving the cells from 250 to 125 nm changes the response by less than 1%",
                change < 0.01,
                f"{100 * change:.2f}%",
            )
        )
        lines += ["### Mesh (relative to 125 nm cells)", "", md_table(mesh, 4, index=False), ""]
        study["mesh"] = mesh.to_dict("records")
    window = studies.sensitivity(
        results,
        {
            "window_w030": "30 um wide",
            "window_w040": "40 um wide",
            "window_w100": "100 um wide",
            "window_sub02": "2 um of substrate",
            "window_sub05": "5 um of substrate",
            "window_sub20": "20 um of substrate",
        },
        "shifter_standard",
    )
    if len(window):
        change = float(window["change"].abs().max())
        checks.append(
            _check(
                "window",
                "the response is insensitive to the size of the window and the depth of substrate",
                change < 0.02,
                f"largest change {100 * change:.2f}% over {len(window)} variants",
            )
        )
        lines += ["### Window (relative to 60 um wide, 10 um of substrate)", "", md_table(window, 4, index=False), ""]
        study["window"] = window.to_dict("records")
    if "solver_cg_standard" in results and "shifter_standard" in results:
        a, b = (next(iter(results[n]["response"]["wg"].values())) for n in ("solver_cg_standard", "shifter_standard"))
        checks.append(
            _check(
                "solvers",
                "conjugate gradients and banded Cholesky give the same response",
                abs(a / b - 1) < 1e-8,
                f"relative difference {abs(a / b - 1):.1e}, {int(results['solver_cg_standard']['metrics']['iterations'])} iterations",
            )
        )
    lh = studies.line_heaters(results, sections)
    if len(lh):
        near = lh[lh["distance over thickness"] <= 2.1]["solver over closed form"]
        checks.append(
            _check(
                "closed form",
                "line heater on a slab: solver within 3% of P/(pi k) ln coth(pi x / 4d) out to two thicknesses",
                bool(((near - 1).abs() < 0.03).all()),
                f"ratios from {near.min():.4f} to {near.max():.4f} at {len(near)} probes",
            )
        )
        lines += ["### Line heater on a slab: solver over closed form", "", md_table(lh, 4, index=False), ""]
        study["line_heaters"] = lh.to_dict("records")
    top = studies.surface(results, sections)
    if len(top):
        lines += ["### Heat leaving through the top surface", "", md_table(top, 5, index=False), ""]
        study["surface"] = top.to_dict("records")

    val = pd.DataFrame(checks)
    counts = val["status"].value_counts().to_dict()
    agreement = pd.DataFrame({"largest relative difference over the library": worst}).rename_axis("quantity")
    report = [
        "# HeatPhase report",
        "",
        f"{len(table)} cross-sections, {int(table['cells'].sum())} cells, {len(val)} checks: {counts.get('PASS', 0)} PASS, {counts.get('FAIL', 0)} FAIL.",
        "",
        "Every cross-section is a text file in `cases/`. Its expectations come from an independent solution of the same mesh or from a closed form, stated in the comments at the top of the file.",
        "Per-section results are in `reports/cases/`, temperature maps in `reports/figures/cases/`.",
        "",
        "## Families",
        "",
        md_table(fam, 0),
        "",
        "## Agreement with the independent solution",
        "",
        "Each cross-section is read, meshed, assembled and solved again by `reference.py` with SciPy's sparse direct solver. The field is compared with one solve with all heaters on, relative to its peak; the printed field has seven digits.",
        "",
        md_table(agreement, 12),
        "",
        "## Studies",
        "",
        *lines,
        "## All cross-sections",
        "",
        md_table(table.drop(columns=["family"]), 3, index=False),
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    (out / "validation.md").write_text("# Validation checks\n\n" + md_table(val, index=False) + "\n")
    metrics = {
        "sections": len(table),
        "cells": int(table["cells"].sum()),
        "families": fam.to_dict("index"),
        "largest_differences": worst,
        "studies": study,
        "sections_table": json.loads(table.set_index("section").to_json(orient="index")),
        "checks": checks,
    }
    (out / "metrics.json").write_text(json.dumps(metrics, indent=1, default=float))
    return metrics
