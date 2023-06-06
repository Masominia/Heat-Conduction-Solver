"""Static figures for the report (matplotlib, PNG for the README)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Validated categorical order (fixed, never cycled) + recessive ink.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
NEUTRAL = "#8a8985"
CRITICAL = "#c4271a"

plt.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "figure.dpi": 150,
        "axes.axisbelow": True,
    }
)


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _extent(result):
    m = result["metrics"]
    return (m["x0"], m["x0"] + m["nx"] * m["dx"], m["y0"], m["y0"] + m["ny"] * m["dy"])


def _view(section, extent):
    """The part of the window around the heaters and probes."""
    boxes = [h["box"] for h in section["heaters"]] + [p["box"] for p in section["probes"]]
    lo, hi = min(b[0] for b in boxes), max(b[1] for b in boxes)
    pad = max(6.0, 0.25 * (hi - lo))
    return max(extent[0], lo - pad), min(extent[1], hi + pad), max(extent[2], -3.0), extent[3]


def _map(ax, result, section):
    extent = _extent(result)
    field, k = result["field"], result["conductivity"]
    image = ax.imshow(
        field, origin="lower", extent=extent, cmap="Oranges", aspect="equal", interpolation="nearest", vmin=0
    )
    ny, nx = k.shape
    xs = extent[0] + (np.arange(nx) + 0.5) * (extent[1] - extent[0]) / nx
    ys = extent[2] + (np.arange(ny) + 0.5) * (extent[3] - extent[2]) / ny
    logk = np.log10(k)
    if np.ptp(logk) > 0:
        levels = sorted({0.5 * (a + b) for a, b in zip(np.unique(logk)[:-1], np.unique(logk)[1:], strict=True)})
        ax.contour(xs, ys, logk, levels=levels, colors=[INK2], linewidths=0.5)
    x0, x1, y0, y1 = _view(section, extent)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.grid(False)
    ax.set(xlabel="x, um", ylabel="y, um")
    return image


def _profile(ax, result, section):
    """Temperature along the row of the first probe."""
    extent = _extent(result)
    m = result["metrics"]
    box = section["probes"][0]["box"] if section["probes"] else section["heaters"][0]["box"]
    row = min(int((0.5 * (box[2] + box[3]) - m["y0"]) / m["dy"]), int(m["ny"]) - 1)
    xs = extent[0] + (np.arange(int(m["nx"])) + 0.5) * m["dx"]
    ax.plot(xs, result["field"][row], color=SERIES[1], label=f"y = {m['y0'] + (row + 0.5) * m['dy']:.2f} um")
    names = list(result["temperature"])
    if 0 < len(names) <= 16:
        centres = [0.5 * (p["box"][0] + p["box"][1]) for p in section["probes"]]
        ax.plot(
            centres, [result["temperature"][n] for n in names], ls="none", marker="o", ms=4, color=INK, label="probes"
        )
    x0, x1, _, _ = _view(section, extent)
    ax.set_xlim(x0, x1)
    ax.set_ylim(bottom=0)
    ax.set(xlabel="x, um", ylabel="temperature rise, K")
    ax.legend(loc="upper right", fontsize=7.5)


def _step(ax, result):
    for i, (name, trace) in enumerate(result["traces"].items()):
        steady = trace[-1] if trace[-1] > 0 else 1.0
        ax.plot(result["time"], trace / steady, color=SERIES[i % len(SERIES)], label=name)
    ax.axhline(1 - np.exp(-1), color=NEUTRAL, lw=0.8, ls=":")
    ax.set_xlim(0, result["time"][-1])
    ax.set_ylim(0, 1.05)
    ax.set(xlabel="time, us", ylabel="step response over its last value")
    ax.legend(loc="lower right", fontsize=7.5)


def case_figure(result, section, title, path):
    fig, (left, right) = plt.subplots(1, 2, figsize=(8.8, 3.2), gridspec_kw={"width_ratios": [1.5, 1]})
    image = _map(left, result, section)
    fig.colorbar(image, ax=left, shrink=0.8, pad=0.02, label="K")
    if "traces" in result and len(result["traces"]) <= 4:
        _step(right, result)
    else:
        _profile(right, result, section)
    fig.suptitle(title if len(title) < 95 else title[:92] + "...", fontsize=9, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=105)
    plt.close(fig)


def overview(gallery, path):
    names = list(gallery)
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.0))
    for ax, name in zip(axes.ravel(), names, strict=False):
        result, section, title = gallery[name]
        image = _map(ax, result, section)
        fig.colorbar(image, ax=ax, shrink=0.7, pad=0.02, label="K")
        ax.set_title(title if len(title) < 56 else title[:53] + "...", fontsize=8.5)
    for ax in axes.ravel()[len(names) :]:
        ax.axis("off")
    _save(fig, path)
