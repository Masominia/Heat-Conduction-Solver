import numpy as np

from heatphase import studies


def pair(own, other):
    return {"response": {"wg1": {"h1": own}, "wg2": {"h1": other}}, "pi_power": {"wg1": {"h1": 20.0}}, "step": {}}


def test_crosstalk_fit_recovers_the_decay_length():
    results = {f"pitch_box2_p{p:02d}": pair(100.0, 100.0 * np.exp(-p / 3.0)) for p in (4, 6, 8, 10, 15, 20, 30)}
    results["shifter_standard"] = {}
    t = studies.crosstalk(results, "pitch_box2_p")
    assert list(t["pitch, um"]) == [4, 6, 8, 10, 15, 20, 30]
    assert abs(t.attrs["decay_um"] - 3.0) < 1e-9 and abs(t["crosstalk"][0] - np.exp(-4 / 3)) < 1e-12
    assert studies.crosstalk(results, "pitch_box9_p").empty


def test_tradeoff_and_trenches():
    def steady(p):
        return {"pi_power": {"wg": {"h": p}}}

    def transient(tau):
        return {"step": {"wg": {"tau": tau, "rise": 2 * tau}}}

    results = {
        "shifter_box1": steady(40.0),
        "transient_box1": transient(6.0),
        "shifter_standard": steady(24.0),
        "transient_standard": transient(10.0),
    }
    t = studies.oxide_tradeoff(results)
    assert list(t["product, nJ"]) == [240.0, 240.0]
    tr = studies.trenches(
        {"trench_none": pair(100.0, 8.0), "trench_box": {**pair(120.0, 1.2), "step": {"wg1": {"tau": 12.0}}}}
    )
    assert list(tr["variant"]) == ["none", "box"] or list(tr["variant"]) == ["box", "none"]
    assert tr.set_index("variant").loc["box", "crosstalk at 10 um"] == 0.01


def test_sensitivity_is_relative_to_the_baseline():
    def one(r):
        return {"response": {"wg": {"h": r}}}

    t = studies.sensitivity({"a": one(110.0), "b": one(100.0)}, {"a": "variant a", "missing": "not there"}, "b")
    assert len(t) == 1 and abs(t["change"][0] - 0.1) < 1e-12
    assert studies.sensitivity({}, {"a": "x"}, "b").empty
