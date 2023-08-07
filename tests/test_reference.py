import numpy as np
import pytest

from heatphase import engine, reference
from heatphase.compare import compare
from heatphase.engine import ROOT


def section(tmp_path, text):
    path = tmp_path / "s.heat"
    path.write_text(text)
    return reference.read_section(path)


def test_reader_understands_every_statement(tmp_path):
    s = section(
        tmp_path,
        "# note\ntitle A B\ngrid dx=0.5 dy=0.125\nwindow width=40 ymin=-5 ymax=6\nmaterial polymer k=0.2 c=1.8e6\nbackground polymer\nlayer si -5 0\n"
        "rect si -0.25 0.25 2 2.25 # core\nheater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.2\nprobe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\n"
        "probe point x=3 y=1\ntop convective h=50\nsolver cg\noptics wavelength=1.31 dndt=1.9e-4\ntransient heater=h1 dt=0.5 steps=100\nexpect peak 30 0.5\n",
    )
    assert (s["title"], s["dx"], s["dy"], s["width"], s["ymin"], s["ymax"]) == ("A B", 0.5, 0.125, 40.0, -5.0, 6.0)
    assert s["background"] == "polymer" and s["materials"]["polymer"] == (0.2, 1.8e6)
    assert s["rects"][0][0] == "si" and s["rects"][1] == ("si", -0.25, 0.25, 2.0, 2.25)
    assert s["heaters"] == [{"name": "h1", "box": (-1.0, 1.0, 4.0, 4.25), "power": 0.2}]
    assert s["probes"][1]["box"] == (3.0, 3.0, 1.0, 1.0)
    assert (s["h"], s["wavelength"], s["dndt"]) == (50.0, 1.31, 1.9e-4)
    assert s["transient"] == {"heater": "h1", "dt": 0.5, "steps": 100}


def test_matrix_is_symmetric_and_conserves_heat(tmp_path):
    s = section(
        tmp_path,
        "grid dx=0.5 dy=0.5\nwindow width=10 ymin=-2 ymax=3\nlayer si -2 0\nheater h xmin=-1 xmax=1 ymin=2 ymax=2.5\n",
    )
    _, _, k, _ = reference.mesh(s)
    a = reference.conduction_matrix(s, k)
    assert abs(a - a.T).max() < 1e-12
    # row sums vanish except in the bottom row, where they are the conductance to the sink
    sums = np.asarray(a.sum(axis=1)).reshape(k.shape)
    assert np.abs(sums[1:]).max() < 1e-9 and np.allclose(sums[0], reference.edge_conductances(s, k)[0])


def test_layered_stack_is_exact(tmp_path):
    s = section(
        tmp_path,
        "grid dx=0.5 dy=0.25\nwindow width=8 ymin=-1 ymax=2\nlayer si -1 0\nheater h xmin=-4 xmax=4 ymin=1.75 ymax=2\nprobe p x=0.25 y=0.875\n",
    )
    e = reference.evaluate(s)
    expected = reference.stack_temperature([(1.0, 148.0), (0.875, 1.38)], 1.0, 8.0)
    assert e["response"]["p"]["h"] == pytest.approx(expected, rel=1e-12)
    assert e["to_sink"] == pytest.approx(e["injected"], rel=1e-12)
    assert np.allclose(e["field"], e["field_all_on"])


def test_closed_forms_agree_with_each_other():
    for x in (0.3, 1.0, 3.0):
        assert reference.line_source_images(0.1, 1.38, 2.0, x) == pytest.approx(
            float(reference.line_source_surface(0.1, 1.38, 2.0, x)), abs=1e-3
        )
    far = reference.line_source_surface(0.1, 1.38, 2.0, np.array([10.0, 12.0]))
    assert far[1] / far[0] == pytest.approx(np.exp(-2.0 / reference.crosstalk_decay_length(2.0)), rel=1e-3)
    tau = reference.slab_time_constant(2.0, 1.38, 1.61e6)
    late = reference.slab_step_response(2.0, 1.38, 1.61e6, [3 * tau, 4 * tau])
    assert (1 - late[0]) / (1 - late[1]) == pytest.approx(np.e, rel=1e-3)
    assert reference.slab_step_response(2.0, 1.38, 1.61e6, [1e-9])[0] == pytest.approx(0.0, abs=2e-3)


def test_crossing():
    assert reference.crossing([1.0, 2.0], [0.2, 0.6], 0.4) == pytest.approx(1.5)
    assert np.isnan(reference.crossing([1.0, 2.0], [0.2, 0.6], 0.9))


@pytest.mark.skipif(not engine.available(), reason="Rust solver not built")
@pytest.mark.parametrize(
    "name",
    [
        "shifter_standard",
        "pitch_box2_p10",
        "trench_both",
        "array_ramp",
        "top_liquid",
        "transient_box1",
        "solver_cg_trench",
    ],
)
def test_solver_and_reference_agree(name):
    path = ROOT / "cases" / f"{name}.heat"
    cmp = compare(path, engine.run(path))
    assert not cmp["outside"], cmp["outside"]
