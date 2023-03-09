import pytest

from heatphase import engine

pytestmark = pytest.mark.skipif(not engine.available(), reason="Rust solver not built")

BASE = (
    "grid dx=0.25 dy=0.25\nwindow width=40 ymin=-5 ymax=5\nlayer si -5 0\nrect si -0.25 0.25 2 2.25\n"
    "heater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.1\nprobe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\nprobe far x=8.125 y=2.125\n"
)


def write(tmp_path, text):
    path = tmp_path / "section.heat"
    path.write_text(text)
    return path


def test_run_returns_fields_and_traces(tmp_path):
    r = engine.run(write(tmp_path, "title Test\n" + BASE + "transient heater=h1 dt=1 steps=100\n"))
    assert r["title"] == "Test"
    assert r["field"].shape == (40, 160) == r["conductivity"].shape
    assert r["field"].max() == pytest.approx(r["metrics"]["peak"], rel=1e-6)
    assert set(r["conductivity"].ravel()) == {1.38, 148.0}
    assert r["temperature"]["wg"] == pytest.approx(0.1 * r["response"]["wg"]["h1"])
    assert r["pi_power"]["wg"]["h1"] == pytest.approx(1.55 / (2 * 1.86e-4 * r["response"]["wg"]["h1"]))
    assert len(r["time"]) == 100 and set(r["traces"]) == {"wg", "far"}
    assert r["step"]["far"]["tau"] > r["step"]["wg"]["tau"] > 0


def test_check_reports_pass_and_fail(tmp_path):
    r = engine.run(write(tmp_path, BASE + "expect hotter wg far\nexpect hotter far wg\n"), "check")
    assert [c["status"] for c in r["checks"]] == ["PASS", "FAIL"]
    assert "field" not in r


def test_bad_file_raises_with_the_line(tmp_path):
    with pytest.raises(RuntimeError, match="line 2"):
        engine.run(write(tmp_path, "grid dx=0.5\nnonsense 1\n"))
    with pytest.raises(RuntimeError, match="no heater"):
        engine.run(write(tmp_path, "window width=40 ymin=-5 ymax=5\n"))
    with pytest.raises(RuntimeError, match="cannot read"):
        engine.run(tmp_path / "missing.heat")
