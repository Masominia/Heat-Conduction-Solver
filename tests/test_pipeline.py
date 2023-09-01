import json

import pytest

from heatphase import engine
from heatphase.pipeline import run


@pytest.mark.skipif(not engine.available(), reason="Rust solver not built")
def test_fast_pipeline(tmp_path):
    m = run(out=tmp_path, fast=True)
    assert not [c for c in m["checks"] if c["status"] == "FAIL"]
    assert m["sections"] >= 25 and m["largest_differences"]["field"] < 2e-6
    assert 2.5 < m["studies"]["pitch_box2"]["decay_um"] < 5 and len(m["studies"]["oxide_tradeoff"]) == 3
    assert (tmp_path / "report.md").exists() and (tmp_path / "figures" / "overview.png").exists()
    records = list((tmp_path / "cases").glob("*.json"))
    assert len(records) == m["sections"] == len(list((tmp_path / "figures" / "cases").glob("*.png")))
    assert "field" in json.loads(records[0].read_text())["against_reference"]
