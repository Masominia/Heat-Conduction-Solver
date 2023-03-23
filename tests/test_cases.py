import json

import pytest

from heatphase import cases, engine
from heatphase.engine import ROOT

ALL = cases.discover(ROOT / "cases")


def test_library_is_present_and_described():
    assert len(ALL) >= 75
    for path in ALL:
        assert len(cases.notes(path)) >= 2, f"{path.name} has no description"
        assert "expect " in path.read_text(), f"{path.name} checks nothing"


def test_families():
    fam = {"crosstalk": ["pitch_"], "numerics": ["mesh_", "solver_"]}
    assert cases.family_of("pitch_box2_p10", fam) == "crosstalk"
    assert cases.family_of("solver_cg_array", fam) == "numerics"
    assert cases.family_of("metal_via", fam) == "other"


@pytest.mark.skipif(not engine.available(), reason="Rust solver not built")
@pytest.mark.parametrize("path", ALL, ids=lambda p: p.stem)
def test_every_section_meets_its_expectations(path):
    result = engine.run(path, "check")
    assert result["checks"], "no expectations evaluated"
    failed = [c for c in result["checks"] if c["status"] != "PASS"]
    assert not failed, failed


@pytest.mark.skipif(not engine.available(), reason="Rust solver not built")
def test_record_round_trip(tmp_path):
    record, result = cases.run_case(ROOT / "cases" / "shifter_standard.heat")
    back = json.loads(cases.write_record(record, tmp_path).read_text())
    assert back["case"] == "shifter_standard" and 30 < back["metrics"]["peak"] < 50
    assert "field" not in back and result["field"].shape == (60, 240) == result["conductivity"].shape
