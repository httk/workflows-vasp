"""The ``vasp-relax-ions`` package, driven end to end through the real manager.

Reuses the fake VASP and the job helpers of ``test_vasp_runners.py``: every test
submits a real job to a real workspace and lets a real :class:`httk.workflow.TaskManager`
run it from this repository's package directory.
"""

import json
from pathlib import Path

import pytest

from test_vasp_runners import _campaign, _failure, _job_state, _payload_of, _runs


def _history(payload: Path) -> dict:
    return json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))


def test_the_ionic_relaxation_runs_its_template_and_collects(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from httk.workflow import collect

    workspace, job_id = _campaign(tmp_path / "ions", monkeypatch, "vasp-relax-ions", publish_data=False)

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    (run,) = _runs(workdir)
    tags = run["tags"]
    assert (tags["ISIF"], tags["IBRION"], tags["NSW"]) == ("2", "2", "99")
    # The template's tight EDIFF and EDIFFG are kept, not derived.
    assert (tags["EDIFF"], tags["EDIFFG"]) == ("1E-7", "-0.001")
    assert (tags["ISMEAR"], tags["SIGMA"], tags["LREAL"], tags["ICHARG"]) == ("0", "0.05", ".FALSE.", "2")
    assert run["INCAR"].count("\nISIF") == 1
    # Without kpoint_grid the mesh is derived from kpoint_density.
    assert (workdir / "KPOINTS").read_text(encoding="utf-8").splitlines()[2:4] == ["Monkhorst-Pack", "11 11 11"]
    assert _job_state(payload)["classification"] == "completed"

    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {"relaxed_structure", "total_energy"}
    assert not item.unfulfilled
    assert item.outputs["total_energy"].value == pytest.approx(-10.5)
    assert len(item.outputs["relaxed_structure"].sites.reduced_coords) == 2


def test_isif_4_reaches_the_incar(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "isif", monkeypatch, "vasp-relax-ions", parameters={"isif": "4"})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    (run,) = _runs(payload / "run")
    assert run["tags"]["ISIF"] == "4"
    assert run["INCAR"].count("\nISIF") == 1


@pytest.mark.parametrize(
    ("parameters", "code"),
    (
        ({"incar_tags": {"ISIF": 3}}, "vasp.unsupported_incar_tag"),
        ({"incar_tags": {"nsw": 5}}, "vasp.unsupported_incar_tag"),
        ({"isif": 3}, "vasp.input_invalid"),
        ({"kpoint_grid": "4 4"}, "vasp.input_invalid"),
        ({"excluded_remedies": "bump_kpoint"}, "vasp.input_invalid"),
    ),
)
def test_a_refused_parameter_fails_before_any_vasp_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict[str, object], code: str
) -> None:
    workspace, job_id = _campaign(tmp_path / "refused", monkeypatch, "vasp-relax-ions", parameters=parameters)

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == code
    if code == "vasp.unsupported_incar_tag":
        assert failure["details"]["source"] == "incar_tags"
        assert failure["details"]["tag"] in ("ISIF", "NSW")
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()


def test_kpoint_grid_fixes_the_mesh(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "grid", monkeypatch, "vasp-relax-ions", parameters={"kpoint_grid": "4 4 4", "centering": "Gamma"}
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    assert (payload / "run" / "KPOINTS").read_text(encoding="utf-8").splitlines()[1:5] == [
        "0",
        "Gamma",
        "4 4 4",
        "0 0 0",
    ]


@pytest.mark.parametrize(("excluded", "step"), (("", 1), ("bump_kpoints", 2)))
def test_the_ladder_never_scales_the_cell_and_skips_excluded_rungs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, excluded: str, step: int
) -> None:
    # The zpotrf ladder: scale_lattice, bump_kpoints, bump_bands, ...
    workspace, job_id = _campaign(
        tmp_path / "zpotrf",
        monkeypatch,
        "vasp-relax-ions",
        fail_once=True,
        parameters={"excluded_remedies": excluded},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    history = _history(payload)
    # The ladder records the next rung to try.
    assert history["attempts"] == {"zpotrf": step + 1}
    assert history["events"][0]["step"] == step
    assert (workdir / "POSCAR").read_text(encoding="utf-8").splitlines()[1] == "1.0"
    first, second = _runs(workdir)
    if step == 2:
        # The first run had VASP's default: max(floor((8 + 2) / 2) + 3, floor(0.6 * 8)) = 8 for two Si (8
        # electrons), plus floor((2 + 1) / 2) = 1 for ISPIN = 2 with the default moments, so 9.
        assert "NBANDS" not in first["tags"]
        assert int(second["tags"]["NBANDS"]) > 9
        assert (workdir / "KPOINTS").read_text(encoding="utf-8").splitlines()[3] == "11 11 11"
    else:
        assert (workdir / "KPOINTS").read_text(encoding="utf-8").splitlines()[3] != "11 11 11"


def test_ionic_nonconvergence_continues_from_the_contcar(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "ionic", monkeypatch, "vasp-relax-ions", fake={"ionic": (0,)})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    assert len(_runs(workdir)) == 2
    history = _history(payload)
    assert history["attempts"] == {"ionic_nonconvergence": 1}
    assert history["events"][0]["files"][0]["path"] == "POSCAR"
    assert _job_state(payload)["remedies"] == 1


def test_a_second_ionic_nonconvergence_has_no_remaining_remedy(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "ionic", monkeypatch, "vasp-relax-ions", fake={"ionic": (0, 1)})

    assert _payload_of(workspace, job_id)[0] == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.failed"
    assert failure["message"] == "VASP nonconverged with no remaining remedy"
    assert failure["details"]["problem"] == "ionic_nonconvergence"
