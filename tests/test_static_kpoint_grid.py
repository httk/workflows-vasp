"""The ``kpoint_grid`` and ``excluded_remedies`` parameters of ``vasp-static``.

Reuses the fake VASP and the job helpers of ``test_vasp_runners.py``.
"""

import json
from pathlib import Path

import pytest

from test_vasp_runners import _campaign, _failure, _payload_of, _runs


def test_kpoint_grid_fixes_the_static_mesh(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "grid", monkeypatch, "vasp-static", parameters={"kpoint_grid": [3, 5, 7]})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    assert (payload / "run" / "KPOINTS").read_text(encoding="utf-8").splitlines()[1:5] == [
        "0",
        "Monkhorst-Pack",
        "3 5 7",
        "0 0 0",
    ]


@pytest.mark.parametrize("grid", ("4 4 0", "4 4 4.5", [4, 4, True], "a b c"))
def test_an_invalid_kpoint_grid_fails_before_any_vasp_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, grid: object
) -> None:
    workspace, job_id = _campaign(tmp_path / "bad", monkeypatch, "vasp-static", parameters={"kpoint_grid": grid})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_invalid"
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()


def test_excluded_remedies_skip_a_static_rung(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # The zpotrf ladder starts with scale_lattice; excluded, bump_kpoints is next.
    workspace, job_id = _campaign(
        tmp_path / "zpotrf",
        monkeypatch,
        "vasp-static",
        fail_once=True,
        parameters={"excluded_remedies": "scale_lattice incar.LREAL"},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    # The ladder records the next rung to try.
    assert history["attempts"] == {"zpotrf": 2}
    assert history["events"][0]["step"] == 1
    assert (workdir / "POSCAR").read_text(encoding="utf-8").splitlines()[1] == "1.0"
    assert len(_runs(workdir)) == 2
