"""The ``vasp-relax-static-hse`` composer, end to end through the real manager, against a model VASP.

The model VASP is that of ``test_vasp_relax_static_coarse.py`` with one change:
a run whose INCAR sets ``LHFCALC`` sees a target metric scaled by
``hse_metric_scale``, so the HSE surface has its minimum at a lattice about 1 %
off the PBE one. The composer calls ``vasp-relax-static-coarse`` as child jobs,
so both packages are installed into the test workspace from this repository.
"""

import json
import math
import time
from pathlib import Path
from typing import Any

import pytest
from httk.workflow import TaskManager, Workspace, _kernel, _store
from httk.workflow.introspection import list_jobs, resolve_job
from httk.workflow.scaffold import new_job

from test_vasp_relax_static_coarse import _FAKE_VASP, _MODEL, _POSCAR, _POTCAR

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-relax-static-hse"
CALLEE = REPO_ROOT / "vasp-relax-static-coarse"

_GSTAR_LINE = 'gstar = [[g * (MODEL["final_metric_scale"] if final else 1.0) for g in row] for row in MODEL["gstar"]]\n'
_HSE_GSTAR_LINE = (
    'hse = tags.get("LHFCALC", "").strip(".").upper().startswith("T")\n'
    'gstar = [[g * (MODEL["final_metric_scale"] if final else 1.0) * (MODEL["hse_metric_scale"] if hse else 1.0)'
    ' for g in row] for row in MODEL["gstar"]]\n'
)
assert _FAKE_VASP.count(_GSTAR_LINE) == 1
_FAKE_HSE_VASP = _FAKE_VASP.replace(_GSTAR_LINE, _HSE_GSTAR_LINE)

# The HSE metric 2 % smaller: lattice lengths about 1 % shorter than PBE's.
_HSE_MODEL = {**_MODEL, "hse_metric_scale": 0.98}
_MESHES = {"coarse_mesh": "3 3 3", "final_mesh": "6 6 6"}


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    parameters: dict[str, object] | None = None,
    sleep: float | None = None,
    end_time_after: float | None = None,
) -> tuple[Workspace, str]:
    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    executable.write_text(_FAKE_HSE_VASP.replace("MODEL_JSON", repr(json.dumps(_HSE_MODEL))), encoding="utf-8")
    executable.chmod(0o755)
    if sleep is not None:
        (root / "fake-sleep").write_text(str(sleep), encoding="utf-8")
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    files = {}
    for name, content in (("POSCAR", _POSCAR), ("POTCAR", _POTCAR)):
        (staged / name).write_text(content, encoding="utf-8")
        files[name] = staged / name
    # The callee first, from this repository, so the composer's declared call resolves to it.
    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        _store.install(workspace, owner, CALLEE)
    finally:
        owner.close()
    job = new_job(
        workspace,
        PACKAGE,
        files=files,
        parameters={"publish_data": True, "centering": "Gamma", **(parameters or {})},
        name="packaged vasp-relax-static-hse",
        install=True,
    )
    _run(workspace, end_time_after)
    return workspace, job.job_id


def _run(workspace: Workspace, end_time_after: float | None = None) -> None:
    end_time = None if end_time_after is None else time.time() + end_time_after
    with TaskManager(workspace, heartbeat_interval=0.01, end_time=end_time, deadline_margin=0.0) as manager:
        manager.run_until_idle(timeout=600.0)


def _job(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _failure(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "state.json").read_text(encoding="utf-8"))["failure"]


def _stages(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "run" / "stages.json").read_text(encoding="utf-8"))


def _child(workspace: Workspace, stages: dict[str, Any], label: str) -> Path:
    state, path = _job(workspace, stages["stages"][label]["job_id"])
    assert state == "succeeded"
    return path


def _runs(workdir: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in (workdir / "fake-vasp-runs.jsonl").read_text(encoding="utf-8").splitlines()]


def _report(workdir: Path) -> dict[str, Any]:
    return json.loads((workdir / "relax-driver" / "report.json").read_text(encoding="utf-8"))


def _lengths(poscar: str) -> list[float]:
    lines = poscar.splitlines()
    return [math.sqrt(sum(float(x) ** 2 for x in lines[i].split())) for i in (2, 3, 4)]


def test_without_pbe_the_composer_relaxes_pbe_first_and_hse_from_its_result(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The meshes are given to the PBE stage only: the HSE stage must take them from the PBE result.
    workspace, job_id = _campaign(tmp_path / "chain", monkeypatch, parameters={"pbe_parameters": _MESHES})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    stages = _stages(payload)
    assert stages["mode"] == "no_pbe"
    pbe, hse = _child(workspace, stages, "pbe"), _child(workspace, stages, "hse")
    pbe_runs, hse_runs = _runs(pbe / "run"), _runs(hse / "run")
    pbe_report, hse_report = _report(pbe / "run"), _report(hse / "run")

    # Each stage is one vasp.relax-static-coarse job: PBE without LHFCALC, HSE with it in every evaluation.
    for record in (stages["stages"]["pbe"], stages["stages"]["hse"]):
        assert record["workflow"]["name"] == "vasp.relax-static-coarse"
        assert record["within_limits"] is True
    assert not any("LHFCALC" in run["tags"] for run in pbe_runs)
    assert all(run["tags"]["LHFCALC"] == ".TRUE." and run["tags"]["ALGO"] == "Damped" for run in hse_runs)

    # The HSE stage starts at the PBE stage's relaxed geometry, on the PBE meshes.
    assert hse_runs[0]["POSCAR"] == (pbe / "run" / "relax" / "CONTCAR").read_text(encoding="utf-8")
    assert hse_report["meshes"] == pbe_report["meshes"] == {"coarse": ["3", "3", "3"], "final": ["6", "6", "6"]}

    # Seeded: the first HSE coarse and final evaluations restart from the PBE WAVECARs of their meshes.
    first = {}
    for record in hse_report["evaluations"]:
        first.setdefault(record["kind"], record)
    assert first["coarse"]["restart_source"] == "seed" and first["coarse"]["restart_observed"] is True
    assert first["coarse"]["scf_steps"] == 6
    assert first["final"]["restart_source"] == "seed" and first["final"]["restart_observed"] is True
    assert set(stages["stages"]["hse"]["seeds"]) == {"coarse_wavecar", "final_wavecar"}
    # Off by default: the HSE stage starts from the isotropic Hessian model.
    assert hse_report["parameters"]["hessian_seed"] == ""

    # The geometry shift is the model's: the HSE lattice about 1 % shorter than the PBE one.
    shift = stages["stages"]["hse"]["shift_from_start"]
    assert all(-1.3 < x < -0.7 for x in shift["lengths_percent"])
    assert -3.5 < shift["volume_percent"] < -2.5

    # The composer's workdir holds the HSE result in the driver's layout, which collect.py reads.
    workdir = payload / "run"
    assert (workdir / "relax" / "CONTCAR").read_text() == (hse / "run" / "relax" / "CONTCAR").read_text()
    assert (workdir / "OUTCAR").read_text() == (hse / "run" / "OUTCAR").read_text()
    assert (payload / "data" / "static" / "OUTCAR").is_file() and (payload / "data" / "stages.json").is_file()
    state = json.loads((payload / ".httk-job" / "state.json").read_text(encoding="utf-8"))
    assert state["within_limits"] is True
    assert state["static_energy"] == pytest.approx(stages["stages"]["hse"]["final_energy"])
    assert _lengths((workdir / "relax" / "CONTCAR").read_text())[0] < _lengths(hse_runs[0]["POSCAR"])[0]


def test_from_a_pbe_relaxed_structure_only_the_hse_stage_runs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "given",
        monkeypatch,
        parameters={"mode": "start_from_pbe", "pbe_provenance": "test:J-0000", **_MESHES},
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    stages = _stages(payload)
    assert stages["stages"]["pbe"]["source"] == "given" and stages["stages"]["pbe"]["provenance"] == "test:J-0000"
    hse = _child(workspace, stages, "hse")
    runs = _runs(hse / "run")
    # The HSE stage starts at the given structure itself.
    assert _lengths(runs[0]["POSCAR"]) == pytest.approx(_lengths(_POSCAR), abs=1e-9)
    # No seeds were given: every first evaluation is a fresh start.
    assert stages["stages"]["hse"]["seeds"] == {}
    assert runs[0]["tags"]["ISTART"] == "0"


def test_an_allocation_that_ends_mid_chain_is_resumed_by_a_later_manager(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Every model run takes 2 s and the allocation ends after 9 s: the PBE stage's time guard defers its next
    # evaluation, and the composer waits on it without holding a worker.
    workspace, job_id = _campaign(
        tmp_path / "split", monkeypatch, parameters={"pbe_parameters": _MESHES}, sleep=2.0, end_time_after=9.0
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "waiting"
    others = [row["job_id"] for row in list_jobs(workspace) if row["job_id"] != job_id]
    assert len(others) == 1
    pbe_kind, pbe = _job(workspace, others[0])
    assert pbe_kind == "ready"
    assert "mintime" in json.loads((pbe / "state.json").read_text(encoding="utf-8"))["resources"]

    # A manager without an allocation end finishes the PBE stage, calls the HSE stage and completes the composer.
    (tmp_path / "split" / "fake-sleep").unlink()
    _run(workspace)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    stages = _stages(payload)
    assert stages["stages"]["pbe"]["within_limits"] is True and stages["stages"]["hse"]["within_limits"] is True


def test_a_pbe_relaxed_start_without_provenance_fails_by_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "noprov", monkeypatch, parameters={"mode": "start_from_pbe", **_MESHES})
    kind, payload = _job(workspace, job_id)
    assert kind == "failed"
    assert _failure(payload)["code"] == "vasp.pbe_provenance_missing"


def test_the_pbe_hessian_seeds_the_hse_stage_when_asked(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "hessian", monkeypatch, parameters={"hessian_from_pbe": True, **_MESHES})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    stages = _stages(payload)
    pbe, hse = _child(workspace, stages, "pbe"), _child(workspace, stages, "hse")
    count = stages["stages"]["pbe"]["evaluations"]
    snapshot = json.loads((pbe / "run" / "relax-driver" / f"driver-{count:03d}.json").read_text(encoding="utf-8"))
    seed = json.loads((hse / "files" / "hessian.pbe.json").read_text(encoding="utf-8"))
    assert seed["hessian"] == snapshot["optimizer"]["B"]
    assert _report(hse / "run")["parameters"]["hessian_seed"] == "files/hessian.pbe.json"
    assert "hessian_seed" in stages["stages"]["hse"]["seeds"]


def test_a_failed_pbe_stage_fails_the_job_by_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "pbefail", monkeypatch, parameters={"pbe_parameters": {**_MESHES, "maximum_evaluations": 3}}
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "failed"
    failure = _failure(payload)
    assert failure["code"] == "vasp.pbe_stage_failed"
    # The child's own failure travels with it, and no HSE stage was called.
    stages = _stages(payload)
    assert stages["stages"]["pbe"]["failure"]["code"] == "vasp.relax_not_converged"
    assert stages["stages"]["pbe"]["kind"] == "failed"
    assert "hse" not in stages["stages"]
