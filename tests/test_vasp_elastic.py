"""The ``vasp-elastic`` composite, end to end through the real manager, against a model VASP.

The model VASP answers a strained cell with the stress of a known cubic stiffness
``C`` plus a small offset: it reads the lattice ``L'`` of its POSCAR, recovers
``F^T = L^-1 L'`` from the reference lattice ``L`` baked into it, the small strain
``eps = sym(F) - I``, and writes ``sigma = C e + sigma_0`` as the OUTCAR ``in kB``
row (VASP order, compressive-positive kbar) of a frame the real ``OutcarFile``
parses. Every child is a real job: ``vasp-relax`` and ``vasp-relax-ions`` are
installed into the test workspace from this repository before the parent.
"""

import json
import runpy
import time
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("httk.analyse")

from httk.workflow import TaskManager, Workspace, _kernel, _store
from httk.workflow.introspection import list_jobs, resolve_job
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POTCAR, _failure, _runs

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-elastic"

_LATTICE = [[3.0, 0.0, 0.0], [0.2, 3.1, 0.0], [0.1, 0.3, 2.9]]
_POSCAR = (
    "silicon, sheared\n1.0\n"
    + "".join(" ".join(f"{x:.10f}" for x in row) + "\n" for row in _LATTICE)
    + "Si\n2\nDirect\n0.0000000000 0.0000000000 0.0000000000\n0.2500000000 0.2600000000 0.2700000000\n"
)
_C11, _C12, _C44 = 250.0, 120.0, 80.0
_STIFFNESS = [
    [_C11 if i == j else _C12 if i < 3 and j < 3 else _C44 if i == j else 0.0 for j in range(6)] for i in range(6)
]
_OFFSET = [0.1, -0.2, 0.05, 0.01, 0.0, -0.02]

_FAKE_VASP = '''#!/usr/bin/env python3
"""A model VASP: the stress of a known stiffness for the strain of the POSCAR lattice against a reference."""

import json
from pathlib import Path

import numpy as np
from httk.codes.vasp import read_incar

MODEL = json.loads(MODEL_JSON)
tags = read_incar("INCAR")
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(dict(INCAR=Path("INCAR").read_text(), tags=tags)) + "\\n")
poscar = Path("POSCAR").read_text().splitlines()
lattice = np.array([[float(x) for x in poscar[i].split()[:3]] for i in (2, 3, 4)])
# Rows are the lattice vectors, so L' = L F^T.
gradient = np.linalg.solve(np.array(MODEL["reference"]), lattice).T
eps = (gradient + gradient.T) / 2 - np.eye(3)
strain = [eps[0, 0], eps[1, 1], eps[2, 2], 2 * eps[1, 2], 2 * eps[0, 2], 2 * eps[0, 1]]
if MODEL["fail_xx_below"] is not None and strain[0] < MODEL["fail_xx_below"]:
    print("the model VASP refuses this cell")
    Path("OUTCAR").write_text(" vasp.6.4.1 fake\\n")
    raise SystemExit(1)
sigma = np.array(MODEL["stiffness"]) @ strain + np.array(MODEL["offset"])
kbar = [-10 * sigma[i] for i in (0, 1, 2, 5, 3, 4)]
ions = sum(int(n) for n in poscar[6].split())
energy = -10.5
outcar = [
    " vasp.6.4.1 fake",
    "   ions per type =               %d" % ions,
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5",
    "   NSW    =     %d    number of steps for IOM" % int(tags.get("NSW", "0")),
    "--------------------------------------- Iteration      1(   1)  ---------------------------------------",
    "  FORCE on cell =-STRESS in cart. coord.  units (eV):",
    "  in kB " + " ".join("%16.8f" % value for value in kbar),
    "  external pressure =       0.00 kB  Pullay stress =        0.00 kB",
    " POSITION                                       TOTAL-FORCE (eV/Angst)",
    " -----------------------------------------------------------------------------------",
]
outcar += ["      0.00000      0.00000      0.00000         0.000000      0.000000      0.000000"] * ions
outcar += [
    " -----------------------------------------------------------------------------------",
    "    total drift:                                0.000000      0.000000      0.000000",
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)",
    "  ---------------------------------------------------",
    "  free  energy   TOTEN  =       %.8f eV" % energy,
    "  energy  without entropy=      %.8f  energy(sigma->0) =      %.8f" % (energy, energy),
    " General timing and accounting informations for this job:",
]
Path("OUTCAR").write_text("\\n".join(outcar) + "\\n")
Path("OSZICAR").write_text(
    "DAV:   1    %.12E   -0.50000E+00   -0.10000E+00   128   0.100E+00\\n"
    "   1 F= %.8E E0= %.8E  d E =-.105000E+02\\n" % (energy, energy, energy)
)
# vasp.relax stages leave the structure where it is but say who wrote it; relax-ions children copy it.
if tags.get("ISIF") != "2":
    poscar[0] = "relaxed by the fake vasp"
Path("CONTCAR").write_text("\\n".join(poscar) + "\\n")
Path("vasprun.xml").write_text("<modeling></modeling>\\n")
if tags.get("LCHARG", "").strip(".").upper().startswith("T"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    parameters: dict[str, object],
    *,
    poscar: str = _POSCAR,
    fail_xx_below: float | None = None,
) -> tuple[Workspace, str]:
    root.mkdir(parents=True)
    model = {"reference": _LATTICE, "stiffness": _STIFFNESS, "offset": _OFFSET, "fail_xx_below": fail_xx_below}
    executable = root / "fake-vasp"
    executable.write_text(_FAKE_VASP.replace("MODEL_JSON", repr(json.dumps(model))), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    files = {}
    for name, content in (("POSCAR", poscar), ("POTCAR", _POTCAR)):
        (staged / name).write_text(content, encoding="utf-8")
        files[name] = staged / name
    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        for callee in ("vasp-relax", "vasp-relax-ions"):
            _store.install(workspace, owner, REPO_ROOT / callee)
    finally:
        owner.close()
    job = new_job(workspace, PACKAGE, files=files, parameters=parameters, name="packaged vasp-elastic", install=True)
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=600.0)
    return workspace, job.job_id


def _job(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _result(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "run" / "elastic.json").read_text(encoding="utf-8"))


def _runner() -> dict[str, Any]:
    return runpy.run_path(str(PACKAGE / "run.py"))


@pytest.mark.parametrize("relax", (False, True))
def test_the_fit_recovers_the_stiffness_of_the_model(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, test_profile, relax: bool
) -> None:
    from httk.codes.vasp import read_incar
    from httk.workflow import collect

    # Normal: two magnitudes, 12 cells; extended: atomate2's default four, 24 cells.
    magnitudes = {} if test_profile.extended else {"strain_magnitudes": [-0.01, 0.01]}
    started = time.monotonic()
    workspace, job_id = _campaign(
        tmp_path / "elastic",
        monkeypatch,
        {"relax": relax, "maximum_remedies": 3, **magnitudes},
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id) if kind == "failed" else kind
    result = _result(payload)
    children = result["children"]
    assert len(children) == (24 if test_profile.extended else 12)
    assert [child["label"] for child in children] == [f"d{i:02d}" for i in range(len(children))]

    # The 5-decimal kbar of a real OUTCAR would still give 1e-3 GPa; this model writes 8.
    fit = result["fit"]
    for row, expected in zip(fit["stiffness"], _STIFFNESS, strict=True):
        assert row == pytest.approx(expected, abs=1e-3)
    assert fit["stress_offset"] == pytest.approx(_OFFSET, abs=1e-3)
    assert fit["stable"] is True
    assert fit["bulk_modulus_voigt"] == pytest.approx((_C11 + 2 * _C12) / 3, abs=1e-3)

    # One mesh, ISIF = 2 and the excluded remedies in every child; the shared parameters forwarded.
    kpoints = set()
    for child in children:
        state, child_payload = _job(workspace, child["job_id"])
        assert state == "succeeded"
        kpoints.add((child_payload / "run" / "KPOINTS").read_text(encoding="utf-8"))
        assert read_incar(child_payload / "run" / "INCAR")["ISIF"] == "2"
        parameters = json.loads((child_payload / "job.json").read_text(encoding="utf-8"))["parameters"]
        assert parameters["excluded_remedies"] == result["excluded_remedies"]
        assert parameters["kpoint_grid"] == result["kpoint_grid"]
        assert parameters["maximum_remedies"] == 3
        assert "strain_magnitudes" not in parameters
    (text,) = kpoints
    assert text.splitlines()[2:4] == ["Gamma", result["kpoint_grid"]]

    relaxed = (payload / "run" / "relaxed" / "POSCAR").read_text(encoding="utf-8")
    first = (payload / "run" / "deformations" / "d00" / "POSCAR").read_text(encoding="utf-8")
    if relax:
        # The relaxation ran its relax rounds tight from this package's INCAR.relax; the pre-stages stayed loose.
        state, relax_payload = _job(workspace, result["relax"]["job_id"])
        assert state == "succeeded"
        runs = _runs(relax_payload / "run")
        assert [run["tags"]["ISIF"] for run in runs[:2]] == ["7", "3"]
        assert all(run["tags"]["EDIFFG"] != "-0.001" for run in runs[:2])
        for run in runs[2:]:
            assert (run["tags"]["EDIFF"], run["tags"]["EDIFFG"], run["tags"]["LREAL"]) == ("1E-7", "-0.001", ".FALSE.")
        assert relaxed == (relax_payload / "run" / "CONTCAR").read_text(encoding="utf-8")
        assert first.splitlines()[0] == "relaxed by the fake vasp"
    else:
        assert result["relax"] is None
        assert relaxed == _POSCAR

    items = {item.record.job_id: item for item in collect(workspace, fail_fast=True, allow_job_collector=True)}
    item = items[job_id]
    assert set(item.outputs) == {
        "relaxed_structure",
        "elastic_result",
        "elastic_tensor",
        "compliance_tensor",
        "bulk_modulus_voigt",
        "bulk_modulus_reuss",
        "bulk_modulus_hill",
        "shear_modulus_voigt",
        "shear_modulus_reuss",
        "shear_modulus_hill",
        "universal_anisotropy_index",
    }
    assert not item.unfulfilled
    bulk = item.outputs["bulk_modulus_voigt"]
    assert bulk.bulk_modulus_voigt == pytest.approx((_C11 + 2 * _C12) / 3, abs=1e-3)
    print(f"vasp-elastic relax={relax}: {len(children)} deformations in {time.monotonic() - started:.1f} s")


def test_a_failed_deformation_fails_the_job_by_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "failed",
        monkeypatch,
        {"relax": False, "strain_magnitudes": "-0.01 0.01", "maximum_remedies": 0},
        fail_xx_below=-0.005,
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.dependency_failed"
    assert set(failure["details"]["failed"]) == {"d00"}
    assert failure["details"]["failed"]["d00"]["code"] == "vasp.failed"
    result = _result(payload)
    assert "fit" not in result
    assert [child["kind"] for child in result["children"]].count("succeeded") == 11


@pytest.mark.parametrize(
    "parameters",
    (
        {"strain_magnitudes": [0.0, 0.01]},
        {"strain_magnitudes": [0.01]},
        {"strain_states": [[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0]] * 2},
        {"strain_states": "[[1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]"},
        {"deformation_kpoint_density": 0},
        {"relax": "perhaps"},
    ),
)
def test_invalid_parameters_fail_before_any_call(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict[str, object]
) -> None:
    workspace, job_id = _campaign(tmp_path / "invalid", monkeypatch, parameters)
    assert _job(workspace, job_id)[0] == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_invalid"
    assert len(list(list_jobs(workspace))) == 1


def test_a_cartesian_structure_is_refused_before_any_call(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "cartesian", monkeypatch, {"relax": False}, poscar=_POSCAR.replace("Direct", "Cartesian")
    )
    assert _job(workspace, job_id)[0] == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_unsupported"
    assert len(list(list_jobs(workspace))) == 1


def test_a_shear_deformation_rewrites_only_the_lattice() -> None:
    runner = _runner()
    strain = [0.0, 0.0, 0.0, 0.02, -0.01, 0.006]
    gradient = runner["deformation_gradient"](strain)
    assert gradient == [[1.0, 0.003, -0.005], [0.003, 1.0, 0.01], [-0.005, 0.01, 1.0]]
    strained = runner["strained_poscar"](_POSCAR, gradient)
    old, new = _POSCAR.splitlines(keepends=True), strained.splitlines(keepends=True)
    assert [line for i, line in enumerate(new) if i not in (2, 3, 4)] == [
        line for i, line in enumerate(old) if i not in (2, 3, 4)
    ]
    # L' = L F^T: each new row is F applied to the old row.
    for i, row in enumerate(_LATTICE):
        expected = [sum(gradient[k][j] * row[j] for j in range(3)) for k in range(3)]
        assert [float(x) for x in new[2 + i].split()] == pytest.approx(expected, abs=1e-15)
    for refused in (
        _POSCAR.replace("Direct", "Cartesian"),
        _POSCAR.replace("Direct", "Selective dynamics\nDirect"),
        _POSCAR.replace("\n1.0\n", "\n-27.0\n"),
        _POSCAR.replace("Si\n2\n", "2\n"),
    ):
        with pytest.raises(ValueError):
            runner["strained_poscar"](refused, gradient)


def test_children_that_ran_with_different_settings_are_inconsistent() -> None:
    inconsistencies = _runner()["inconsistencies"]

    def child(label: str, **tags: str) -> dict[str, Any]:
        return {"label": label, "incar": {"ENCUT": "600", "ISMEAR": "0", **tags}, "kpoints": "mesh\n"}

    # A remedy may change the electronic path (ALGO, NBANDS, a tighter EDIFF, NEDOS, ...) of one child; values
    # compare normalized.
    assert inconsistencies([child("d00"), child("d01", ALGO="All", NBANDS="40"), child("d02", ENCUT=" 600 ")]) == {}
    assert inconsistencies([child("d00", EDIFF="1E-7"), child("d01", EDIFF="1E-8", NEDOS="2001")]) == {}
    differences = inconsistencies([child("d00"), child("d01", ISMEAR="-5"), child("d02", POTIM="0.1")])
    assert differences == {
        "incar_tags": {
            "ISMEAR": {"d00": "0", "d01": "-5", "d02": "0"},
            "POTIM": {"d00": "", "d01": "", "d02": "0.1"},
        }
    }
    other = {**child("d01"), "kpoints": "other\n"}
    assert inconsistencies([child("d00"), other])["kpoints"] == {"d00": "mesh\n", "d01": "other\n"}
