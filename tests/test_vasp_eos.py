"""The ``vasp-eos`` composite, end to end through the real manager, against a model VASP.

The model VASP reports the energy of a third-order Birch–Murnaghan curve of known
parameters at the volume of its POSCAR, in an OUTCAR the real ``OutcarFile`` parses;
a single point (``NSW = 0``) reports it shifted by ``static_shift``, so a test sees
which energies were fitted. A ``vasp.relax`` stage (``ISIF = 3`` or ``7``) relaxes the
cell isotropically to the model's minimum; every other run keeps its POSCAR as CONTCAR.
The composite calls ``vasp-relax``, ``vasp-relax-ions`` and ``vasp-static`` as child
jobs, so those packages are installed into the test workspace from this repository.
"""

import json
import runpy
from pathlib import Path
from typing import Any

import pytest
from httk.workflow import TaskManager, Workspace, _kernel, _store
from httk.workflow.introspection import list_jobs, resolve_job
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POSCAR, _POTCAR

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-eos"
CALLEES = ("vasp-relax", "vasp-relax-ions", "vasp-static")
EXCLUDED = "bump_kpoints equal_kpoints centering incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice"
GPA_PER_EV_PER_A3 = 160.21766208
V_REF = 8.0  # the volume of _POSCAR
_MODEL = {
    "v0": V_REF * 1.01,
    "e0": -10.8,
    "b0": 100.0 / GPA_PER_EV_PER_A3,
    "b1": 4.5,
    "static_shift": -0.01,
    "fail_below": 0.0,
}

_FAKE_VASP = '''#!/usr/bin/env python3
"""A model VASP: the energy of a Birch–Murnaghan curve at the volume of POSCAR."""

import json
from pathlib import Path

MODEL = json.loads(MODEL_JSON)

tags = {}
for line in Path("INCAR").read_text().splitlines():
    for part in line.split("#")[0].split("!")[0].split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            tags[key.strip().upper()] = value.strip()
poscar = Path("POSCAR").read_text()
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps({"INCAR": Path("INCAR").read_text(), "tags": tags, "POSCAR": poscar}) + "\\n")
lines = poscar.splitlines()
scale = float(lines[1].split()[0])
m = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
volume = abs(
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
)
if tags.get("ISIF") == "4" and volume < MODEL["fail_below"]:
    print("LAPACK: Routine ZPOTRF failed!")
    Path("OUTCAR").write_text(" fake vasp 6.4.1\\n   NELM   =     60\\n   NSW    =     99\\n")
    raise SystemExit(1)
x = (MODEL["v0"] / volume) ** (2 / 3) - 1
energy = MODEL["e0"] + 9 * MODEL["v0"] * MODEL["b0"] / 16 * (2 * x * x + (MODEL["b1"] - 4) * x**3)
nsw = int(tags.get("NSW", "0"))
if nsw == 0:
    energy += MODEL["static_shift"]
ions = sum(int(n) for n in lines[6].split())
positions = "".join("      0.00000      0.00000      0.00000         0.000000      0.000000      0.000000\\n" for _ in range(ions))
Path("OUTCAR").write_text(
    " vasp.6.4.1 fake\\n"
    "   ions per type =               %d\\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\\n"
    "   NSW    =     %d    number of steps for IOM\\n"
    "--------------------------------------- Iteration      1(   1)  ---------------------------------------\\n"
    "  FORCE on cell =-STRESS in cart. coord.  units (eV):\\n"
    "  in kB       0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\\n"
    "  external pressure =        0.00 kB  Pullay stress =        0.00 kB\\n"
    " POSITION                                       TOTAL-FORCE (eV/Angst)\\n"
    " -----------------------------------------------------------------------------------\\n"
    "%s"
    " -----------------------------------------------------------------------------------\\n"
    "    total drift:                                0.000000      0.000000      0.000000\\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\\n"
    "  ---------------------------------------------------\\n"
    "  free  energy   TOTEN  =       %.8f eV\\n"
    "  energy  without entropy=      %.8f  energy(sigma->0) =      %.8f\\n"
    " General timing and accounting informations for this job:\\n" % (ions, nsw, positions, energy, energy, energy)
)
Path("OSZICAR").write_text(
    "       N       E                     dE             d eps       ncg     rms\\n"
    "DAV:   1    %.12E   -0.50000E+00   -0.10000E+00   128   0.100E+00\\n"
    "   1 F= %.8E E0= %.8E  d E =-.105000E+02\\n" % (energy, energy, energy)
)
contcar = list(lines)
if tags.get("ISIF") in ("3", "7"):
    factor = (MODEL["v0"] / volume) ** (1 / 3)
    for i in (2, 3, 4):
        contcar[i] = "  " + "  ".join("%.16f" % (factor * x) for x in m[i - 2])
Path("CONTCAR").write_text("\\n".join(contcar) + "\\n")
Path("vasprun.xml").write_text('<modeling><structure name="finalpos"></structure></modeling>\\n')
if tags.get("LCHARG", "").strip(".").upper().startswith("T"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''

_CARTESIAN = _POSCAR.replace("Direct", "Cartesian")


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    parameters: dict[str, object],
    *,
    model: dict[str, float] | None = None,
    poscar: str = _POSCAR,
) -> tuple[Workspace, str]:
    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    model_json = repr(json.dumps({**_MODEL, **(model or {})}))
    executable.write_text(_FAKE_VASP.replace("MODEL_JSON", model_json), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    files = {}
    for name, content in (("POSCAR", poscar), ("POTCAR", _POTCAR)):
        (staged / name).write_text(content, encoding="utf-8")
        files[name] = staged / name
    # The callees first, from this repository, so the composite's declared calls resolve to them.
    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        for callee in CALLEES:
            _store.install(workspace, owner, REPO_ROOT / callee)
    finally:
        owner.close()
    job = new_job(workspace, PACKAGE, files=files, parameters=parameters, name="packaged vasp-eos", install=True)
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=600.0)
    return workspace, job.job_id


def _job(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _runs(payload: Path) -> list[dict[str, Any]]:
    path = payload / "run" / "fake-vasp-runs.jsonl"
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _eos(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "run" / "eos.json").read_text(encoding="utf-8"))


def _volume(poscar: str) -> float:
    lines = poscar.splitlines()
    a, b, c = ([float(x) for x in lines[i].split()] for i in (2, 3, 4))
    return abs(
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _state_failure(workspace: Workspace, job_id: str) -> dict[str, Any]:
    state, payload = _job(workspace, job_id)
    assert state == "failed"
    return json.loads((payload / "state.json").read_text(encoding="utf-8"))["failure"]


@pytest.mark.parametrize(("relax", "static"), ((False, False), (False, True), (True, False)))
def test_the_fit_recovers_the_model_equation_of_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, test_profile, relax: bool, static: bool
) -> None:
    from httk.workflow import collect

    frames = 6 if test_profile.extended else 5
    # The static case also publishes, so collect reads the published copies there.
    parameters = {"relax": relax, "static": static, "number_of_frames": frames, "publish_data": static}
    workspace, job_id = _campaign(tmp_path / "eos", monkeypatch, parameters)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    eos = _eos(payload)
    shift = _MODEL["static_shift"] if static else 0.0
    fit = eos["fit"]
    assert fit["equilibrium_volume"] == pytest.approx(_MODEL["v0"], rel=1e-4)
    assert fit["equilibrium_energy"] == pytest.approx(_MODEL["e0"] + shift, abs=1e-5)
    assert fit["bulk_modulus"] == pytest.approx(100.0, abs=0.1)
    assert fit["bulk_modulus_derivative"] == pytest.approx(4.5, abs=0.05)
    assert eos["energy_source"] == ("static" if static else "relax_ions")

    # The reference: the input structure, or the tight relaxation's CONTCAR (at the model's minimum).
    reference = (payload / "run" / "relaxed" / "POSCAR").read_text(encoding="utf-8")
    if relax:
        _, relax_payload = _job(workspace, eos["children"]["relax"]["job_id"])
        assert reference == (relax_payload / "run" / "CONTCAR").read_text(encoding="utf-8")
        tight = [run["tags"] for run in _runs(relax_payload) if run["tags"].get("IBRION") == "1"]
        assert tight and all(
            (t["EDIFF"], t["EDIFFG"], t["LREAL"], t["ENCUT"]) == ("1E-7", "-0.001", ".FALSE.", "600") for t in tight
        )
        assert _volume(reference) == pytest.approx(_MODEL["v0"], rel=1e-12)
    else:
        assert reference == _POSCAR and "relax" not in eos["children"]
    v_ref = _volume(reference)
    assert eos["reference_volume"] == pytest.approx(v_ref, rel=1e-12)

    # Every volume child: ISIF = 4 on one mesh with the excluded remedies, at V_ref (1 + e)^3.
    strains = [-0.05 + 0.1 * i / (frames - 1) for i in range(frames)]
    assert [point["linear_strain"] for point in eos["points"]] == pytest.approx(strains)
    meshes = set()
    for point, strain in zip(eos["points"], strains, strict=True):
        assert point["volume"] == pytest.approx(v_ref * (1 + strain) ** 3, rel=1e-12)
        _, child = _job(workspace, eos["children"][point["label"]]["job_id"])
        (run,) = _runs(child)
        assert run["tags"]["ISIF"] == "4"
        assert run["POSCAR"].splitlines()[5:] == reference.splitlines()[5:]
        job = json.loads((child / "job.json").read_text(encoding="utf-8"))
        assert job["parameters"]["excluded_remedies"] == EXCLUDED
        assert job["parameters"]["kpoint_grid"] == eos["kpoint_grid"]
        meshes.add((child / "run" / "KPOINTS").read_text(encoding="utf-8"))
        if static:
            _, single = _job(workspace, eos["children"][point["static_label"]]["job_id"])
            (static_run,) = _runs(single)
            assert static_run["tags"]["NSW"] == "0"
            assert static_run["POSCAR"] == (child / "run" / "CONTCAR").read_text(encoding="utf-8")
            meshes.add((single / "run" / "KPOINTS").read_text(encoding="utf-8"))
        else:
            assert "static_label" not in point
    assert len(meshes) == 1 and eos["kpoint_grid"] in meshes.pop()
    assert len(list_jobs(workspace)) == 1 + frames * (2 if static else 1) + relax

    items = {item.record.job_id: item for item in collect(workspace, fail_fast=True, allow_job_collector=True)}
    item = items[job_id]
    assert set(item.outputs) == {
        "relaxed_structure",
        "equilibrium_volume",
        "equilibrium_energy",
        "bulk_modulus",
        "bulk_modulus_pressure_derivative",
        "eos_result",
    }
    assert not item.unfulfilled
    assert item.outputs["bulk_modulus"].value == pytest.approx(fit["bulk_modulus"])
    assert item.outputs["eos_result"].name == "eos.json"
    assert ("/data/" in item.outputs["eos_result"].url) == static


def test_a_failed_volume_child_fails_the_job_by_its_label(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Only the smallest volume (strain -5 %) fails, and with no remedy budget it stays failed.
    workspace, job_id = _campaign(
        tmp_path / "failed",
        monkeypatch,
        {"relax": False, "number_of_frames": 5, "maximum_remedies": 0},
        model={"fail_below": V_REF * 0.9},
    )
    failure = _state_failure(workspace, job_id)
    assert failure["code"] == "vasp.dependency_failed"
    assert "v00" in failure["message"]
    (failed,) = failure["details"]["failed"]
    assert failed["label"] == "v00" and failed["failure"]["code"] == "vasp.failed"


@pytest.mark.parametrize(
    ("parameters", "poscar", "code"),
    (
        ({"number_of_frames": 3}, _POSCAR, "vasp.input_invalid"),
        ({"linear_strain": [0.01, 0.05]}, _POSCAR, "vasp.input_invalid"),
        ({"relax": "perhaps"}, _POSCAR, "vasp.input_invalid"),
        ({"relax": False}, _CARTESIAN, "vasp.input_unsupported"),
    ),
)
def test_refused_inputs_fail_before_any_call(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict[str, object], poscar: str, code: str
) -> None:
    workspace, job_id = _campaign(tmp_path / "refused", monkeypatch, parameters, poscar=poscar)
    assert _state_failure(workspace, job_id)["code"] == code
    assert len(list_jobs(workspace)) == 1


def test_the_consistency_check_ignores_only_the_electronic_path_tags() -> None:
    inconsistency = runpy.run_path(str(PACKAGE / "run.py"))["inconsistency"]
    mesh = "mesh\n0\nGamma\n4 4 4\n0 0 0\n"
    a = {"label": "a", "kpoints": mesh, "incar": {"ENCUT": "600", "ALGO": "Normal"}, "remedies": []}
    remedied = {"label": "b", "kpoints": mesh, "incar": {"ENCUT": "600", "ALGO": "All"}, "remedies": ["incar.ALGO"]}
    assert inconsistency([a, remedied]) is None
    # scale_ediff tightens EDIFF and the dentet rung raises NEDOS: neither changes the energy surface.
    tightened = {"label": "e", "kpoints": mesh, "incar": {**a["incar"], "EDIFF": "1E-8", "NEDOS": "2001"}}
    assert inconsistency([a | {"incar": {**a["incar"], "EDIFF": "1E-7"}}, tightened]) is None
    other_cutoff = {"label": "c", "kpoints": mesh, "incar": {"ENCUT": "520", "ALGO": "Normal"}}
    assert inconsistency([a, other_cutoff]) == {"ENCUT": {"a": "600", "c": "520"}}
    other_mesh = {"label": "d", "kpoints": mesh.replace("4 4 4", "5 5 5"), "incar": a["incar"]}
    assert set(inconsistency([a, other_mesh])) == {"KPOINTS"}
