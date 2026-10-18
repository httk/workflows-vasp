"""The ``vasp-dielectric`` package, driven end to end through the real manager with a fake VASP.

The fake is the one of ``test_vasp_runners.py`` cut down to a single point; when ``LEPSILON`` is
true it appends the dielectric, Born charge and piezoelectric blocks in the formats
``httk.atomistic``'s OUTCAR parser accepts.
"""

import json
from pathlib import Path

import pytest

pytest.importorskip("httk.atomistic")
pytest.importorskip("httk.analyse")

from httk.workflow import TaskManager, Workspace, collect
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POSCAR, _POTCAR, _failure, _payload_of, _runs

PACKAGE = Path(__file__).resolve().parent.parent / "vasp-dielectric"

_FAKE_VASP = '''#!/usr/bin/env python3
import json
from pathlib import Path

from httk.codes.vasp import read_incar

tags = read_incar("INCAR")
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(dict(tags=tags, KPOINTS=Path("KPOINTS").read_text())) + "\\n")
text = (
    " fake vasp 6.4.1\\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\\n"
    "   NSW    =      1    number of steps for IOM\\n"
    "   number of ions     NIONS =      2\\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\\n"
    "   free  energy   TOTEN  =       -10.80000000 eV\\n"
    "   energy  without entropy=      -10.80000000  energy(sigma->0) =      -10.80000000\\n"
    " General timing and accounting informations for this job:\\n"
)
if tags.get("LEPSILON") and not {blocks_absent}:
    text += """
 MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
 ------------------------------------------------------
          10.000000     0.500000     0.000000
           0.500000    11.000000     0.000000
           0.000000     0.000000    12.000000
 ------------------------------------------------------

 MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
 ------------------------------------------------------
           2.000000     0.000000     0.000000
           0.000000     2.000000     0.000000
           0.000000     0.000000     2.000000
 ------------------------------------------------------

 BORN EFFECTIVE CHARGES (in e, cummulative output)
 -------------------------------------------------
 ion    1
    1     2.00000     0.30000     0.00000
    2     0.00000     2.00000     0.00000
    3     0.00000     0.00000     2.00000
 ion    2
    1    -2.00000     0.00000     0.00000
    2     0.00000    -2.00000     0.00000
    3     0.00000     0.00000    -2.00000

 PIEZOELECTRIC TENSOR  for field in x, y, z        (C/m^2)
  XX           YY           ZZ           XY           YZ           ZX
  --------------------------------------------------------------------------------
  x     0.10000     0.20000     0.30000     0.40000     0.50000     0.60000
  y     0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  z     0.00000     0.00000     0.00000     0.00000     0.00000     0.00000

 PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        (C/m^2)
  XX           YY           ZZ           XY           YZ           ZX
  --------------------------------------------------------------------------------
  x     1.00000     1.00000     1.00000     1.00000     1.00000     1.00000
  y     0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
  z     0.00000     0.00000     0.00000     0.00000     0.00000     0.00000
"""
Path("OUTCAR").write_text(text)
Path("OSZICAR").write_text("   1 F= -.10800000E+02 E0= -.10800000E+02  d E =-.108000E+02\\n")
Path("CONTCAR").write_text(Path("POSCAR").read_text())
Path("vasprun.xml").write_text("<modeling></modeling>\\n")
'''


def _campaign(
    root: Path, monkeypatch: pytest.MonkeyPatch, *, parameters=None, blocks_absent: bool = False, poscar: str = _POSCAR
):
    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    executable.write_text(_FAKE_VASP.format(blocks_absent=blocks_absent), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    stage = {}
    for name, content in (("POSCAR", poscar), ("POTCAR", _POTCAR)):
        stage[name] = staged / name
        stage[name].write_text(content, encoding="utf-8")
    job = new_job(
        workspace, PACKAGE, files=stage, parameters=parameters or {}, name="packaged vasp-dielectric", install=True
    )
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=300.0)
    return workspace, job.job_id


def test_the_dielectric_run_writes_the_response_and_collects_every_role(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id = _campaign(tmp_path / "ok", monkeypatch, parameters={"publish_data": True})
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id)
    workdir = payload / "run"

    (run,) = _runs(workdir)
    expected = {
        "IBRION": "8",
        "LEPSILON": ".TRUE.",
        "LPEAD": ".TRUE.",
        "NSW": "1",
        "ISMEAR": "-5",
        "LREAL": ".FALSE.",
        "EDIFF": "1E-7",
    }
    assert {tag: run["tags"][tag] for tag in expected} == expected
    assert run["KPOINTS"].splitlines()[2] == "Gamma"

    document = json.loads((workdir / "dielectric.json").read_text(encoding="utf-8"))
    assert document["electronic"][0] == [10.0, 0.5, 0.0]
    assert document["ionic"][2] == [0.0, 0.0, 2.0]
    assert document["static"][2][2] == 14.0
    assert len(document["born_charges"]) == 2 and document["born_charges"][0][0] == [2.0, 0.3, 0.0]
    assert document["piezoelectric_electronic"][0] == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    assert document["piezoelectric_column_order"] == ["XX", "YY", "ZZ", "XY", "YZ", "ZX"]
    assert document["lattice"] == [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    assert document["positions"] == [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    assert (document["species"], document["counts"], document["energy"]) == (["Si"], [2], pytest.approx(-10.8))
    assert (payload / "data" / "vasp" / "dielectric.json").is_file() and (
        payload / "data" / "vasp" / "POSCAR"
    ).is_file()

    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {
        "dielectric_cell",
        "total_energy",
        "dielectric_response",
        "static_relative_permittivity",
        "high_frequency_relative_permittivity",
        "born_effective_charges",
        "piezoelectric_tensor",
    }
    assert not item.unfulfilled
    assert item.outputs["high_frequency_relative_permittivity"].value == pytest.approx(11.0)
    assert item.outputs["static_relative_permittivity"].value == pytest.approx(13.0)
    assert item.outputs["total_energy"].value == pytest.approx(-10.8)
    from httk.analyse.property_records import BornEffectiveChargesRecord, PiezoelectricTensorRecord

    assert isinstance(item.outputs["born_effective_charges"], BornEffectiveChargesRecord)
    assert isinstance(item.outputs["piezoelectric_tensor"], PiezoelectricTensorRecord)
    assert len(item.outputs["born_effective_charges"].value) == 2
    # Voigt order xx yy zz yz xz xy: VASP's XY (0.4 + 1) moves from column 3 to 5, ZX (0.6 + 1) to 4.
    assert item.outputs["piezoelectric_tensor"].value[0] == pytest.approx([1.1, 1.2, 1.3, 1.5, 1.6, 1.4])
    assert item.outputs["dielectric_response"].media_type == "application/json"


@pytest.mark.parametrize(
    ("parameters", "code"),
    (
        ({"incar_tags": {"npar": 2}}, "vasp.unsupported_incar_tag"),
        ({"incar_tags": {"LEPSILON": False}}, "vasp.unsupported_incar_tag"),
        ({"parallel_tag": "NCORE", "parallel_value": 2}, "vasp.unsupported_parallel_tag"),
        ({"parallel_tag": "NPAR", "parallel_value": 2}, "vasp.unsupported_parallel_tag"),
    ),
)
def test_unsupported_parallelization_and_defining_tags_fail_before_vasp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict, code: str
) -> None:
    workspace, job_id = _campaign(tmp_path / "refused", monkeypatch, parameters=parameters)
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == code
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()


def test_an_outcar_without_the_response_blocks_fails(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "missing", monkeypatch, blocks_absent=True)
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.dielectric_missing"
    assert not (payload / "run" / "dielectric.json").exists()


def test_a_cartesian_poscar_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "cartesian", monkeypatch, poscar=_POSCAR.replace("Direct", "Cartesian"))
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_unsupported"
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()
