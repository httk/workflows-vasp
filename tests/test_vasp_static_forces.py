"""``vasp.static`` collects ``cartesian_site_forces`` only when the run kept the staged geometry."""

import json
import logging
import runpy
from pathlib import Path, PurePosixPath

import pytest

pytest.importorskip("httk.atomistic")

from httk.codes.vasp import job_remedy_history_path
from httk.core.property_records import CartesianSiteForcesRecord
from httk.workflow.collecting import JobRecord
from httk.workflow.packages import load_workflow_package

REPO_ROOT = Path(__file__).resolve().parent.parent
PAYLOAD = PurePosixPath("jobs/job--12345678-1234-4234-8234-123456789abc")

# Species are not contiguous on purpose: the forces must keep the POSCAR's site order.
_POSCAR = """rocksalt fragment
1.0
4.0 0.0 0.0
0.0 4.0 0.0
0.0 0.0 4.0
Na Cl Na
1 1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
0.98 0.5 0.0
"""
# The OUTCAR writes the third site at its periodic image x = -0.08 Å.
_POSITIONS = ((0.0, 0.0, 0.0), (2.0, 2.0, 2.0), (-0.08, 2.0, 0.0))
_FORCES = ((0.1, -0.2, 0.3), (-0.4, 0.5, -0.6), (0.3, -0.3, 0.3))


def _outcar(*frames: tuple[tuple[tuple[float, ...], ...], tuple[tuple[float, ...], ...]], edge: float = 4.0) -> str:
    ions = "  ".join("1" for _ in frames[0][0])
    text = f" vasp.6.4.1 synthetic\n   ions per type =  {ions}\n   NSW    =     0    number of steps for IOM\n"
    # The cubic cell of edge *edge*, as VASP prints it beside the reciprocal vectors.
    lattice = "".join(
        "     " + "  ".join(f"{edge if i == j else 0.0:.9f}" for j in range(3)) + "     0.0 0.0 0.0\n" for i in range(3)
    )
    for positions, forces in frames:
        rows = "".join(
            "".join(f" {x:12.5f}" for x in r) + "   " + "".join(f" {x:14.8f}" for x in f) + "\n"
            for r, f in zip(positions, forces, strict=True)
        )
        text += (
            "  FORCE on cell =-STRESS in cart. coord.  units (eV):\n"
            "  in kB       0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            "       direct lattice vectors                 reciprocal lattice vectors\n"
            + lattice
            + " POSITION                                       TOTAL-FORCE (eV/Angst)\n"
            " " + "-" * 83 + "\n" + rows + " " + "-" * 83 + "\n"
            "    total drift:                                0.000000      0.000000      0.000000\n"
            "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\n"
            "  free  energy   TOTEN  =       -10.80000000 eV\n"
            "  energy  without entropy=      -10.80000000  energy(sigma->0) =      -10.80000000\n"
        )
    return text + " General timing and accounting informations for this job:\n"


def _collect(
    tmp_path: Path, outcar: str, parameters: dict[str, object] | None = None, poscar: str = _POSCAR
) -> dict[str, object]:
    (tmp_path / "run").mkdir()
    (tmp_path / "run" / "OUTCAR").write_text(outcar, encoding="utf-8")
    (tmp_path / PAYLOAD / "files").mkdir(parents=True, exist_ok=True)
    (tmp_path / PAYLOAD / "files" / "POSCAR").write_text(poscar, encoding="utf-8")
    record = JobRecord(
        workspace_root=tmp_path,
        workspace_id="ws",
        job_id="12345678-1234-4234-8234-123456789abc",
        job_key=PAYLOAD.name,
        job={"parameters": parameters or {}},
        runner_provenance=None,
        state="succeeded",
        failure=None,
        placement=PurePosixPath("jobs"),
        payload_path=PAYLOAD,
        workdir_path=PurePosixPath("run"),
        data_path=None,
        provenance={},
        runner_steps=None,
        children={},
        declarations={},
    )
    return runpy.run_path(str(REPO_ROOT / "vasp-static" / "collect.py"))["collect"](record)


def test_forces_of_the_kept_geometry_are_collected_in_poscar_order(tmp_path: Path) -> None:
    outputs = _collect(tmp_path, _outcar((_POSITIONS, _FORCES)))
    forces = outputs["cartesian_site_forces"]
    assert isinstance(forces, CartesianSiteForcesRecord)
    assert forces.value == [list(row) for row in _FORCES]
    assert outputs["total_energy"].value == pytest.approx(-10.8)


def test_several_ionic_configurations_leave_the_forces_out(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.INFO)
    outputs = _collect(tmp_path, _outcar((_POSITIONS, _FORCES), (_POSITIONS, _FORCES)))
    assert set(outputs) == {"total_energy"}
    assert "exactly one ionic configuration" in caplog.text


def test_a_rattling_remedy_leaves_the_forces_out(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.INFO)
    history = job_remedy_history_path(tmp_path / PAYLOAD)
    history.parent.mkdir(parents=True)
    history.write_text(json.dumps({"events": [{"changes": [{"operation": "bump_bands", "value": 1.1}]}]}))
    outputs = _collect(tmp_path, _outcar((_POSITIONS, _FORCES)), {"rattle_amplitude": 0.01})
    assert set(outputs) == {"total_energy"}
    assert "a remedy rattled the structure" in caplog.text


def test_moved_positions_leave_the_forces_out(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.INFO)
    moved = (_POSITIONS[0], (2.0, 2.0, 2.001), _POSITIONS[2])
    outputs = _collect(tmp_path, _outcar((moved, _FORCES)))
    assert set(outputs) == {"total_energy"}
    assert "position of site 2 differs" in caplog.text


def test_a_scaled_cell_leaves_the_forces_out(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    # One site at the origin keeps its position under scale_lattice: only the cell check sees the change.
    caplog.set_level(logging.INFO)
    poscar = "one site\n1.0\n4.0 0.0 0.0\n0.0 4.0 0.0\n0.0 0.0 4.0\nNa\n1\nDirect\n0.0 0.0 0.0\n"
    outcar = _outcar((((0.0, 0.0, 0.0),), ((0.0, 0.0, 0.0),)), edge=4.04)
    outputs = _collect(tmp_path, outcar, {"rattle_amplitude": 0}, poscar=poscar)
    assert set(outputs) == {"total_energy"}
    assert "cell differs" in caplog.text


def test_the_declaration_matches_the_manifest() -> None:
    provider = load_workflow_package(REPO_ROOT / "vasp-static", register=False)
    assert provider is not None
