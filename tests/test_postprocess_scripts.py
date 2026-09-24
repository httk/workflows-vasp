"""Postprocess scripts of the packaged ``vasp-relax`` workflow.

Ported from httk-workflow's ``tests/test_vasp_collect.py`` (``git show
15b0b36:tests/test_vasp_collect.py``), keeping only its last three tests: the
ones that run this package's ``scripts/relaxation_report`` and
``scripts/relaxation_plot`` against a collected job record. The earlier tests
of that file exercised the ``collect_vasp_*`` collector functions themselves,
which stayed in httk-workflow's ``httk.workflow.vasp.collect`` module (still
covered by that module's own tests there); this file is only the packaged,
directory-sourced side.

``registered_workflow("vasp-relax")`` is gone along with the built-in
provider, so the workflow here is resolved straight from the checked-out
directory with :func:`httk.workflow.packages.load_workflow_package`.
"""

import json
import re
import stat
from pathlib import Path, PurePosixPath

import pytest

pytest.importorskip("httk.atomistic")

from httk.workflow.collecting import JobRecord
from httk.workflow.packages import load_workflow_package
from httk.workflow.postprocessing import run_postprocess_script

REPO_ROOT = Path(__file__).resolve().parent.parent
_VASP_RELAX = REPO_ROOT / "vasp-relax"

_OUTCAR = """ vasp.5.2.12 synthetic
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -26.00000000 eV
   energy  without entropy=      -26.00000000  energy(sigma->0) =      -26.00000000
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -27.00000000 eV
   energy  without entropy=      -27.00000000  energy(sigma->0) =      -27.00000000
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -27.09328752 eV
   energy  without entropy=      -27.09328752  energy(sigma->0) =      -27.09328752
  General timing and accounting informations for this job:
"""
_CONTCAR = """silicon
1.0
2.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 2.0
Si
2
Direct
0.0 0.0 0.0
0.5 0.5 0.5
"""


def _record(root: Path) -> JobRecord:
    return JobRecord(
        workspace_root=root,
        workspace_id="ws",
        job_id="12345678-1234-4234-8234-123456789abc",
        job_key="job--12345678-1234-4234-8234-123456789abc",
        job={"workflow": "vasp.relax"},
        runner_provenance=None,
        state="succeeded",
        failure=None,
        placement=PurePosixPath("jobs"),
        payload_path=PurePosixPath("jobs/job--12345678-1234-4234-8234-123456789abc"),
        workdir_path=PurePosixPath("run"),
        data_path=PurePosixPath("data"),
        data_generation=1,
        provenance={},
        runner_steps=None,
        children={},
        declarations={},
    )


def _write(root: Path, name: str) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / name).write_text(_CONTCAR if name == "CONTCAR" else _OUTCAR, encoding="utf-8")


def test_packaged_relaxation_report_runs_from_published_data(tmp_path: Path) -> None:
    _write(tmp_path / "data" / "vasp", "CONTCAR")
    _write(tmp_path / "data" / "vasp", "OUTCAR")
    (tmp_path / "run").mkdir()
    workflow = load_workflow_package(_VASP_RELAX, register=False)
    script = _VASP_RELAX / "scripts" / "relaxation_report"
    assert script.stat().st_mode & stat.S_IXUSR

    result = run_postprocess_script(workflow, "relaxation-report", _record(tmp_path))
    assert result.returncode == 0
    report = json.loads((result.output_dir / "relaxation_report.json").read_text(encoding="utf-8"))
    assert report["final_energy"] == pytest.approx(-27.09328752)
    assert report["structure_files"] == ["vasp/CONTCAR"]
    assert "vasp/CONTCAR" in (result.output_dir / "relaxation_report.txt").read_text(encoding="utf-8")


def test_packaged_relaxation_plot_runs_from_published_data(tmp_path: Path) -> None:
    _write(tmp_path / "data" / "vasp", "OUTCAR")
    (tmp_path / "run").mkdir()
    workflow = load_workflow_package(_VASP_RELAX, register=False)

    result = run_postprocess_script(workflow, "relaxation-plot", _record(tmp_path))
    assert result.returncode == 0
    svg = (result.output_dir / "relaxation_energies.svg").read_text(encoding="utf-8")
    points = re.search(r'<polyline points="([^"]+)"', svg)
    assert points is not None
    assert len(points.group(1).split()) == 3
    assert "relaxation_energies.svg" in result.stdout


def test_packaged_relaxation_plot_tolerates_missing_outcar(tmp_path: Path) -> None:
    (tmp_path / "run").mkdir()
    workflow = load_workflow_package(_VASP_RELAX, register=False)

    result = run_postprocess_script(workflow, "relaxation-plot", _record(tmp_path))
    assert result.returncode == 0
    assert "no OUTCAR" in result.stdout
    assert not (result.output_dir / "relaxation_energies.svg").exists()
