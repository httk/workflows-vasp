"""Each package's ``collect.py`` reads its outputs from where its runner leaves them.

Ported from httk-workflow-vasp's ``tests/test_vasp_collect.py``, which tested
the per-workflow collectors while they lived in ``httk.codes.vasp.collect``.
The layouts are now decided by the hooks here, so the layout matrix is tested
here; the helpers the hooks build on are tested in httk-workflow-vasp.
"""

import runpy
from pathlib import Path, PurePosixPath

import pytest

pytest.importorskip("httk.atomistic")

import httk.core
from httk.workflow.collecting import JobRecord

REPO_ROOT = Path(__file__).resolve().parent.parent

_POSCAR = """silicon
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
_OUTCAR = """ vasp.5.2.12 synthetic
   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
   free  energy   TOTEN  =       -27.09328752 eV
   energy  without entropy=      -27.09328752  energy(sigma->0) =      -27.09328752
  General timing and accounting informations for this job:
"""


def _record(root: Path, parameters: dict[str, str], *, transactional: bool) -> JobRecord:
    return JobRecord(
        workspace_root=root,
        workspace_id="ws",
        job_id="12345678-1234-4234-8234-123456789abc",
        job_key="job--12345678-1234-4234-8234-123456789abc",
        job={"parameters": parameters},
        runner_provenance=None,
        state="succeeded",
        failure=None,
        placement=PurePosixPath("jobs"),
        payload_path=PurePosixPath("jobs/job--12345678-1234-4234-8234-123456789abc"),
        workdir_path=PurePosixPath("run"),
        data_path=PurePosixPath("data") if transactional else None,
        data_generation=1 if transactional else None,
        provenance={},
        runner_steps=None,
        children={},
        declarations={},
    )


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


@pytest.mark.parametrize("transactional", (False, True))
@pytest.mark.parametrize("prefix", (None, "", "custom/results"))
@pytest.mark.parametrize("package", ("vasp-relax", "vasp-relax-bash", "vasp-static", "vasp-relax-static"))
def test_collect_hooks_follow_their_runner_layouts(
    tmp_path: Path, transactional: bool, prefix: str | None, package: str
) -> None:
    relax_static = package == "vasp-relax-static"
    parameters = {} if prefix is None else {"data_prefix": prefix}
    if transactional:
        # The runners' own data_prefix defaults: vasp-relax-static publishes its stages at the top.
        default_prefix = "" if relax_static else "vasp"
        root = tmp_path / "data" / (default_prefix if prefix is None else prefix)
        structure = root / "relax" / "CONTCAR" if relax_static else root / "CONTCAR"
        energy = root / "static" / "OUTCAR" if relax_static else root / "OUTCAR"
    else:
        root = tmp_path / "run"
        structure = root / "relax" / "CONTCAR" if relax_static else root / "CONTCAR"
        energy = root / "OUTCAR"
    _write(structure, _POSCAR)
    _write(energy, _OUTCAR)
    if relax_static:
        # The archived relaxation energy must not become the static result.
        _write(root / "relax" / "OUTCAR", _OUTCAR.replace("-27.09328752", "-20.0"))

    collect = runpy.run_path(str(REPO_ROOT / package / "collect.py"))["collect"]
    outputs = collect(_record(tmp_path, parameters, transactional=transactional))

    assert set(outputs) == ({"total_energy"} if package == "vasp-static" else {"relaxed_structure", "total_energy"})
    assert isinstance(outputs["total_energy"], httk.core.DataRecord)
    assert outputs["total_energy"].value == pytest.approx(-27.09328752)
