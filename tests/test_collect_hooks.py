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
@pytest.mark.parametrize(
    "package",
    (
        "vasp-relax",
        "vasp-relax-bash",
        "vasp-static",
        "vasp-relax-static",
        "vasp-relax-static-coarse",
        "vasp-relax-static-hse",
        "vasp-general",
    ),
)
def test_collect_hooks_follow_their_runner_layouts(
    tmp_path: Path, transactional: bool, prefix: str | None, package: str
) -> None:
    # vasp-relax-static-coarse and -hse keep vasp-relax-static's layout: relax/ for the structure, the static at the top.
    relax_static = package in ("vasp-relax-static", "vasp-relax-static-coarse", "vasp-relax-static-hse")
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
    if package == "vasp-relax-static-coarse":
        # The driver's report, next to relax/: this job's final static is within the limits.
        _write(root / "relax-driver" / "report.json", '{"within_limits": true}')

    collect = runpy.run_path(str(REPO_ROOT / package / "collect.py"))["collect"]
    outputs = collect(_record(tmp_path, parameters, transactional=transactional))

    structure_role = "final_structure" if package == "vasp-general" else "relaxed_structure"
    assert set(outputs) == ({"total_energy"} if package == "vasp-static" else {structure_role, "total_energy"})
    assert isinstance(outputs["total_energy"], httk.core.DataRecord)
    assert outputs["total_energy"].value == pytest.approx(-27.09328752)


@pytest.mark.parametrize("transactional", (False, True))
def test_the_coarse_driver_harvests_nothing_from_a_final_static_outside_the_limits(
    tmp_path: Path, transactional: bool
) -> None:
    root = tmp_path / ("data" if transactional else "run")
    _write(root / "relax" / "CONTCAR", _POSCAR)
    _write(root / ("static/OUTCAR" if transactional else "OUTCAR"), _OUTCAR)
    _write(root / "relax-driver" / "report.json", '{"within_limits": false}')

    collect = runpy.run_path(str(REPO_ROOT / "vasp-relax-static-coarse" / "collect.py"))["collect"]

    assert collect(_record(tmp_path, {}, transactional=transactional)) == {}


@pytest.mark.parametrize("contcar", ("absent", "empty"))
def test_the_general_collect_hook_leaves_out_a_missing_structure(tmp_path: Path, contcar: str) -> None:
    _write(tmp_path / "run" / "OUTCAR", _OUTCAR)
    if contcar == "empty":
        _write(tmp_path / "run" / "CONTCAR", "")

    collect = runpy.run_path(str(REPO_ROOT / "vasp-general" / "collect.py"))["collect"]
    outputs = collect(_record(tmp_path, {}, transactional=False))

    assert set(outputs) == {"total_energy"}


@pytest.mark.parametrize("transactional", (False, True))
@pytest.mark.parametrize("prefix", (None, "", "custom/results"))
def test_the_band_structure_collect_hook_follows_its_runner_layout(
    tmp_path: Path, transactional: bool, prefix: str | None
) -> None:
    from httk.core.digests import sha256_file

    parameters = {} if prefix is None else {"data_prefix": prefix}
    if transactional:
        # The runner publishes its stages under scf/ and bands/ below the (default empty) prefix.
        root = tmp_path / "data" / (prefix or "")
        bands, scf = root / "bands", root / "scf"
    else:
        bands = tmp_path / "run"
        scf = bands / "scf"
    _write(bands / "POSCAR", _POSCAR)
    _write(bands / "EIGENVAL", "fake eigenvalues\n")
    _write(bands / "band_path.json", '{"convention": "HPKOT"}\n')
    _write(scf / "OUTCAR", _OUTCAR)
    # Only the scf run's energy is the result; a top-level OUTCAR is the band run's.
    _write(bands / "OUTCAR", _OUTCAR.replace("-27.09328752", "-20.0"))

    collect = runpy.run_path(str(REPO_ROOT / "vasp-band-structure" / "collect.py"))["collect"]
    outputs = collect(_record(tmp_path, parameters, transactional=transactional))

    assert set(outputs) == {"band_structure_cell", "total_energy", "band_eigenvalues", "band_path"}
    assert outputs["total_energy"].value == pytest.approx(-27.09328752)
    assert len(outputs["band_structure_cell"].sites.reduced_coords) == 2
    for role, name, media_type in (
        ("band_eigenvalues", "EIGENVAL", "application/octet-stream"),
        ("band_path", "band_path.json", "application/json"),
    ):
        record = outputs[role]
        assert isinstance(record, httk.core.FileRecord)
        assert record.url == (bands / name).relative_to(tmp_path).as_posix()
        assert (record.name, record.media_type, record.size) == (name, media_type, (bands / name).stat().st_size)
        assert record.sha256 == sha256_file(bands / name)
