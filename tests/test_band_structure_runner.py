"""The ``vasp-band-structure`` package, driven end to end through the real manager.

Like ``test_vasp_runners.py``, every job runs in a real workspace under a real
:class:`httk.workflow.TaskManager`, with a fake ``vasp``. This fake is the one of
``test_vasp_runners.py`` cut down to what a band structure needs and extended:
it logs the KPOINTS and POSCAR of every run besides the INCAR, writes a CHGCAR
when ``LCHARG`` is true, and an ``EIGENVAL`` along the line-mode path when
``ICHARG = 11``. The input is a full-precision conventional diamond Si cell with
a little noise, so the symmetry recognition the runner does is exercised.
"""

import json
import runpy
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

pytest.importorskip("httk.atomistic")
pytest.importorskip("spglib")

import httk.core
from httk.atomistic import band_path
from httk.workflow import TaskManager, Workspace
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POTCAR, _failure, _files, _job_state, _payload_of, _runs

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-band-structure"

# A relaxed conventional diamond Si cell, written at full precision with ~1e-8 noise:
# the precision-derived tolerance would find it triclinic (aP2), the default 1e-3 Å finds cF2.
_SI = """Si diamond conventional, relaxed
1.0
   5.4687280000000001   0.0000000000000000   0.0000000000000000
   0.0000000000000000   5.4687280000000001   0.0000000000000000
   0.0000000000000000   0.0000000000000000   5.4687280000000001
Si
8
Direct
  0.0000000012345678  0.9999999987654321  0.0000000043210000
  0.0000000000000000  0.5000000012345678  0.4999999987654322
  0.5000000034567891  0.0000000000000000  0.5000000000000000
  0.4999999976543210  0.5000000000000000  0.0000000011111111
  0.2500000012345678  0.2499999987654321  0.2500000000000000
  0.2500000000000000  0.7500000034567890  0.7500000000000000
  0.7499999987654321  0.2500000000000000  0.7500000012345678
  0.7500000000000000  0.7500000000000000  0.2499999965432100
"""

_FAKE_VASP = '''#!/usr/bin/env python3
"""A fake VASP that writes plausible scf and band outputs, scripted by the test that installed it."""

import json
from pathlib import Path

from httk.codes.vasp import read_incar

attempts = Path("fake-vasp-attempts")
count = int(attempts.read_text()) if attempts.is_file() else 0
attempts.write_text(str(count + 1))
tags = read_incar("INCAR")
chgcar = Path("CHGCAR").is_file() and Path("CHGCAR").stat().st_size > 0
record = dict(tags=tags, KPOINTS=Path("KPOINTS").read_text(), POSCAR=Path("POSCAR").read_text(), CHGCAR=chgcar)
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(record) + "\\n")

if {fail_when}:
    print({failure!r})
    Path("OUTCAR").write_text(" fake vasp 6.4.1\\n   NELM   =     60\\n   NSW    =      0\\n")
    Path("OSZICAR").write_text("DAV:   1    -0.100000000000E+02\\n")
    raise SystemExit(1)

Path("OUTCAR").write_text(
    " fake vasp 6.4.1\\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\\n"
    "   NSW    =      0    number of steps for IOM\\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\\n"
    "   free  energy   TOTEN  =       -10.80000000 eV\\n"
    "   energy  without entropy=      -10.80000000  energy(sigma->0) =      -10.80000000\\n"
    " General timing and accounting informations for this job:\\n"
)
Path("OSZICAR").write_text(
    "DAV:   1    -0.100000000000E+02   -0.10000E+02   -0.30000E+01   128   0.500E+01\\n"
    "DAV:   2    -0.108000000000E+02   -0.50000E+00   -0.10000E+00   128   0.100E+00\\n"
    "   1 F= -.10800000E+02 E0= -.10800000E+02  d E =-.108000E+02\\n"
)
Path("CONTCAR").write_text(Path("POSCAR").read_text())
Path("vasprun.xml").write_text("<modeling></modeling>\\n")
if tags.get("ICHARG") == "11":
    # EIGENVAL along the line-mode path: the segment endpoints, interpolated.
    lines = Path("KPOINTS").read_text().splitlines()
    divisions = int(lines[1])
    ends = [[float(value) for value in line.split("!")[0].split()] for line in lines[4:] if line.strip()]
    points = []
    for start, end in zip(ends[0::2], ends[1::2]):
        for step in range(divisions):
            points.append([a + (b - a) * step / (divisions - 1) for a, b in zip(start, end)])
    eigenval = ["    2    2    1    1", "  0.4E+02  0.3E-09  0.3E-09  0.3E-09  0.5E-15", "  1.0E-004", "  CAR", " Si"]
    eigenval.append("      8  %d  4" % len(points))
    for k in points:
        eigenval += ["", "  %.7E  %.7E  %.7E  %.7E" % (*k, 0.0)]
        eigenval += ["    %d  %.6f  %.6f  1.0  1.0" % (band, -5.0 + 3 * band, -4.9 + 3 * band) for band in range(1, 5)]
    Path("EIGENVAL").write_text("\\n".join(eigenval) + "\\n")
else:
    Path("EIGENVAL").write_text("scf eigenvalues on the uniform mesh\\n")
if tags.get("LCHARG", "").strip(".").upper().startswith("T"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    parameters: dict[str, object] | None = None,
    extra_files: dict[str, str] | None = None,
    fail_when: str = "False",
    failure: str = "",
    poscar: str = _SI,
) -> tuple[Workspace, str]:
    """Submit and run one band-structure job, and return where it landed."""

    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    executable.write_text(_FAKE_VASP.format(fail_when=fail_when, failure=failure), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    stage: dict[str, Path] = {}
    for name, content in {"POSCAR": poscar, "POTCAR": _POTCAR, **(extra_files or {})}.items():
        stage[name] = staged / name
        stage[name].write_text(content, encoding="utf-8")
    job = new_job(
        workspace, PACKAGE, files=stage, parameters=parameters or {}, name="packaged vasp-band-structure", install=True
    )
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=300.0)
    return workspace, job.job_id


def _no_vasp_ran(payload: Path) -> bool:
    return not (payload / "run" / "fake-vasp-runs.jsonl").exists()


@pytest.mark.parametrize("publish_data", (False, True), ids=("default-workdir", "publish-data-opt-in"))
def test_the_band_structure_runs_scf_then_bands_on_the_hpkot_cell(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, publish_data: bool
) -> None:
    from httk.codes.vasp import read_poscar_header
    from httk.core import FileRecord
    from httk.core.digests import sha256_file
    from httk.workflow import collect

    workspace, job_id = _campaign(tmp_path / "bands", monkeypatch, parameters={"publish_data": publish_data})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"

    # The fixture needs the tolerance: with the precision derived from its digits it is triclinic.
    (tmp_path / "noisy.POSCAR").write_text(_SI, encoding="utf-8")
    assert band_path(httk.core.load(str(tmp_path / "noisy.POSCAR"))).bravais_lattice_extended == "aP2"
    # The noisy conventional cell was recognized with the default 1e-3 Å tolerance.
    document = json.loads((workdir / "band_path.json").read_text(encoding="utf-8"))
    assert document["convention"] == "HPKOT"
    assert (document["bravais_lattice"], document["bravais_lattice_extended"], document["spacegroup"]) == (
        "cF",
        "cF2",
        227,
    )
    assert document["symmetry_tolerance"] == 1e-3
    assert document["points"]["GAMMA"] == ["0", "0", "0"]
    assert document["points"]["X"] == ["1/2", "0", "1/2"]
    assert document["path"][0] == ["GAMMA", "X"]
    assert all(len(row) == 3 and all(isinstance(value, int) for value in row) for row in document["transform"])
    assert sum(read_poscar_header(workdir / "POSCAR").counts) == 2

    scf, bands = _runs(workdir)
    # The scf stage ran on an automatic mesh and wrote the density.
    assert scf["KPOINTS"].splitlines()[1:3] == ["0", "Monkhorst-Pack"]
    assert (scf["tags"]["LCHARG"], scf["tags"]["ICHARG"], scf["tags"]["ISMEAR"]) == (".TRUE.", "2", "-5")
    # The band stage read it along the line-mode path, on the same cell.
    bands_kpoints = (workdir / "KPOINTS.bands").read_text(encoding="utf-8")
    assert bands["KPOINTS"] == bands_kpoints
    assert bands_kpoints.splitlines()[:4] == ["HPKOT cF2 band path (httk)", "40", "Line-mode", "Reciprocal"]
    assert (bands["tags"]["ICHARG"], bands["tags"]["LCHARG"], bands["CHGCAR"]) == ("11", ".FALSE.", True)
    assert scf["POSCAR"] == bands["POSCAR"]
    assert all(scf["tags"].get(tag) == bands["tags"].get(tag) for tag in ("PREC", "ENCUT", "ENAUG", "MAGMOM", "NBANDS"))

    # scf/ archived the self-consistent run; the CHGCAR stayed at the top.
    for name in ("INCAR", "KPOINTS", "POSCAR", "OUTCAR", "OSZICAR", "vasprun.xml", "vasp-run-report.json"):
        assert (workdir / "scf" / name).is_file(), name
    assert (workdir / "CHGCAR").is_file() and not (workdir / "scf" / "CHGCAR").exists()
    state = _job_state(payload)
    assert (state["bravais_lattice_extended"], state["spacegroup"]) == ("cF2", 227)
    assert state["scf_energy"] == pytest.approx(-10.8)
    assert state["scf_incar"]["LCHARG"] == ".TRUE."
    assert (state["scf_classification"], state["bands_classification"]) == ("completed", "completed")

    if publish_data:
        published = _files(payload / "data")
        for name in ("scf/OUTCAR", "bands/POSCAR", "bands/EIGENVAL", "bands/band_path.json", "bands/KPOINTS"):
            assert name in published, name
    else:
        assert not (payload / "data").exists()

    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {"band_structure_cell", "total_energy", "band_eigenvalues", "band_path"}
    assert not item.unfulfilled
    assert item.outputs["total_energy"].value == pytest.approx(-10.8)
    assert len(item.outputs["band_structure_cell"].sites.reduced_coords) == 2
    where = payload / "data" / "bands" if publish_data else workdir
    for role, name in (("band_eigenvalues", "EIGENVAL"), ("band_path", "band_path.json")):
        record = item.outputs[role]
        assert isinstance(record, FileRecord)
        # collect prefixes a workspace-relative URL with data/<workspace name>/ (unregistered: its directory's name).
        relative = (where / name).resolve().relative_to(workspace.root.resolve()).as_posix()
        assert record.url == f"data/workspace/{relative}"
        assert record.sha256 == sha256_file(where / name)
    assert item.outputs["band_path"].media_type == "application/json"


def test_a_replayed_scf_icharg_remedy_leaves_the_band_stage_non_self_consistent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id = _campaign(
        tmp_path / "fexcf",
        monkeypatch,
        fail_when="count == 0",
        failure="FEXCF: supplied exchange-correlation table is too small",
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    _first, scf, bands = _runs(payload / "run")
    # The fexcf remedy set ICHARG = 2 with its mixing in scf; the band stage replays the
    # mixing but runs with ICHARG = 11.
    assert (scf["tags"]["ICHARG"], scf["tags"]["AMIX"]) == ("2", "0.1")
    assert (bands["tags"]["ICHARG"], bands["tags"]["AMIX"]) == ("11", "0.1")
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"fexcf": 1}


def test_a_geometry_remedy_is_never_applied(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "zpotrf", monkeypatch, fail_when="count == 0", failure="LAPACK: Routine ZPOTRF failed!"
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    first, scf, bands = _runs(payload / "run")
    # The zpotrf ladder's first rung scales the lattice; it is skipped, and the next
    # rung, bump_kpoints, raised the scf mesh by one instead.
    assert first["POSCAR"] == scf["POSCAR"] == bands["POSCAR"]
    assert [int(n) + 1 for n in first["KPOINTS"].splitlines()[3].split()] == [
        int(n) for n in scf["KPOINTS"].splitlines()[3].split()
    ]
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"zpotrf": 2}
    assert [change["operation"] for change in history["events"][0]["changes"]] == ["bump_kpoints"]
    # The band stage still runs on its line-mode path.
    assert bands["KPOINTS"] == (payload / "run" / "KPOINTS.bands").read_text(encoding="utf-8")


def test_a_band_override_with_another_density_grid_fails_before_the_band_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    template = (PACKAGE / "templates" / "INCAR.bands").read_text(encoding="utf-8")
    workspace, job_id = _campaign(
        tmp_path / "incompatible", monkeypatch, extra_files={"INCAR.bands": template.replace("ENCUT=600", "ENCUT=520")}
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.charge_density_incompatible"
    assert failure["details"]["differences"] == {"ENCUT": {"scf": "600", "bands": "520"}}
    assert len(_runs(payload / "run")) == 1


@pytest.mark.parametrize(
    ("source", "parameters", "extra_files", "poscar"),
    (
        ("incar_tags", {"incar_tags": {"magmom": "2*1"}}, None, _SI),
        ("the comment of POSCAR", None, None, _SI.replace("relaxed", "relaxed [MAGMOM=8*1]", 1)),
        ("files/INCAR.scf", None, {"INCAR.scf": "NELECT = 7\n"}, _SI),
        ("files/INCAR.bands", None, {"INCAR.bands": "ICHARG = 11\nNGX = 24\n"}, _SI),
    ),
    ids=("incar-tags", "poscar-comment", "payload-incar", "payload-fft-grid"),
)
def test_a_cell_dependent_override_is_refused_before_any_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source: str,
    parameters: dict[str, object] | None,
    extra_files: dict[str, str] | None,
    poscar: str,
) -> None:
    workspace, job_id = _campaign(
        tmp_path / "refused", monkeypatch, parameters=parameters, extra_files=extra_files, poscar=poscar
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.cell_dependent_tag_unsupported"
    assert failure["details"]["source"] == source
    assert "HPKOT primitive cell" in failure["message"]
    assert _no_vasp_ran(payload)


def test_a_hybrid_functional_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "hybrid", monkeypatch, parameters={"incar_tags": {"LHFCALC": ".TRUE."}})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.unsupported_calculation"
    assert failure["details"] == {"source": "incar_tags", "tags": {"LHFCALC": ".TRUE."}}
    assert _no_vasp_ran(payload)


def test_an_icharg_override_fails_by_name_before_any_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # incar_tags apply to both stages, so any ICHARG there is wrong for one of them.
    workspace, job_id = _campaign(tmp_path / "icharg", monkeypatch, parameters={"incar_tags": {"ICHARG": 2}})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.invalid_band_inputs"
    assert failure["details"] == {"source": "incar_tags", "ICHARG": "2"}
    assert _no_vasp_ran(payload)


def test_an_scf_that_writes_no_charge_density_fails_the_band_stage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    template = (PACKAGE / "templates" / "INCAR.scf").read_text(encoding="utf-8")
    workspace, job_id = _campaign(
        tmp_path / "nochgcar",
        monkeypatch,
        extra_files={"INCAR.scf": template.replace("LCHARG=.TRUE.", "LCHARG=.FALSE.")},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.charge_density_missing"
    assert len(_runs(payload / "run")) == 1


def test_the_density_comparison_normalizes_values(tmp_path: Path) -> None:
    band_inputs_problem = runpy.run_path(str(PACKAGE / "run.py"))["band_inputs_problem"]
    (tmp_path / "CHGCAR").write_text("density\n", encoding="utf-8")
    scf = {"ICHARG": "2", "LASPH": ".TRUE.", "ENCUT": " 600 "}

    def problem(incar: str) -> Any:
        (tmp_path / "INCAR").write_text(incar, encoding="utf-8")
        return band_inputs_problem(SimpleNamespace(workdir=tmp_path, state={"scf_incar": scf}))

    # Case and surrounding whitespace do not matter, and ICHARG is not a density tag.
    assert problem("ICHARG = 11\nLASPH = .true.\nENCUT =   600\n") is None
    # An absent tag equals only an absent one.
    code, _message, details = problem("ICHARG = 11\nLASPH = .TRUE.\nENCUT = 600\nGGA = PE\n")
    assert (code, details) == (
        "vasp.charge_density_incompatible",
        {"differences": {"GGA": {"scf": None, "bands": "PE"}}},
    )
    code, _message, details = problem("ICHARG = 11\nENCUT = 600\n")
    assert details == {"differences": {"LASPH": {"scf": ".TRUE.", "bands": None}}}


def test_a_rerun_scf_step_repeats_only_the_scf_work(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Rerunning run_scf (called directly, as test_vasp_runners.py calls a stage step) on
    # restored scf inputs repeats the scf work, re-archiving scf/, and advances to
    # stage_bands, never straight to publish. No crash is injected: this checks that the
    # step's outcome is fixed by the step, not by job state a later step wrote.
    workspace, job_id = _campaign(tmp_path / "rerun", monkeypatch)
    _, payload = _payload_of(workspace, job_id)
    workdir = payload / "run"
    for name in ("INCAR", "KPOINTS"):
        shutil.copyfile(workdir / "scf" / name, workdir / name)
    (workdir / "scf" / "OUTCAR").unlink()

    run_scf = runpy.run_path(str(PACKAGE / "run.py"))["run_scf"]
    advanced: list[tuple[str, dict[str, Any]]] = []
    attempt = SimpleNamespace(
        state=_job_state(payload),
        payload=payload,
        workdir=workdir,
        parameter=lambda name, default=None: default,
        setting=lambda name, default=None: str(tmp_path / "rerun" / "fake-vasp") if name == "vasp.command" else default,
        log=SimpleNamespace(append=lambda *arguments: None),
        advance=lambda step, state=None: advanced.append((step, dict(state or {}))),
        fail=lambda *arguments, **keywords: pytest.fail(f"run_scf failed: {arguments}"),
    )
    run_scf(attempt)

    assert [step for step, _ in advanced] == ["stage_bands"]
    assert advanced[0][1]["scf_incar"]["ICHARG"] == "2"
    runs = _runs(workdir)
    assert len(runs) == 3
    assert runs[-1]["tags"]["ICHARG"] == "2" and runs[-1]["KPOINTS"] == runs[0]["KPOINTS"]
    assert (workdir / "scf" / "OUTCAR").is_file()


def _prepare(tmp_path: Path, parameters: dict[str, object]) -> tuple[list[str], list[tuple[Any, ...]]]:
    """Run the prepare step directly on the payload *tmp_path*; return its advances and failures."""

    prepare = runpy.run_path(str(PACKAGE / "run.py"))["prepare"]
    workdir = tmp_path / "run"
    workdir.mkdir()
    advanced: list[str] = []
    failed: list[tuple[Any, ...]] = []
    prepare(
        SimpleNamespace(
            payload=tmp_path,
            workdir=workdir,
            parameter=lambda name, default=None: parameters.get(name, default),
            log=SimpleNamespace(append=lambda *arguments: None),
            advance=lambda step, state=None: advanced.append(step),
            fail=lambda code, message, details=None: failed.append((code, details)),
        )
    )
    return advanced, failed


@pytest.mark.parametrize("species", (("Si", "C"), ("C", "Si")))
def test_the_band_cell_keeps_the_input_species_order(tmp_path: Path, species: tuple[str, str]) -> None:
    from httk.codes.vasp import read_poscar_header

    # Not named POSCAR: the structure is read by format, whatever its file name.
    (tmp_path / "files").mkdir()
    (tmp_path / "files" / "ZnS.relaxed").write_text(
        "zincblende\n1.0\n4.36 0 0\n0 4.36 0\n0 0 4.36\n"
        f"{species[0]} {species[1]}\n4 4\nDirect\n"
        "0 0 0\n0 0.5 0.5\n0.5 0 0.5\n0.5 0.5 0\n"
        "0.25 0.25 0.25\n0.25 0.75 0.75\n0.75 0.25 0.75\n0.75 0.75 0.25\n",
        encoding="utf-8",
    )
    advanced, failed = _prepare(tmp_path, {"poscar": "files/ZnS.relaxed"})

    assert (advanced, failed) == (["stage_scf"], [])
    header = read_poscar_header(tmp_path / "run" / "POSCAR")
    assert (header.species, header.counts) == (species, (1, 1))
    assert json.loads((tmp_path / "run" / "band_path.json").read_text(encoding="utf-8"))["spacegroup"] == 216


@pytest.mark.parametrize("divisions", (40.7, True, 1, "40"))
def test_invalid_band_divisions_fail_by_name(tmp_path: Path, divisions: object) -> None:
    (tmp_path / "files").mkdir()
    (tmp_path / "files" / "POSCAR").write_text(_SI, encoding="utf-8")
    advanced, failed = _prepare(tmp_path, {"band_divisions": divisions})

    assert (advanced, failed) == ([], [("vasp.invalid_band_inputs", {"band_divisions": divisions})])
    assert not (tmp_path / "run" / "KPOINTS.bands").exists()


def test_the_band_template_keeps_the_scf_density_settings() -> None:
    from httk.codes.vasp import read_incar

    scf = read_incar(PACKAGE / "templates" / "INCAR.scf")
    bands = read_incar(PACKAGE / "templates" / "INCAR.bands")
    density = runpy.run_path(str(PACKAGE / "run.py"))["DENSITY_TAGS"]
    assert {tag: scf.get(tag) for tag in density} == {tag: bands.get(tag) for tag in density}
    assert (scf["ICHARG"], scf["LCHARG"], bands["ICHARG"], bands["LCHARG"]) == ("2", ".TRUE.", "11", ".FALSE.")
    assert "LORBIT" not in bands


@pytest.mark.parametrize(
    ("name", "incar", "expected"),
    (
        ("INCAR.bands", "ICHARG = 2\n", ("vasp.invalid_band_inputs", {"source": "files/INCAR.bands", "ICHARG": "2"})),
        (
            "INCAR.bands",
            "ICHARG = 11\nMETAGGA = SCAN\n",
            ("vasp.unsupported_calculation", {"source": "files/INCAR.bands", "tags": {"METAGGA": "SCAN"}}),
        ),
        (
            "INCAR.scf",
            "LSORBIT = T\n",
            ("vasp.unsupported_calculation", {"source": "files/INCAR.scf", "tags": {"LSORBIT": "T"}}),
        ),
    ),
    ids=("bands-icharg", "bands-metagga", "scf-lsorbit"),
)
def test_a_payload_incar_is_checked_in_prepare(
    tmp_path: Path, name: str, incar: str, expected: tuple[str, dict[str, object]]
) -> None:
    (tmp_path / "files").mkdir()
    (tmp_path / "files" / "POSCAR").write_text(_SI, encoding="utf-8")
    (tmp_path / "files" / name).write_text(incar, encoding="utf-8")
    advanced, failed = _prepare(tmp_path, {})

    assert (advanced, failed) == ([], [expected])
