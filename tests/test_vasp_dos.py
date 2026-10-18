"""The ``vasp-dos`` package, driven end to end through the real manager.

Like ``test_band_structure_runner.py``, every job runs in a real workspace under a
real :class:`httk.workflow.TaskManager`, with a fake ``vasp``. This fake is that
file's, adapted: for an ``ICHARG = 11`` run it writes an EIGENVAL with
occupations on three k-points of the KPOINTS mesh, with a known indirect gap
(valence maximum at the first k-point, conduction minimum at the last) and a
larger direct gap, and a DOSCAR of ``NEDOS`` rows with the Fermi energy 5 eV; and
a CHGCAR when ``LCHARG`` is true. ``SCRIPT`` shifts the spin-down valence and
conduction bands, or half-fills one conduction state (a metal).
"""

import json
import re
import runpy
from pathlib import Path

import pytest

pytest.importorskip("httk.atomistic")
pytest.importorskip("httk.analyse")

from httk.workflow import TaskManager, Workspace
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POSCAR, _POTCAR, _failure, _files, _job_state, _payload_of, _runs

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-dos"

_FAKE_VASP = '''#!/usr/bin/env python3
"""A fake VASP that writes plausible scf and nscf outputs, scripted by the test that installed it."""

import json
from pathlib import Path

from httk.codes.vasp import read_incar

SCRIPT = {script!r}
FERMI = 5.0

tags = read_incar("INCAR")
chgcar = Path("CHGCAR").is_file() and Path("CHGCAR").stat().st_size > 0
record = dict(tags=tags, KPOINTS=Path("KPOINTS").read_text(), POSCAR=Path("POSCAR").read_text(), CHGCAR=chgcar)
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(record) + "\\n")

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
    ispin = int(tags.get("ISPIN", "2"))
    n1 = int(Path("KPOINTS").read_text().splitlines()[3].split()[0])
    kpoints = [(i / n1, 0.0, 0.0) for i in range(3)]
    # Spin up, bands by k-points: VBM 4.5 at k0, CBM 5.7 at k2, smallest direct gap 1.7 at k2.
    up = [[1.0, 1.0, 1.0], [4.5, 4.2, 4.0], [6.5, 6.0, 5.7], [8.0, 8.0, 8.0]]
    occupations = [[1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0]]
    if SCRIPT["metal"]:
        occupations[2][2] = 0.5
    if SCRIPT["tetrahedron_metal"]:
        # Bloechl-corrected tetrahedron occupations at the Fermi level of a metal.
        occupations[1][2], occupations[2][2], occupations[2][1] = 1.0012, 0.5, -0.0017
    valence, conduction = SCRIPT["down_shift"]
    down = [[e + (valence if band < 2 else conduction) for e in row] for band, row in enumerate(up)]
    lines = ["    2    2    1    %d" % ispin, "  0.4E+02  0.3E-09  0.3E-09  0.3E-09  0.5E-15", "  1.0E-004", "  CAR"]
    lines += [" Si", "      8  3  4"]
    for k, point in enumerate(kpoints):
        lines += ["", "  %.7E  %.7E  %.7E  %.7E" % (*point, 1 / 3)]
        for band in range(4):
            energies = [up[band][k], down[band][k]][:ispin]
            values = ["%.6f" % e for e in energies] + ["%.4f" % occupations[band][k]] * ispin
            lines.append("    %d  " % (band + 1) + "  ".join(values))
    Path("EIGENVAL").write_text("\\n".join(lines) + "\\n")
    nedos = int(tags.get("NEDOS", "301"))
    doscar = ["    2    2    1    0", "  0.4E+02  0.3E-09  0.3E-09  0.3E-09  0.5E-15", "  1.0E-004", "  CAR", " Si"]
    doscar.append("  10.00000000  -5.00000000  %d  %.8f  1.00000000" % (nedos, FERMI))
    for row in range(nedos):
        density = ["0.5000"] * ispin + ["%.4f" % (0.5 * row)] * ispin
        doscar.append("  %.4f  %s" % (-5.0 + 15.0 * row / (nedos - 1), "  ".join(density)))
    Path("DOSCAR").write_text("\\n".join(doscar) + "\\n")
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
    down_shift: tuple[float, float] = (0.0, 0.0),
    metal: bool = False,
    tetrahedron_metal: bool = False,
) -> tuple[Workspace, str]:
    """Submit and run one DOS job, and return where it landed."""

    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    script = {"down_shift": list(down_shift), "metal": metal, "tetrahedron_metal": tetrahedron_metal}
    executable.write_text(_FAKE_VASP.format(script=script), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    stage: dict[str, Path] = {}
    for name, content in {"POSCAR": _POSCAR, "POTCAR": _POTCAR, **(extra_files or {})}.items():
        stage[name] = staged / name
        stage[name].write_text(content, encoding="utf-8")
    job = new_job(workspace, PACKAGE, files=stage, parameters=parameters or {}, name="packaged vasp-dos", install=True)
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=300.0)
    return workspace, job.job_id


def _mesh(kpoints: str) -> list[int]:
    return [int(n) for n in kpoints.splitlines()[3].split()]


def _edges(payload: Path) -> dict[str, object]:
    return json.loads((payload / "run" / "band_edges.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("publish_data", (False, True), ids=("default-workdir", "publish-data-opt-in"))
def test_the_dos_runs_scf_then_nscf_and_collects_the_gaps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, publish_data: bool
) -> None:
    from httk.core import FileRecord
    from httk.workflow import collect

    workspace, job_id = _campaign(
        tmp_path / "dos", monkeypatch, parameters={"publish_data": publish_data, "incar_tags": {"ISPIN": 1}}
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    scf, nscf = _runs(workdir)
    assert (scf["tags"]["ICHARG"], scf["tags"]["LCHARG"]) == ("2", ".TRUE.")
    assert (nscf["tags"]["ICHARG"], nscf["tags"]["LCHARG"], nscf["CHGCAR"]) == ("11", ".FALSE.", True)
    assert (nscf["tags"]["ISMEAR"], nscf["tags"]["NEDOS"], nscf["tags"]["LORBIT"]) == ("-5", "3001", "11")
    # Both meshes Gamma-centred; the nscf one denser.
    assert scf["KPOINTS"].splitlines()[2] == nscf["KPOINTS"].splitlines()[2] == "Gamma"
    assert all(d >= s for s, d in zip(_mesh(scf["KPOINTS"]), _mesh(nscf["KPOINTS"]), strict=True))
    assert sum(_mesh(nscf["KPOINTS"])) > sum(_mesh(scf["KPOINTS"]))
    assert scf["POSCAR"] == nscf["POSCAR"]

    edges = _edges(payload)
    assert edges["metallic"] is False and edges["given_density"] is False
    assert edges["indirect_gap"] == pytest.approx(1.2) and edges["direct_gap"] == pytest.approx(1.7)
    assert (edges["vbm_kpoint"], edges["cbm_kpoint"]) == (0, 2)
    assert edges["vbm_energy"] == pytest.approx(-0.5) and edges["energy_reference"] == 5.0
    assert (edges["ispin"], edges["nkpoints"], edges["nbands"], edges["occupation_tolerance"]) == (1, 3, 4, 1e-3)
    assert edges["occupations_out_of_range"] == 0
    assert "per_spin" not in edges
    state = _job_state(payload)
    assert (state["metallic"], state["fermi_energy"]) == (False, 5.0)
    assert state["band_gap"] == pytest.approx(1.2)

    if publish_data:
        published = _files(payload / "data")
        for name in ("scf/OUTCAR", "nscf/DOSCAR", "nscf/EIGENVAL", "nscf/band_edges.json"):
            assert name in published, name
    else:
        assert not (payload / "data").exists()

    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {
        *("total_energy", "electronic_density_of_states", "fermi_energy", "band_gap", "direct_band_gap"),
        *("dos", "eigenvalues", "band_edges"),
    }
    assert not item.unfulfilled
    assert item.outputs["total_energy"].value == pytest.approx(-10.8)
    assert item.outputs["band_gap"].band_gap == pytest.approx(1.2)
    assert item.outputs["direct_band_gap"].direct_band_gap == pytest.approx(1.7)
    assert item.outputs["fermi_energy"].fermi_energy == 5.0
    for role in ("dos", "eigenvalues", "band_edges"):
        assert isinstance(item.outputs[role], FileRecord)
    assert item.outputs["band_edges"].media_type == "application/json"
    assert item.outputs["dos"].media_type is None


def test_a_spin_polarized_gap_takes_both_channels(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from httk.workflow import collect

    # Spin down: valence 0.2 eV lower (VBM 4.3), conduction 0.4 eV lower (CBM 5.3).
    workspace, job_id = _campaign(tmp_path / "spin", monkeypatch, down_shift=(-0.2, -0.4))

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    edges = _edges(payload)
    assert edges["ispin"] == 2
    # min CBM over both spins (5.3, down) minus max VBM over both (4.5, up); band indices count up first.
    assert edges["indirect_gap"] == pytest.approx(0.8)
    assert (edges["vbm_band"], edges["cbm_band"]) == (1, 6)
    # Spin-conserving: down's 1.5 at k2, not the cross-spin 5.3 - 4.0 = 1.3.
    assert edges["direct_gap"] == pytest.approx(1.5)
    up, down = edges["per_spin"]
    assert (up["indirect_gap"], up["direct_gap"]) == (pytest.approx(1.2), pytest.approx(1.7))
    assert (down["indirect_gap"], down["direct_gap"]) == (pytest.approx(1.0), pytest.approx(1.5))
    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert item.outputs["band_gap"].band_gap == pytest.approx(0.8)
    assert item.outputs["direct_band_gap"].direct_band_gap == pytest.approx(1.5)


@pytest.mark.parametrize(("variant", "out_of_range"), (("metal", 0), ("tetrahedron_metal", 2 * 2)))
def test_a_metal_has_zero_gaps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, variant: str, out_of_range: int
) -> None:
    from httk.workflow import collect

    # tetrahedron_metal: occupations 1.0012 and -0.0017 (beyond the 1e-3 tolerance) next to the 0.5, in both spins.
    workspace, job_id = _campaign(tmp_path / variant, monkeypatch, **{variant: True})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    edges = _edges(payload)
    assert (edges["metallic"], edges["occupations_out_of_range"]) == (True, out_of_range)
    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert (item.outputs["band_gap"].band_gap, item.outputs["direct_band_gap"].direct_band_gap) == (0.0, 0.0)


def test_a_given_charge_density_skips_the_scf_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from httk.workflow import collect

    # The density, INCAR and OUTCAR of an earlier job's scf run, as a parent passes them.
    workspace, job_id = _campaign(tmp_path / "source", monkeypatch)
    _, source = _payload_of(workspace, job_id)
    given = {
        "CHGCAR": (source / "run" / "CHGCAR").read_text(encoding="utf-8"),
        "scf-INCAR": (source / "run" / "scf" / "INCAR").read_text(encoding="utf-8"),
        "scf-OUTCAR": (source / "run" / "scf" / "OUTCAR").read_text(encoding="utf-8"),
    }

    workspace, job_id = _campaign(tmp_path / "given", monkeypatch, extra_files=given)
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    (nscf,) = _runs(payload / "run")
    assert (nscf["tags"]["ICHARG"], nscf["CHGCAR"]) == ("11", True)
    assert (payload / "run" / "scf" / "OUTCAR").is_file()
    assert _edges(payload)["given_density"] is True
    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert "total_energy" not in item.outputs
    assert item.outputs["band_gap"].band_gap == pytest.approx(1.2)

    # The density's INCAR must be known.
    workspace, job_id = _campaign(tmp_path / "alone", monkeypatch, extra_files={"CHGCAR": given["CHGCAR"]})
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_missing"
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()

    # A density on another grid is refused before the nscf run.
    other = {**given, "scf-INCAR": re.sub(r"ENCUT\s*=\s*600", "ENCUT = 520", given["scf-INCAR"])}
    workspace, job_id = _campaign(tmp_path / "other", monkeypatch, extra_files=other)
    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.charge_density_incompatible"
    assert failure["details"]["differences"] == {"ENCUT": {"scf": "520", "nscf": "600"}}
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()


def test_an_icharg_override_fails_by_name_before_any_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "icharg", monkeypatch, parameters={"incar_tags": {"ICHARG": 2}})

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.invalid_dos_inputs"
    assert failure["details"] == {"source": "incar_tags", "ICHARG": "2"}
    assert not (payload / "run" / "fake-vasp-runs.jsonl").exists()


def test_the_templates_keep_the_scf_density_settings() -> None:
    from httk.codes.vasp import read_incar

    # A copy of vasp-band-structure's, kept identical here.
    scf_text = (PACKAGE / "templates" / "INCAR.scf").read_bytes()
    assert scf_text == (REPO_ROOT / "vasp-band-structure" / "templates" / "INCAR.scf").read_bytes()
    scf = read_incar(PACKAGE / "templates" / "INCAR.scf")
    nscf = read_incar(PACKAGE / "templates" / "INCAR.nscf")
    density = runpy.run_path(str(PACKAGE / "run.py"))["DENSITY_TAGS"]
    assert {tag: scf.get(tag) for tag in density} == {tag: nscf.get(tag) for tag in density}
    assert (nscf["ICHARG"], nscf["LCHARG"], nscf["ISMEAR"], nscf["NEDOS"]) == ("11", ".FALSE.", "-5", "3001")
