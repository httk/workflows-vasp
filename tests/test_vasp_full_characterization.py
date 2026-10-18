"""The ``vasp-full-characterization`` parent, end to end through the real manager, against one synthetic material.

One model VASP answers every run of every child with the same material: a diamond-structure
crystal of conventional edge ``A`` whose HPKOT primitive cell ``L0`` the test computes up front,
as the parent does. Each run writes an OUTCAR the real ``OutcarFile`` parses, an OSZICAR, a
CONTCAR equal to its POSCAR (so the relaxation keeps the structure and converges), an EIGENVAL
with occupations (a 1 eV indirect gap, or a partially occupied band for a metal) and a DOSCAR,
and a CHGCAR when ``LCHARG`` is true. The energy is a Birch–Murnaghan curve of the volume per
atom; an ``ISIF = 2`` run reports the stress ``C e`` of a cubic stiffness for its strain
against ``L0``; a single point (``NSW = 0``) puts the Einstein force ``F = -k u`` on every atom
(``u`` its displacement from the 1/8 grid of the primitive cell, as in ``test_vasp_phonons.py``);
a ``LEPSILON`` run appends the dielectric and Born charge blocks. ``fail_on`` makes every run
that sets that tag fail.

Every child is a real job of a package installed into the test workspace from this repository.
Each job's runs are logged in its own workdir, so the test sees which calculations ran where.
"""

import json
import time
from pathlib import Path
from typing import Any

import pytest

pytest.importorskip("httk.atomistic")
pytest.importorskip("httk.analyse")
pytest.importorskip("phonopy")
pytest.importorskip("spglib")

from httk.workflow import TaskManager, Workspace, _kernel, _store, collect
from httk.workflow.introspection import resolve_job
from httk.workflow.scaffold import new_job

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-full-characterization"
# Callees before callers, so every package's declared calls resolve when it is installed.
CALLEES = (
    "vasp-relax",
    "vasp-relax-ions",
    "vasp-static",
    "vasp-band-structure",
    "vasp-dos",
    "vasp-dielectric",
    "vasp-elastic",
    "vasp-eos",
    "vasp-phonons",
)
A = 5.43
C = f"{A:.16f}"
GPA_PER_EV_PER_A3 = 160.21766208
_POTCAR = "  TITEL  = PAW_PBE Si 05Jan2001\n   ZVAL   =    4.000    mass and valenz\n"
# The conventional cubic cell of diamond silicon, its sites off by up to 2e-5 (fractional, 1e-4 Angstrom):
# within the parent's symmetry tolerance, so its HPKOT cell is the clean two-atom primitive cell.
_NOISE = (1e-5, -2e-5, 0.0, 2e-5, -1e-5, 1e-5, 0.0, -1e-5)
_SITES = ("0 0 0", "0 .5 .5", ".5 0 .5", ".5 .5 0", ".25 .25 .25", ".25 .75 .75", ".75 .25 .75", ".75 .75 .25")
_CONVENTIONAL = f"Si conventional\n1.0\n{C} 0 0\n0 {C} 0\n0 0 {C}\nSi\n8\nDirect\n" + "".join(
    " ".join(f"{float(x) + d:.6f}" for x in site.split()) + "\n" for site, d in zip(_SITES, _NOISE, strict=True)
)
# Zincblende SiC on the same lattice: Si on the first four sites, C on the last four. Re-standardizing its HPKOT
# cell gives another setting of it, so the band-structure and phonons children see other bases than the parent.
_SIC = _CONVENTIONAL.replace("Si conventional", "SiC conventional").replace("Si\n8\n", "Si C\n4 4\n")
_POTCAR_SIC = _POTCAR + "  TITEL  = PAW_PBE C 08Apr2002\n   ZVAL   =    4.000    mass and valenz\n"
_C11, _C12, _C44 = 160.0, 60.0, 80.0
_STIFFNESS = [
    [_C11 if i == j else _C12 if i < 3 and j < 3 else _C44 if i == j else 0.0 for j in range(6)] for i in range(6)
]
_B0 = 100.0  # GPa

_FAKE_VASP = '''#!/usr/bin/env python3
"""One synthetic material for every run of the characterization (see the test module)."""

import json
import math
import sys
from pathlib import Path

import numpy as np

MODEL = json.loads(MODEL_JSON)

tags = {}
for line in Path("INCAR").read_text().splitlines():
    for part in line.split("#")[0].split("!")[0].split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            tags[key.strip().upper()] = value.strip()
poscar = Path("POSCAR").read_text()
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps({"tags": tags, "POSCAR": poscar, "CHGCAR": Path("CHGCAR").is_file()}) + "\\n")


def true(tag):
    return tags.get(tag, "").strip(".").upper().startswith("T")


if MODEL["fail_on"] is not None and MODEL["fail_on"] in tags:
    print("the model VASP refuses this run")
    Path("OUTCAR").write_text(" vasp.6.4.1 fake\\n")
    sys.exit(1)
lines = poscar.splitlines()
scale = float(lines[1].split()[0])
lattice = np.array([[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
counts = [int(n) for n in lines[6].split()]
ions = sum(counts)
fractional = np.array([[float(x) for x in lines[8 + i].split()[:3]] for i in range(ions)])
volume = abs(np.linalg.det(lattice))
x = (MODEL["v0"] / (volume / ions)) ** (2 / 3) - 1
energy = ions * (MODEL["e0"] + 9 * MODEL["v0"] * MODEL["b0"] / 16 * (2 * x * x + (MODEL["b1"] - 4) * x**3))
nsw = int(tags.get("NSW", "0"))
kbar = [0.0] * 6
if tags.get("ISIF") == "2":
    # Rows are the lattice vectors, so L' = L0 F^T.
    gradient = np.linalg.solve(np.array(MODEL["reference"]), lattice).T
    eps = (gradient + gradient.T) / 2 - np.eye(3)
    strain = [eps[0, 0], eps[1, 1], eps[2, 2], 2 * eps[1, 2], 2 * eps[0, 2], 2 * eps[0, 1]]
    sigma = np.array(MODEL["stiffness"]) @ strain
    kbar = [-10 * sigma[i] for i in (0, 1, 2, 5, 3, 4)]
forces = np.zeros((ions, 3))
if nsw == 0:
    # The supercell multiples: each lattice row over the fcc primitive length; u is what rounding to the 1/8 grid
    # of the primitive cell takes away.
    multiples = np.array([round(np.linalg.norm(row) / (MODEL["a"] / math.sqrt(2))) for row in lattice])
    scaled = fractional * multiples
    du = (scaled - np.round(scaled * 8) / 8) / multiples
    forces = -MODEL["k"] * (du @ lattice)
cartesian = fractional @ lattice
rows = "".join(
    "  %12.5f %12.5f %12.5f    %14.8f %14.8f %14.8f\\n" % (*r, *f) for r, f in zip(cartesian, forces)
)
outcar = (
    " vasp.6.4.1 fake\\n"
    "   ions per type =  %s\\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\\n"
    "   NSW    =     %d    number of steps for IOM\\n"
    "--------------------------------------- Iteration      1(   1)  ---------------------------------------\\n"
    "  FORCE on cell =-STRESS in cart. coord.  units (eV):\\n"
    "  in kB %s\\n"
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
    " General timing and accounting informations for this job:\\n"
    % ("  ".join(map(str, counts)), nsw, " ".join("%16.8f" % v for v in kbar), rows, energy, energy, energy)
)
if true("LEPSILON"):
    # Isotropic, +Z on the first species and -Z on the second (diamond Si: Z = 0, as its symmetry requires).
    charges = [MODEL["born"] if s == 0 else -MODEL["born"] for s, n in enumerate(counts) for _ in range(n)]
    born = "".join(
        " ion %4d\\n    1 %11.5f     0.00000     0.00000\\n    2     0.00000 %11.5f     0.00000\\n"
        "    3     0.00000     0.00000 %11.5f\\n" % (i + 1, z, z, z)
        for i, z in enumerate(charges)
    )
    outcar += (
        "\\n MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)\\n"
        " ------------------------------------------------------\\n"
        "          12.000000     0.000000     0.000000\\n"
        "           0.000000    12.000000     0.000000\\n"
        "           0.000000     0.000000    12.000000\\n"
        " ------------------------------------------------------\\n\\n"
        " MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION\\n"
        " ------------------------------------------------------\\n"
        "           0.000000     0.000000     0.000000\\n"
        "           0.000000     0.000000     0.000000\\n"
        "           0.000000     0.000000     0.000000\\n"
        " ------------------------------------------------------\\n\\n"
        " BORN EFFECTIVE CHARGES (in e, cummulative output)\\n"
        " -------------------------------------------------\\n" + born
    )
Path("OUTCAR").write_text(outcar)
Path("OSZICAR").write_text(
    "       N       E                     dE             d eps       ncg     rms\\n"
    "DAV:   1    %.12E   -0.50000E+00   -0.10000E+00   128   0.100E+00\\n"
    "   1 F= %.8E E0= %.8E  d E =-.105000E+02\\n" % (energy, energy, energy)
)
Path("CONTCAR").write_text(poscar)
Path("vasprun.xml").write_text('<modeling><structure name="finalpos"></structure></modeling>\\n')
# Every run: three k-points, four bands, VBM 0 at k0 and CBM 1 at k2 (a 1 eV indirect gap), or a metal.
ispin = int(tags.get("ISPIN", "1"))
bands = [[-5.0, -5.0, -5.0], [0.0, -0.2, -0.4], [1.5, 1.2, 1.0], [3.0, 3.0, 3.0]]
occupations = [[1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0]]
if MODEL["metal"]:
    # A partially occupied state at the Fermi level, as the tetrahedron method (ISMEAR = -5) writes them: outside
    # [0, 1].
    occupations[1][2], occupations[2][2] = 1.0012, -0.0015
eigenval = ["  %d  %d    1    %d" % (ions, ions, ispin), "  0.4E+02  0.3E-09  0.3E-09  0.3E-09  0.5E-15", "  1.0E-004"]
eigenval += ["  CAR", " Si", "      8  3  4"]
for k in range(3):
    eigenval += ["", "  %.7E  %.7E  %.7E  %.7E" % (k / 4, 0.0, 0.0, 1 / 3)]
    for band in range(4):
        values = [bands[band][k]] * ispin + [occupations[band][k]] * ispin
        eigenval.append("    %d  " % (band + 1) + "  ".join("%.6f" % value for value in values))
Path("EIGENVAL").write_text("\\n".join(eigenval) + "\\n")
nedos = int(tags.get("NEDOS", "301"))
doscar = ["  %d  %d    1    0" % (ions, ions), "  0.4E+02  0.3E-09  0.3E-09  0.3E-09  0.5E-15", "  1.0E-004", "  CAR"]
doscar += [" Si", "  10.00000000  -5.00000000  %d  %.8f  1.00000000" % (nedos, 0.5)]
for row in range(nedos):
    density = ["0.5000"] * ispin + ["%.4f" % (0.5 * row)] * ispin
    doscar.append("  %.4f  %s" % (-5.0 + 15.0 * row / (nedos - 1), "  ".join(density)))
Path("DOSCAR").write_text("\\n".join(doscar) + "\\n")
if true("LCHARG"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''


def _standardized_lattice(root: Path, conventional: str) -> list[list[float]]:
    """The HPKOT cell of the input, as the parent's plan step builds it."""

    import httk.atomistic
    import httk.core

    root.mkdir(parents=True)
    (root / "POSCAR").write_text(conventional, encoding="utf-8")
    structure = httk.core.load(str(root / "POSCAR"), format="vasp-poscar", precision=1e-3)
    httk.core.save(httk.atomistic.band_path(structure, tolerance=1e-3).structure, str(root / "HPKOT.vasp"))
    return _lattice((root / "HPKOT.vasp").read_text(encoding="utf-8"))


def _lattice(poscar: str) -> list[list[float]]:
    lines = poscar.splitlines()
    return [[float(lines[1]) * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    parameters: dict[str, object],
    *,
    metal: bool = False,
    fail_on: str | None = None,
    sic: bool = False,
) -> tuple[Workspace, str, list[list[float]]]:
    """Run the parent until idle; return the workspace, its job id and the standardized lattice."""

    conventional, potcar = (_SIC, _POTCAR_SIC) if sic else (_CONVENTIONAL, _POTCAR)
    reference = _standardized_lattice(root / "hpkot", conventional)
    # The volume of the HPKOT cell per atom, a little off the model's minimum.
    v0 = _volume(reference)
    model = {
        "reference": reference,
        "stiffness": _STIFFNESS,
        "a": A,
        "k": 10.0,
        "v0": v0 / 2 * 1.01,
        "e0": -5.4,
        "b0": _B0 / GPA_PER_EV_PER_A3,
        "b1": 4.5,
        "metal": metal,
        "fail_on": fail_on,
        "born": 2.0 if sic else 0.0,
    }
    executable = root / "fake-vasp"
    executable.write_text(_FAKE_VASP.replace("MODEL_JSON", repr(json.dumps(model))), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    files = {}
    for name, content in (("POSCAR", conventional), ("POTCAR", potcar)):
        (staged / name).write_text(content, encoding="utf-8")
        files[name] = staged / name
    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        for callee in CALLEES:
            _store.install(workspace, owner, REPO_ROOT / callee)
    finally:
        owner.close()
    job = new_job(
        workspace, PACKAGE, files=files, parameters=parameters, name="vasp-full-characterization", install=True
    )
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=1800.0)
    return workspace, job.job_id, reference


def _jobs(workspace: Workspace) -> dict[str, tuple[Path, dict[str, Any]]]:
    """Every job of the workspace: job id -> (job directory, job.json)."""

    jobs = {}
    for path in workspace.root.glob("jobs/*/*/job.json"):
        document = json.loads(path.read_text(encoding="utf-8"))
        jobs[document["id"]] = (path.parent, document)
    return jobs


def _children(jobs: dict[str, tuple[Path, dict[str, Any]]], job_id: str) -> dict[str, tuple[Path, dict[str, Any]]]:
    """The children of one job by label."""

    return {
        document["tag"]: (path, document)
        for path, document in jobs.values()
        if (document.get("parent") or {}).get("job_id") == job_id
    }


def _workflow(document: dict[str, Any]) -> str:
    return document["workflow"]["name"]


def _runs(directory: Path) -> list[dict[str, Any]]:
    path = directory / "run" / "fake-vasp-runs.jsonl"
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _failure(workspace: Workspace, job_id: str) -> dict[str, Any]:
    ref = resolve_job(workspace, job_id)
    assert ref.state == "failed", ref.state
    return json.loads((ref.path / "state.json").read_text(encoding="utf-8"))["failure"]


def _result(workspace: Workspace, job_id: str) -> dict[str, Any]:
    return json.loads(
        (resolve_job(workspace, job_id).path / "run" / "characterization.json").read_text(encoding="utf-8")
    )


def _volume(m: list[list[float]]) -> float:
    return abs(
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _same_cell(poscar: str, lattice: list[list[float]]) -> bool:
    """Whether *poscar* is a two-atom cell of *lattice*'s volume: the standardized cell in any setting."""

    counts = sum(int(n) for n in poscar.splitlines()[6].split())
    return counts == 2 and abs(_volume(_lattice(poscar)) - _volume(lattice)) < 1e-6


# Case -> (metal, fail_on, SiC, normal profile?).
_CASES = {
    "insulator": (False, None, False, True),
    "dielectric_failed": (False, "LEPSILON", False, True),
    "metal": (True, None, False, False),
    "polar_insulator": (False, None, True, False),
}


@pytest.mark.parametrize("case", tuple(_CASES))
def test_full_characterization(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, test_profile, case: str) -> None:
    metal, fail_on, sic, normal = _CASES[case]
    if not normal and not test_profile.extended:
        pytest.skip("extended profile only (HTTK_TEST_PROFILE=extended)")
    extended = test_profile.extended
    parameters: dict[str, object] = {
        "maximum_remedies": 2,
        # Normal: two strain magnitudes (12 cells), five volumes, a 2x2x2 supercell; extended: the defaults, 3x3x3.
        "elastic_parameters": {} if extended else {"strain_magnitudes": [-0.01, 0.01]},
        "eos_parameters": {} if extended else {"number_of_frames": 5},
        "phonons_parameters": {"min_length": 9.0 if extended else 6.0},
        # A failing child fails at once, without remedies.
        "dielectric_parameters": {"maximum_remedies": 0},
    }
    started = time.monotonic()
    workspace, job_id, reference = _campaign(
        tmp_path / case, monkeypatch, parameters, metal=metal, fail_on=fail_on, sic=sic
    )
    elapsed = time.monotonic() - started
    jobs = _jobs(workspace)
    children = _children(jobs, job_id)
    result = _result(workspace, job_id)
    decisions, parts = result["decisions"], result["parts"]
    expected = {"relax", "band_structure", "dos", "elastic", "eos", "phonons"} | (
        {"dielectric"} if not metal else set()
    )
    assert set(children) == set(parts) == expected
    assert {label: _workflow(document) for label, (_, document) in children.items()} == {
        label: "vasp." + label.replace("_", "-") for label in expected
    }
    assert decisions["insulator"] is (not metal) and decisions["band_gap_source"] == "relax EIGENVAL"
    assert decisions["band_gap"] == pytest.approx(0.0 if metal else 1.0)
    # One relaxation in the whole workspace: elastic, eos and phonons ran none of their own.
    assert [_workflow(document) for _, document in jobs.values()].count("vasp.relax") == 1
    for label in ("elastic", "eos", "phonons"):
        assert children[label][1]["parameters"]["relax"] is False
        assert children[label][1]["parameters"]["maximum_remedies"] == 2
    assert children["phonons"][1]["parameters"]["standardize"] is False
    if not metal:
        # <label>_parameters win over the shared parameters.
        assert children["dielectric"][1]["parameters"]["maximum_remedies"] == 0
    assert children["elastic"][1]["parameters"].get("strain_magnitudes") == (None if extended else [-0.01, 0.01])

    # Exactly one self-consistent ground-state run on the standardized cell in the whole tree: band structure's.
    scf = [
        (path, run)
        for path, _ in jobs.values()
        for run in _runs(path)
        if run["tags"].get("NSW", "0") == "0"
        and run["tags"].get("ICHARG") != "11"
        and "LEPSILON" not in run["tags"]
        and _same_cell(run["POSCAR"], reference)
    ]
    assert [path for path, _ in scf] == [children["band_structure"][0]]
    # dos ran VASP once, the nscf stage, on the band-structure child's density.
    (dos_run,) = _runs(children["dos"][0])
    assert dos_run["tags"]["ICHARG"] == "11" and dos_run["CHGCAR"]
    edges = json.loads((children["dos"][0] / "run" / "band_edges.json").read_text(encoding="utf-8"))
    assert edges["given_density"] is True and edges["metallic"] is metal
    phonons = json.loads((children["phonons"][0] / "run" / "phonons.json").read_text(encoding="utf-8"))
    assert decisions["dos_band_gap"]["band_gap"] == pytest.approx(0.0 if metal else 1.0)

    if case == "dielectric_failed":
        failure = _failure(workspace, job_id)
        # The phonons succeeded, but without the correction a polar insulator needs: incomplete too.
        assert failure["code"] == "vasp.characterization_incomplete" and set(failure["details"]["failed"]) == {
            "dielectric",
            "phonons (without the non-analytical correction: dielectric did not succeed)",
        }
        assert parts["dielectric"]["kind"] == "failed" and parts["phonons"]["kind"] == "succeeded"
        assert phonons["nac"] == {"used": False} and decisions["phonons_nac"] is False
        assert "dielectric failed: phonons without non-analytical correction" in decisions["notes"]
        return
    assert resolve_job(workspace, job_id).state == "succeeded", _failure(workspace, job_id)
    assert all(entry["kind"] == "succeeded" for entry in parts.values())
    assert phonons["nac"]["used"] is (not metal)
    # A metal's phonons need no correction: the job succeeds without it.
    assert decisions["phonons_nac"] is (not metal)
    # The metal's relaxation EIGENVAL has tetrahedron occupations outside [0, 1]; they count as partial.
    assert (decisions["occupations_out_of_range"] > 0) is metal
    if metal:
        assert "metal: no dielectric child" in decisions["notes"]
    else:
        assert len(_runs(children["dielectric"][0])) == 1

    items = {item.record.job_id: item for item in collect(workspace, fail_fast=True, allow_job_collector=True)}
    assert set(items[job_id].outputs) == {"relaxed_structure", "standardized_structure", "characterization"}
    assert items[job_id].outputs["characterization"].media_type == "application/json"
    outputs = {label: items[document["id"]].outputs for label, (_, document) in children.items()}
    assert outputs["dos"]["band_gap"].band_gap == pytest.approx(0.0 if metal else 1.0)
    assert outputs["eos"]["bulk_modulus"].value == pytest.approx(_B0, rel=1e-3)
    assert outputs["elastic"]["bulk_modulus_voigt"].value == pytest.approx((_C11 + 2 * _C12) / 3, abs=1e-3)
    assert "elastic_tensor" in outputs["elastic"] and "zero_point_energy" in outputs["phonons"]
    if not metal:
        assert outputs["dielectric"]["high_frequency_relative_permittivity"].value == pytest.approx(12.0)
    print(f"vasp-full-characterization {case} ({test_profile.name}): {len(jobs)} jobs in {elapsed:.1f} s")


def test_parts_eos_runs_only_the_relaxation_and_the_eos(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id, _ = _campaign(
        tmp_path / "eos", monkeypatch, {"parts": "eos", "eos_parameters": {"number_of_frames": 5}}
    )
    assert resolve_job(workspace, job_id).state == "succeeded", _failure(workspace, job_id)
    jobs = _jobs(workspace)
    assert set(_children(jobs, job_id)) == {"relax", "eos"}
    workflows = sorted(_workflow(document) for _, document in jobs.values())
    assert workflows == sorted(["vasp.full-characterization", "vasp.relax", "vasp.eos", *["vasp.relax-ions"] * 5])


@pytest.mark.parametrize(
    ("parameters", "code"),
    (
        ({"parts": "eos bands"}, "vasp.input_invalid"),
        ({"dos_parameters": "fast"}, "vasp.input_invalid"),
        ({"incar_tags": {"MAGMOM": "8*1"}}, "vasp.cell_dependent_tag_unsupported"),
        ({"phonons_parameters": {"incar_tags": {"nupdown": "0"}}}, "vasp.cell_dependent_tag_unsupported"),
    ),
)
def test_refused_inputs_fail_before_any_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict[str, object], code: str
) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "refused", monkeypatch, parameters)
    assert _failure(workspace, job_id)["code"] == code
    assert len(_jobs(workspace)) == 1
