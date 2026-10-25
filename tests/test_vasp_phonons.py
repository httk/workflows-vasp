"""The ``vasp-phonons`` composite, end to end through the real manager, against an Einstein-model VASP.

The model VASP of a single point (``NSW = 0``) puts a force ``F = -k u`` on every atom, ``u``
its displacement from the ideal site. The ideal sites are recovered from the POSCAR alone: the
test cells are fcc lattices of conventional edge ``A`` (any HPKOT standardization of them keeps
that primitive length), whose sites all lie on the 1/8 grid of the primitive cell's fractional
coordinates, so the supercell multiples are the lattice row lengths over the primitive length and
``u`` is what rounding the unit-cell fractional coordinates to that grid takes away. The force
constants are then the self term ``k`` alone, so every branch is flat at ``sqrt(k / m)``, which
the real OUTCAR parser and the real phonopy must recover. A relaxation run reports zero forces
and keeps its POSCAR as CONTCAR. The composite calls ``vasp-relax`` and ``vasp-static`` as child
jobs, so those packages are installed into the test workspace from this repository.
"""

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import numpy
import pytest
from httk.workflow import TaskManager, Workspace, _kernel, _store
from httk.workflow.introspection import list_jobs, resolve_job
from httk.workflow.scaffold import new_job

from test_vasp_runners import _POTCAR

pytest.importorskip("phonopy")
pytest.importorskip("spglib")
pytest.importorskip("httk.analyse")

from phonopy.physical_units import get_physical_units
from phonopy.structure.atoms import PhonopyAtoms

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-phonons"
CALLEES = ("vasp-relax", "vasp-static")
EXCLUDED = "bump_kpoints equal_kpoints centering incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice"
A = 5.43
H_EV_PER_THZ = 4.135667696e-3
K = 10.0  # eV/A^2
# Written in full, as VASP and httk write POSCARs: the reader takes a short number as uncertain in its last digit.
H = f"{A / 2:.16f}"
C = f"{A:.16f}"

_SI = f"""Si diamond, HPKOT
1.0
0.0 {H} {H}
{H} 0.0 {H}
{H} {H} 0.0
Si
2
Direct
0.125 0.125 0.125
0.875 0.875 0.875
"""
_SIC = f"""SiC zincblende, HPKOT
1.0
0.0 {H} {H}
{H} 0.0 {H}
{H} {H} 0.0
Si C
1 1
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""
# The conventional cubic cell of diamond silicon, its sites off by up to 2e-5 (fractional, 1.1e-4 A; at 5e-5 the
# recognition at symmetry_tolerance 1e-3 already finds a lower symmetry).
_NOISE = (2e-5, -2e-5, 0.0, 2e-5, -1e-5, 1e-5, 0.0, -2e-5)
_SITES = ("0 0 0", "0 .5 .5", ".5 0 .5", ".5 .5 0", ".25 .25 .25", ".25 .75 .75", ".75 .25 .75", ".75 .75 .25")
_CONVENTIONAL = f"Si conventional\n1.0\n{C} 0 0\n0 {C} 0\n0 0 {C}\nSi\n8\nDirect\n" + "".join(
    " ".join(f"{float(x) + (-d if k == 1 else d):.8f}" for k, x in enumerate(site.split())) + "\n"
    for site, d in zip(_SITES, _NOISE, strict=True)
)
# The HPKOT vectors permuted cyclically: a unimodular change of basis of the same lattice, still right-handed (a
# left-handed one is negated, and httk-atomistic's OUTCAR reader drops POSITION rows that start with a minus sign).
_PERMUTED = _SI.replace(f"0.0 {H} {H}\n{H} 0.0 {H}\n{H} {H} 0.0\n", f"{H} 0.0 {H}\n{H} {H} 0.0\n0.0 {H} {H}\n")
_POTCAR_SIC = _POTCAR + "  TITEL  = PAW_PBE C 08Apr2002\n   ZVAL   =    4.000    mass and valenz\n"

_FAKE_VASP = '''#!/usr/bin/env python3
"""An Einstein-model VASP: F = -k u on every atom of a single point, zero forces otherwise."""

import json
import math
from fractions import Fraction
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
    stream.write(json.dumps({"tags": tags, "POSCAR": poscar}) + "\\n")
lines = poscar.splitlines()
scale = float(lines[1].split()[0])
m = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
species = lines[5].split()
counts = [int(n) for n in lines[6].split()]
symbols = [s for s, n in zip(species, counts) for _ in range(n)]
fractional = [[float(x) for x in lines[8 + i].split()[:3]] for i in range(len(symbols))]
# The supercell multiples: each lattice row over the fcc primitive length.
multiples = [round(math.sqrt(sum(x * x for x in row)) / (MODEL["a"] / math.sqrt(2))) for row in m]
nsw = int(tags.get("NSW", "0"))
rows = []
for symbol, s in zip(symbols, fractional):
    du = [0.0, 0.0, 0.0]
    if nsw == 0:
        for i in range(3):
            x = s[i] * multiples[i]
            du[i] = (x - round(x * 8) / 8) / multiples[i]
    u = [sum(du[i] * m[i][j] for i in range(3)) for j in range(3)]
    r = [sum(s[i] * m[i][j] for i in range(3)) for j in range(3)]
    f = [-MODEL["k"][symbol] * x for x in u]
    rows.append("  %12.5f %12.5f %12.5f    %14.8f %14.8f %14.8f\\n" % (*r, *f))
energy = -10.8
Path("OUTCAR").write_text(
    " vasp.6.4.1 fake\\n"
    "   ions per type =  %s\\n"
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
    " General timing and accounting informations for this job:\\n"
    % ("  ".join(map(str, counts)), nsw, "".join(rows), energy, energy, energy)
)
Path("OSZICAR").write_text(
    "       N       E                     dE             d eps       ncg     rms\\n"
    "DAV:   1    %.12E   -0.50000E+00   -0.10000E+00   128   0.100E+00\\n"
    "   1 F= %.8E E0= %.8E  d E =-.105000E+02\\n" % (energy, energy, energy)
)
Path("CONTCAR").write_text(poscar)
Path("vasprun.xml").write_text('<modeling><structure name="finalpos"></structure></modeling>\\n')
if tags.get("LCHARG", "").strip(".").upper().startswith("T"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    parameters: dict[str, object],
    *,
    poscar: str = _SI,
    k: dict[str, float] | None = None,
    dielectric: dict[str, Any] | None = None,
) -> tuple[Workspace, str]:
    root.mkdir(parents=True)
    executable = root / "fake-vasp"
    model = repr(json.dumps({"a": A, "k": k or {"Si": K, "C": K}}))
    executable.write_text(_FAKE_VASP.replace("MODEL_JSON", model), encoding="utf-8")
    executable.chmod(0o755)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    files = {}
    contents = {"POSCAR": poscar, "POTCAR": _POTCAR_SIC if " C\n" in poscar else _POTCAR}
    if dielectric is not None:
        contents["dielectric.json"] = json.dumps(dielectric)
    for name, content in contents.items():
        (staged / name).write_text(content, encoding="utf-8")
        files[name] = staged / name
    # The callees first, from this repository, so the composite's declared calls resolve to them.
    owner = _kernel.register_owner(workspace, kind="cli", label="test", allocation=None, advertised={})
    try:
        for callee in CALLEES:
            _store.install(workspace, owner, REPO_ROOT / callee)
    finally:
        owner.close()
    job = new_job(workspace, PACKAGE, files=files, parameters=parameters, name="packaged vasp-phonons", install=True)
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=600.0)
    return workspace, job.job_id


def _job(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _phonons(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "run" / "phonons.json").read_text(encoding="utf-8"))


def _runs(payload: Path) -> list[dict[str, Any]]:
    path = payload / "run" / "fake-vasp-runs.jsonl"
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _failure(workspace: Workspace, job_id: str) -> dict[str, Any]:
    state, payload = _job(workspace, job_id)
    assert state == "failed"
    return json.loads((payload / "state.json").read_text(encoding="utf-8"))["failure"]


def _einstein(symbol: str, k: float = K) -> float:
    """The Einstein frequency sqrt(k / m) in THz, with phonopy's mass and VASP unit factor."""

    (mass,) = PhonopyAtoms(
        symbols=[symbol], cell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], scaled_positions=[[0, 0, 0]]
    ).masses
    return math.sqrt(k / mass) * get_physical_units().DefaultToTHz


def _sic_dielectric(poscar: str, *, lattice_scale: float = 1.0, born: float | None = 2.0) -> dict[str, Any]:
    """A ``dielectric.json`` of the cell of *poscar*, as ``vasp.dielectric`` writes it."""

    lines = poscar.splitlines()
    lattice = [[lattice_scale * float(x) for x in lines[i].split()] for i in (2, 3, 4)]
    z = [[2.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    return {
        "electronic": [[10.0 if i == j else 0.0 for j in range(3)] for i in range(3)],
        "ionic": None,
        "static": None,
        "born_charges": None if born is None else [z, [[-x for x in row] for row in z]],
        "species": ["Si", "C"],
        "counts": [1, 1],
        "poscar_sha256": "0" * 64,
        "lattice": lattice,
        "positions": [[float(x) for x in lines[i].split()] for i in (8, 9)],
        "energy": -10.8,
    }


def _standardized(poscar: str, root: Path) -> str:
    """Return *poscar* standardized as the runner's ``standardize = true`` path writes it."""

    import httk.atomistic
    import httk.core

    root.mkdir(parents=True)
    (root / "POSCAR").write_text(poscar, encoding="utf-8")
    structure = httk.core.load(str(root / "POSCAR"), format="vasp-poscar", precision=1e-3)
    httk.core.save(httk.atomistic.band_path(structure, tolerance=1e-3).structure, str(root / "HPKOT.vasp"))
    return (root / "HPKOT.vasp").read_text(encoding="utf-8")


# primitive: an HPKOT input, standardized again; relax: the same after the tight relaxation; conventional: a noisy
# conventional cell; standardized: a cell this package's standardize path wrote, given with standardize = false;
# permuted: an HPKOT cell with its vectors permuted, given with standardize = false (a unimodular change of basis).
@pytest.mark.parametrize("case", ("primitive", "relax", "conventional", "standardized", "permuted"))
def test_the_einstein_model_gives_flat_bands_and_its_zero_point_energy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, test_profile, case: str
) -> None:
    from httk.workflow import collect

    # 2x2x2 (16 atoms), extended 3x3x3 (54 atoms), of the two-atom cell (perpendicular width 3.135 A).
    min_length = 9.0 if test_profile.extended else 6.0
    standardize = case not in ("standardized", "permuted")
    parameters = {"relax": case == "relax", "standardize": standardize, "min_length": min_length}
    poscar = {"conventional": _CONVENTIONAL, "permuted": _PERMUTED}.get(case, _SI)
    if case == "standardized":
        poscar = _standardized(_CONVENTIONAL, tmp_path / "hpkot")
    workspace, job_id = _campaign(tmp_path / case, monkeypatch, parameters, poscar=poscar)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id)
    phonons = _phonons(payload)
    multiple = 3 if test_profile.extended else 2
    assert phonons["supercell_matrix"] == [[multiple if i == j else 0 for j in range(3)] for i in range(3)]
    assert phonons["supercell_atoms"] == 2 * multiple**3
    # The two-atom HPKOT cell, also from the noisy conventional cell.
    assert phonons["unit_cell"]["counts"] == [2]
    assert phonons["band_path"]["bravais_lattice_extended"] == "cF2"
    hpkot = ([0, A / 2, A / 2], [A / 2, 0, A / 2], [A / 2, A / 2, 0])
    if case == "permuted":
        hpkot = (hpkot[1], hpkot[2], hpkot[0])
    for row, expected in zip(phonons["unit_cell"]["lattice"], hpkot, strict=True):
        assert row == pytest.approx(expected, abs=1e-6)
    band_path = json.loads((payload / "run" / "unitcell" / "band_path.json").read_text(encoding="utf-8"))
    # band_path's cell is the HPKOT basis: L_path = T L_unit. band_path keeps an input that already is a standard
    # primitive cell, and the cyclic permutation of the cubic HPKOT vectors is one (a symmetry of the cube).
    assert band_path["basis_change"] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

    nu = _einstein("Si")
    # Flat bands at the Einstein frequency, on the path and on the mesh (phonopy diagonalizes in float64).
    bands = phonons["band_structure"]
    assert len(bands["frequencies"]) == 6
    assert all(f == pytest.approx(nu, abs=1e-6) for mode in bands["frequencies"] for f in mode)
    assert len(bands["qpoints"]) == 51 * len(bands["qpoint_labels"]) // 2
    assert bands["qpoint_labels"][:2] == ["GAMMA", "X"] and bands["qpoint_label_indices"][:2] == [0, 50]
    assert phonons["imaginary_modes"]["count"] == 0
    assert phonons["imaginary_modes"]["lowest_frequency"] == pytest.approx(nu, abs=1e-6)
    dos = phonons["dos"]
    assert dos["frequencies"][dos["density"].index(max(dos["density"]))] == pytest.approx(nu, abs=1e-6)
    thermodynamics = phonons["thermodynamics"]
    assert thermodynamics["zero_point_energy"] == pytest.approx(3 * 2 * H_EV_PER_THZ * nu / 2, abs=1e-6)
    assert thermodynamics["temperatures"] == [float(t) for t in range(0, 1001, 10)]
    assert phonons["nac"] == {"used": False}

    # One vasp.static per displacement, on one mesh, with the phonon INCAR and the excluded remedies.
    assert len(phonons["children"]) == len(phonons["displacements"]) == 1
    for record in phonons["children"]:
        _, child = _job(workspace, record["job_id"])
        (run,) = _runs(child)
        tags = run["tags"]
        assert (tags["NSW"], tags["IBRION"], tags["ISMEAR"], tags["SIGMA"], tags["EDIFF"]) == (
            "0",
            "-1",
            "0",
            "0.05",
            "1E-7",
        )
        assert (tags["LREAL"], tags["ADDGRID"], tags["LCHARG"], tags["LWAVE"]) == (
            ".FALSE.",
            ".TRUE.",
            ".FALSE.",
            ".FALSE.",
        )
        job = json.loads((child / "job.json").read_text(encoding="utf-8"))
        assert job["parameters"]["excluded_remedies"] == EXCLUDED
        assert job["parameters"]["kpoint_grid"] == phonons["kpoint_grid"]
        assert phonons["kpoint_grid"] in record["kpoints"]
        assert len(record["forces"]) == phonons["supercell_atoms"]
    if case == "relax":
        _, relax = _job(workspace, phonons["relax"]["job_id"])
        tight = [run["tags"] for run in _runs(relax) if run["tags"].get("IBRION") == "1"]
        assert tight and all((t["EDIFF"], t["EDIFFG"], t["LREAL"]) == ("1E-7", "-0.001", ".FALSE.") for t in tight)
    assert len(list_jobs(workspace)) == 2 + (case == "relax")
    assert (payload / "run" / "FORCE_SETS").is_file() and (payload / "run" / "phonopy_disp.yaml").is_file()

    items = {item.record.job_id: item for item in collect(workspace, fail_fast=True, allow_job_collector=True)}
    item = items[job_id]
    expected = {
        "phonon_cell",
        "phonon_band_structure",
        "phonon_density_of_states",
        "vibrational_thermodynamics",
        "zero_point_energy",
        "phonons_result",
        "phonopy_disp",
        "force_sets",
    }
    assert set(item.outputs) == expected | ({"relaxed_structure"} if case == "relax" else set())
    assert item.outputs["zero_point_energy"].value == pytest.approx(thermodynamics["zero_point_energy"])
    # Typed records, which a store serves, not generic ones.
    assert type(item.outputs["phonon_band_structure"]).__name__ == "PhononBandStructureRecord"
    assert type(item.outputs["phonon_density_of_states"]).__name__ == "PhononDensityOfStatesRecord"
    assert item.outputs["phonons_result"].name == "phonons.json"


def _cartesian(point: list[object], lattice: list[list[float]]) -> numpy.ndarray:
    """A k-point's Cartesian coordinates (1/A, without 2 pi) from its fractional ones in the reciprocal basis."""

    return numpy.array([float(Fraction(str(x))) for x in point]) @ numpy.linalg.inv(numpy.array(lattice)).T


# standardize = true: the input is standardized, and the dielectric document is of that cell. standardize = false:
# the input is the cell the standardize path writes, which band_path keeps (it is idempotent), so both are its basis.
@pytest.mark.parametrize("standardize", (True, False))
def test_born_charges_from_dielectric_json_split_lo_and_to_at_gamma(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, standardize: bool
) -> None:
    import httk.atomistic
    import httk.core

    hpkot = _standardized(_SIC, tmp_path / "hpkot")
    workspace, job_id = _campaign(
        tmp_path / "nac",
        monkeypatch,
        {"relax": False, "standardize": standardize, "min_length": 6.0},
        poscar=_SIC if standardize else hpkot,
        dielectric=_sic_dielectric(hpkot),
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id)
    phonons = _phonons(payload)
    assert phonons["nac"]["used"] and phonons["nac"]["factor"] == pytest.approx(14.399652, abs=1e-6)
    bands = phonons["band_structure"]
    assert bands["qpoint_labels"][:2] == ["GAMMA", "X"]
    at_gamma = [mode[0] for mode in bands["frequencies"]]
    # Without the correction: three flat Si and three flat C branches. With it, the LO mode at Gamma rises.
    assert at_gamma[:2] == pytest.approx([_einstein("Si")] * 2, abs=1e-6)
    assert at_gamma[-1] > _einstein("C") + 0.5
    assert (payload / "run" / "unitcell" / "dielectric.json").is_file()

    # The path's X, in the unit cell's reciprocal basis, is band_path's X of its own cell, in Cartesian coordinates.
    band_path = json.loads((payload / "run" / "unitcell" / "band_path.json").read_text(encoding="utf-8"))
    assert band_path["basis_change"] == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    # band_path of the cell the job gave it: the input, standardized or not.
    given = payload / "run" / "relaxed" / "POSCAR"
    (tmp_path / "given.vasp").write_bytes(given.read_bytes())
    own = httk.atomistic.band_path(
        httk.core.load(str(tmp_path / "given.vasp"), format="vasp-poscar", precision=1e-3), tolerance=1e-3
    )
    own_lattice = own.structure.cell.basis.to_floats()
    unit_lattice = phonons["unit_cell"]["lattice"]
    expected = _cartesian(own.points["X"].to_fractions(), own_lattice)
    assert _cartesian(band_path["points"]["X"], unit_lattice) == pytest.approx(expected, abs=1e-9)
    assert _cartesian(bands["qpoints"][50], unit_lattice) == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("mismatch", ("lattice", "born_charges"))
def test_a_mismatched_dielectric_json_fails_before_any_child(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mismatch: str
) -> None:
    hpkot = _standardized(_SIC, tmp_path / "hpkot")
    if mismatch == "lattice":
        dielectric = _sic_dielectric(hpkot, lattice_scale=1.01)
    else:
        dielectric = _sic_dielectric(hpkot, born=None)
    workspace, job_id = _campaign(
        tmp_path / "mismatch", monkeypatch, {"relax": False, "min_length": 6.0}, poscar=_SIC, dielectric=dielectric
    )
    assert _failure(workspace, job_id)["code"] == "vasp.dielectric_mismatch"
    assert len(list_jobs(workspace)) == 1


@pytest.mark.parametrize("imaginary", ("raise", "omit"))
def test_imaginary_modes_omit_the_thermodynamics_or_their_modes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, imaginary: str
) -> None:
    from httk.workflow import collect

    # A negative curvature on the C sites: three imaginary branches at -sqrt(|k| / m_C).
    workspace, job_id = _campaign(
        tmp_path / imaginary,
        monkeypatch,
        {"relax": False, "min_length": 6.0, "imaginary": imaginary},
        poscar=_SIC,
        k={"Si": K, "C": -K},
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id)
    phonons = _phonons(payload)
    assert phonons["imaginary_modes"]["count"] > 0
    assert phonons["imaginary_modes"]["lowest_frequency"] == pytest.approx(-_einstein("C"), abs=1e-6)
    item = {i.record.job_id: i for i in collect(workspace, fail_fast=True, allow_job_collector=True)}[job_id]
    if imaginary == "raise":
        assert phonons["thermodynamics"] is None and "negative" in phonons["thermodynamics_omitted"]
        assert "zero_point_energy" not in item.outputs and "vibrational_thermodynamics" not in item.outputs
    else:
        # Only the three real Si branches remain.
        zpe = phonons["thermodynamics"]["zero_point_energy"]
        assert zpe == pytest.approx(3 * H_EV_PER_THZ * _einstein("Si") / 2, abs=1e-6)
        assert item.outputs["zero_point_energy"].value == pytest.approx(zpe)
    assert "phonon_band_structure" in item.outputs


@pytest.mark.parametrize(
    ("parameters", "poscar", "code"),
    (
        ({"max_atoms": 15, "min_length": 6.0}, _SI, "vasp.supercell_too_large"),
        ({"incar_tags": {"MAGMOM": "2*1"}}, _SI, "vasp.cell_dependent_tag_unsupported"),
        ({"incar_tags": {"NSW": "5"}}, _SI, "vasp.unsupported_incar_tag"),
        ({"imaginary": "ignore"}, _SI, "vasp.input_invalid"),
        # The conventional cell spans a sublattice of the HPKOT primitive lattice: no unimodular basis change.
        ({"relax": False, "standardize": False}, _CONVENTIONAL, "vasp.input_invalid"),
    ),
)
def test_refused_inputs_fail_before_any_call(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, parameters: dict[str, object], poscar: str, code: str
) -> None:
    workspace, job_id = _campaign(tmp_path / "refused", monkeypatch, {"relax": False, **parameters}, poscar=poscar)
    failure = _failure(workspace, job_id)
    assert failure["code"] == code
    if code == "vasp.supercell_too_large":
        assert failure["details"] == {"multiples": [2, 2, 2], "atoms": 16}
    assert len(list_jobs(workspace)) == 1


def test_a_null_magmom_is_vasps_default_and_not_a_cell_dependent_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    parameters = {"relax": False, "min_length": 6.0, "incar_tags": {"MAGMOM": None}}
    workspace, job_id = _campaign(tmp_path / "null", monkeypatch, parameters)

    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(workspace, job_id)
    children = [_job(workspace, record["job_id"])[1] for record in _phonons(payload)["children"]]
    assert children and all("MAGMOM" not in run["tags"] for child in children for run in _runs(child))
