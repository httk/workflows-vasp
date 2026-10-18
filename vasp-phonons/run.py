#!/usr/bin/env python3
"""Harmonic phonons of a structure with phonopy, by calling ``vasp.relax`` and ``vasp.static``.

This runner runs no VASP itself. It composes child jobs through httk's workflow
calls (``[workflow.calls]``, ``Attempt.call`` and ``Attempt.gather``):

* ``start`` validates the parameters and, with ``relax`` (default true), calls
  ``vasp.relax`` (child ``relax``) with this package's tight
  ``templates/INCAR.relax`` for its relax rounds; with ``relax = false`` the
  input structure is taken as relaxed.
* ``displace`` copies the relaxed structure to ``relaxed/POSCAR`` and, with
  ``standardize`` (default true), standardizes it to its HPKOT primitive cell
  with ``httk.atomistic.band_path``, as ``vasp.band-structure`` does; that cell,
  ``unitcell/POSCAR``, is the cell phonopy diagonalizes, and the q-point path of
  ``unitcell/band_path.json`` refers to its reciprocal basis. phonopy builds the
  diagonal supercell ``n_i = ceil(min_length / w_i)`` (``w_i`` the cell's width
  perpendicular to the other two vectors) and its displacements, written to
  ``phonopy_disp.yaml`` and ``displacements/<NNN>/POSCAR``; ``vasp.static`` runs
  once per displacement (``templates/INCAR.phonon``), all on one k-point mesh of
  the undisplaced supercell and with the numerics-changing remedies excluded.
* ``analyse`` records every child's INCAR, KPOINTS, remedies and forces in
  ``phonons.json``, refuses children that ran with different settings, and lets
  phonopy produce the force constants (``FORCE_SETS``), the non-analytical
  correction when the payload has ``files/dielectric.json``, the mesh, the
  density of states and the band structure along the HPKOT path; httk-analyse's
  ``harmonic_from_phonopy`` gives the thermodynamics.

``phonons.json`` holds every input of the analysis and its results, and
``collect.py`` reads only it, ``unitcell/POSCAR`` and ``relaxed/POSCAR``. The job
parameters are documented in this repository's README.
"""

import json
import math
import shutil
from fractions import Fraction
from pathlib import Path

import httk.atomistic
import httk.core
from httk.codes.vasp import (
    KPOINT_CENTERINGS,
    automatic_kpoint_grid,
    job_remedy_history_path,
    normalize_poscar_handedness,
    read_incar,
    read_poscar_header,
)
from httk.workflow import Attempt, Runner

TEMPLATES = Path(__file__).resolve().parent / "templates"
RESULT = "phonons.json"
RELAXED = Path("relaxed") / "POSCAR"
UNITCELL = Path("unitcell") / "POSCAR"
BAND_PATH = Path("unitcell") / "band_path.json"
DIELECTRIC = Path("unitcell") / "dielectric.json"
SUPERCELL = Path("supercell") / "POSCAR"
# What publish_data puts into data/, besides nothing else: what collect.py reads and phonopy's own inputs.
PUBLISHED = (RESULT, RELAXED.as_posix(), UNITCELL.as_posix(), "phonopy_disp.yaml", "FORCE_SETS")
# Forwarded to every child when this job sets them; nothing else of this job's parameters is.
SHARED = (
    "incar_tags",
    "accuracy_per_atom",
    "maximum_remedies",
    "pseudopotential_library",
    "parallel_tag",
    "parallel_value",
    "timeout",
    "vasp_command",
    "remedy_policy",
)
# The remedies that would change the numerics of one child but not the others.
EXCLUDED_REMEDIES = "bump_kpoints equal_kpoints centering incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice"
# Remedies may change these between children: they do not change the converged energy surface.
ELECTRONIC_TAGS = frozenset(
    {
        *("ALGO", "TIME", "AMIX", "BMIX", "AMIN", "ICHARG", "ISTART", "NELM", "NELMIN"),
        *("NBANDS", "NPAR", "NCORE", "KPAR", "ISYM", "SYMPREC"),
        # Ladder rungs neither excluded nor geometric: scale_ediff tightens EDIFF, dentet changes NEDOS.
        *("EDIFF", "NEDOS"),
    }
)
# Per-site, cell-total or cell-sized tags whose meaning changes with the standardized cell and the supercells.
CELL_DEPENDENT_TAGS = ("MAGMOM", "NELECT", "NUPDOWN", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF")
# The tags that make the displacement runs single points; incar_tags may not change them.
DEFINING_TAGS = ("NSW", "IBRION")
# The dielectric document must describe the unit cell within this (Angstrom, fractional).
CELL_TOLERANCE = 1e-5
DEFAULTS: dict[str, object] = {
    "relax": True,
    "standardize": True,
    "min_length": 20.0,
    "max_atoms": 300,
    "displacement": 0.01,
    "phonon_kpoint_density": 20.0,
    "centering": "Gamma",
    "mesh_length": 100.0,
    "band_points": 51,
    "temperatures": [float(t) for t in range(0, 1001, 10)],
    "cutoff_frequency": 0.05,
    "imaginary": "raise",
    "imaginary_tolerance": 0.05,
    "symmetry_tolerance": 1e-3,
}

run = Runner("vasp.phonons")


# -- parameters -------------------------------------------------------------------------------------------------------


def flag(value: object, what: str) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line."""

    text = str(value).strip().lower()
    if text not in ("true", "1", "yes", "false", "0", "no"):
        raise ValueError(f"{what} must be true or false, not {value!r}")
    return text in ("true", "1", "yes")


def real(value: object, what: str, *, positive: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{what} must be a number, not {value!r}")  # noqa: TRY004 (a job parameter error)
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{what} must be a number, not {value!r}") from None
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        raise ValueError(f"{what} must be a finite {'positive' if positive else 'nonnegative'} number, not {value!r}")
    return number


def integer(value: object, what: str, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, str)) or not str(value).strip().isdigit():
        raise ValueError(f"{what} must be an integer of at least {minimum}, not {value!r}")
    if int(value) < minimum:
        raise ValueError(f"{what} must be an integer of at least {minimum}, not {value!r}")
    return int(value)


def settings(a: Attempt) -> dict:
    """Return this job's own parameters, validated, with their defaults.

    :raises ValueError: Naming the first invalid parameter.
    """

    def get(name: str) -> object:
        value = a.parameter(name, DEFAULTS[name])
        return DEFAULTS[name] if value is None else value

    temperatures = get("temperatures")
    words = temperatures.split() if isinstance(temperatures, str) else temperatures
    if not isinstance(words, (list, tuple)) or not words:
        raise ValueError(f"temperatures must be a nonempty list of numbers, not {temperatures!r}")
    imaginary = get("imaginary")
    if imaginary not in ("raise", "omit"):
        raise ValueError(f"imaginary must be raise or omit, not {imaginary!r}")
    centering = get("centering")
    if centering not in KPOINT_CENTERINGS:
        raise ValueError(f"centering must be one of {', '.join(KPOINT_CENTERINGS)}, not {centering!r}")
    return {
        "relax": flag(get("relax"), "relax"),
        "standardize": flag(get("standardize"), "standardize"),
        "min_length": real(get("min_length"), "min_length"),
        "max_atoms": integer(get("max_atoms"), "max_atoms", 1),
        "displacement": real(get("displacement"), "displacement"),
        "phonon_kpoint_density": real(get("phonon_kpoint_density"), "phonon_kpoint_density"),
        "centering": centering,
        "mesh_length": real(get("mesh_length"), "mesh_length"),
        "band_points": integer(get("band_points"), "band_points", 2),
        "temperatures": [real(t, "a temperature", positive=False) for t in words],
        "cutoff_frequency": real(get("cutoff_frequency"), "cutoff_frequency", positive=False),
        "imaginary": imaginary,
        "imaginary_tolerance": real(get("imaginary_tolerance"), "imaginary_tolerance", positive=False),
        "symmetry_tolerance": real(get("symmetry_tolerance"), "symmetry_tolerance"),
    }


def shared_parameters(a: Attempt) -> dict:
    parameters = dict(a.parameters)
    return {name: parameters[name] for name in SHARED if parameters.get(name) is not None}


def write_json(path: Path, document: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(document, indent=1), encoding="utf-8")
    temporary.replace(path)


def read_result(a: Attempt) -> dict:
    return json.loads((a.workdir / RESULT).read_text(encoding="utf-8"))


# -- cells ------------------------------------------------------------------------------------------------------------


def cell_of(poscar: Path) -> tuple[list[list[float]], list[list[float]]]:
    """Return the lattice rows (Angstrom) and fractional site positions of a direct-coordinate POSCAR."""

    from phonopy.interface.vasp import read_vasp

    cell = read_vasp(poscar)
    return cell.cell.tolist(), cell.scaled_positions.tolist()


def dielectric_problem(document: dict, header, lattice: list, positions: list) -> str | None:
    """Return why a ``dielectric.json`` does not describe the unit cell, or ``None`` when it does."""

    if document.get("born_charges") is None or document.get("electronic") is None:
        return "it has no Born effective charges or no electronic dielectric tensor"
    if list(document.get("species") or []) != list(header.species) or list(document.get("counts") or []) != list(
        header.counts
    ):
        return f"its species {document.get('species')} {document.get('counts')} are not the unit cell's"
    if len(document["born_charges"]) != len(positions):
        return f"it has {len(document['born_charges'])} Born charges for {len(positions)} sites"
    given_lattice, given_positions = document.get("lattice") or [], document.get("positions") or []
    if len(given_lattice) != 3 or len(given_positions) != len(positions):
        return "its lattice or site positions do not match the unit cell's in shape"
    if any(abs(x - y) > CELL_TOLERANCE for row, ref in zip(given_lattice, lattice) for x, y in zip(row, ref)):
        return "its lattice differs from the unit cell's"
    for site, ref in zip(given_positions, positions):
        if any(abs((x - y) - round(x - y)) > CELL_TOLERANCE for x, y in zip(site, ref)):
            return "its site positions differ from the unit cell's"
    return None


# -- children ---------------------------------------------------------------------------------------------------------


def child_record(child) -> dict:
    """What phonons.json keeps of one displacement child: identity, outcome, settings and forces."""

    record: dict = {
        "label": child.label,
        "job_key": child.job_key,
        "job_id": child.job_id,
        "kind": child.kind,
        "failure": None if child.failure is None else {"code": child.failure.code, "message": child.failure.message},
        "incar": None,
        "kpoints": None,
        "remedies": [],
        "forces": None,
    }
    if child.payload is not None:
        history = job_remedy_history_path(child.payload)
        if history.is_file():
            events = json.loads(history.read_text(encoding="utf-8")).get("events", [])
            record["remedies"] = [change["operation"] for event in events for change in event.get("changes", [])]
    workdir = child.workdir
    if workdir is None:
        return record
    if (workdir / "INCAR").is_file():
        record["incar"] = read_incar(workdir / "INCAR")
    if (workdir / "KPOINTS").is_file():
        record["kpoints"] = (workdir / "KPOINTS").read_text(encoding="utf-8")
    if child.succeeded and (workdir / "OUTCAR").is_file():
        from httk.atomistic.integrations.vasp.io import OutcarFile

        frame = OutcarFile(workdir / "OUTCAR").last_frame
        forces = None if frame is None else frame.forces_floats()
        record["forces"] = None if forces is None else [list(row) for row in forces]
    return record


def inconsistencies(records: list[dict]) -> dict:
    """Return how the children's settings differ: KPOINTS texts, or INCAR tags other than the electronic-path ones."""

    differences: dict = {}
    kpoints = {record["label"]: record["kpoints"] for record in records}
    if len(set(kpoints.values())) > 1:
        differences["kpoints"] = kpoints
    tags = {tag for record in records for tag in record["incar"] or {}} - ELECTRONIC_TAGS
    differing = {}
    for tag in sorted(tags):
        values = {record["label"]: str((record["incar"] or {}).get(tag, "")).strip().lower() for record in records}
        if len(set(values.values())) > 1:
            differing[tag] = values
    if differing:
        differences["incar_tags"] = differing
    return differences


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def start(a: Attempt) -> None:
    """Validate the parameters, then call the tight relaxation, or with ``relax = false`` go on to the displacements."""

    try:
        own = settings(a)
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail("vasp.input_missing", f"the starting structure {poscar.name} is not in this payload")
        return
    explicit = {str(tag).strip().upper() for tag in dict(a.parameter("incar_tags", {}) or {})}
    for tag in DEFINING_TAGS:
        if tag in explicit:
            a.fail(
                "vasp.unsupported_incar_tag",
                f"{tag} in incar_tags is not supported: the displacement runs are single points",
                details={"tag": tag},
            )
            return
    # Checked before the relaxation: the relaxed CONTCAR keeps the input POSCAR's comment.
    comment = {"MAGMOM"} if "[MAGMOM=" in read_poscar_header(poscar).comment else set()
    for source, tags in (("incar_tags", explicit), (f"the comment of {poscar.name}", comment)):
        for tag in CELL_DEPENDENT_TAGS:
            if tag in tags:
                a.fail(
                    "vasp.cell_dependent_tag_unsupported",
                    f"{tag} from {source} is not supported: the phonon supercells' sites and size differ from the "
                    "input structure's (magnetic phonons are out of scope)",
                    details={"tag": tag, "source": source},
                )
                return
    if not own["relax"]:
        (a.workdir / RELAXED).parent.mkdir(exist_ok=True)
        shutil.copyfile(poscar, a.workdir / RELAXED)
        a.log.append("note", "relax = false: the input structure is taken as relaxed")
        a.advance("displace")
        return
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    files = {"POSCAR": poscar, "INCAR.relax": TEMPLATES / "INCAR.relax"}
    override = a.payload / "files" / "INCAR.relax"
    if override.is_file():
        files["INCAR.relax"] = override
    if potcar.is_file():
        files["POTCAR"] = potcar
    a.call("relax", label="relax", files=files, parameters=shared_parameters(a))
    a.log.append(
        "note", f"called the tight relaxation with {files['INCAR.relax'].name} from {files['INCAR.relax'].parent}"
    )
    a.gather("displace", when="all_succeeded", on_impossible="relax_failed")


@run.step
def relax_failed(a: Attempt) -> None:
    """Fail by name with the relaxation child's failure."""

    child = a.children.get("relax")
    failure = None if child is None or child.failure is None else child.failure
    a.fail(
        "vasp.relax_failed",
        "the tight relaxation did not succeed",
        details={
            "job_key": None if child is None else child.job_key,
            "job_id": None if child is None else child.job_id,
            "failure": None if failure is None else {"code": failure.code, "message": failure.message},
        },
    )


@run.step
def displace(a: Attempt) -> None:
    """Fix the phonon unit cell and its q-point path, and call one single point per displaced supercell."""

    own = settings(a)
    relaxed = a.workdir / RELAXED
    child = a.children.get("relax")
    if child is not None:
        if child.workdir is None:
            a.fail("vasp.child_workdir_missing", f"the relaxation {child.job_key} has no workdir")
            return
        if not (child.workdir / "CONTCAR").is_file():
            a.fail("vasp.relax_failed", f"the relaxation {child.job_key} left no CONTCAR")
            return
        relaxed.parent.mkdir(exist_ok=True)
        shutil.copyfile(child.workdir / "CONTCAR", relaxed)
    header = read_poscar_header(relaxed)
    unitcell = a.workdir / UNITCELL
    unitcell.parent.mkdir(exist_ok=True)
    tolerance = own["symmetry_tolerance"]
    try:
        # From the file as given: the reader takes the digits written as its precision, and the handedness
        # normalization below rewrites them. The unit cell is band_path's cell, or the given one.
        structure = httk.core.load(str(relaxed), format="vasp-poscar", precision=tolerance)
        result = httk.atomistic.band_path(structure, tolerance=tolerance)
        if own["standardize"]:
            httk.core.save(result.structure, str(unitcell))
        else:
            shutil.copyfile(relaxed, unitcell)
        # Once, on the unit cell, so every displaced supercell inherits it; the basis change below absorbs it.
        normalize_poscar_handedness(unitcell)
    except ImportError as exception:
        a.fail("vasp.band_path_failed", f"cannot recognize the symmetry of the structure: {exception}")
        return
    except ValueError as exception:
        a.fail("vasp.band_path_failed", f"cannot find the band path of the structure: {exception}")
        return
    # A supplied POTCAR is per species in the input's order.
    unit_header = read_poscar_header(unitcell)
    if unit_header.species != header.species:
        a.fail(
            "vasp.input_mismatch",
            f"the HPKOT cell lists its species as {' '.join(unit_header.species)}, not in the input's order",
            details={"input": list(header.species), "written": list(unit_header.species)},
        )
        return
    try:
        from phonopy import Phonopy
        from phonopy.interface.vasp import read_vasp, write_vasp
    except ImportError as exception:
        a.fail("vasp.analysis_unavailable", f"phonons need phonopy: {exception} (pip install phonopy)")
        return
    lattice, positions = cell_of(unitcell)
    path_lattice = [[float(x) for x in row] for row in result.structure.cell.basis.to_floats()]
    # band_path's transform starts from httk's primitive cell, not the unit cell, so the cells are compared: with
    # rows as vectors, L_path = T L_unit, and the unit cell is accepted when T is a unimodular integer matrix (the
    # same lattice in the same frame; band_path need not return the very basis it was given, e.g. for zincblende).
    # Standardized, T is the identity, or -I when the handedness was normalized.
    basis_change = basis_change_of(path_lattice, lattice)
    if basis_change is None:
        a.fail(
            "vasp.input_invalid" if not own["standardize"] else "vasp.band_path_failed",
            "the given unit cell is not a basis of its HPKOT lattice, so the q-point path would not refer to it: "
            "pass standardize = true",
            details={"unit_cell": lattice, "hpkot_cell": path_lattice},
        )
        return
    symbol = result.bravais_lattice_extended
    write_json(
        a.workdir / BAND_PATH,
        {
            "convention": "HPKOT",
            "bravais_lattice": result.bravais_lattice,
            "bravais_lattice_extended": symbol,
            "spacegroup": result.spacegroup.it_number,
            "symmetry_tolerance": tolerance,
            # Fractional coordinates of the unit cell's reciprocal basis, transformed from band_path's cell.
            "points": {
                label: [str(v) for v in to_unit_cell(point.to_fractions(), basis_change)]
                for label, point in result.points.items()
            },
            "path": [list(segment) for segment in result.path],
            "transform": [[int(v) for v in row] for row in result.transform.to_fractions()],
            "basis_change": basis_change,
        },
    )
    given = a.payload / "files" / "dielectric.json"
    if given.is_file():
        problem = dielectric_problem(json.loads(given.read_text(encoding="utf-8")), unit_header, lattice, positions)
        if problem is not None:
            a.fail(
                "vasp.dielectric_mismatch",
                f"files/dielectric.json does not describe the phonon unit cell: {problem}",
                details={"unit_cell": {"species": list(unit_header.species), "lattice": lattice}},
            )
            return
        shutil.copyfile(given, a.workdir / DIELECTRIC)

    # The diagonal supercell by perpendicular widths w_i = V / |a_j x a_k|.
    a1, a2, a3 = lattice
    volume = abs(_dot(a1, _cross(a2, a3)))
    widths = [volume / math.sqrt(_dot(c, c)) for c in (_cross(a2, a3), _cross(a3, a1), _cross(a1, a2))]
    # ponytail: diagonal supercells only; atomate2's cubic-like (non-diagonal) search is the upgrade.
    multiples = [math.ceil(own["min_length"] / w - 1e-9) for w in widths]
    atoms = len(positions) * multiples[0] * multiples[1] * multiples[2]
    if atoms > own["max_atoms"]:
        a.fail(
            "vasp.supercell_too_large",
            f"the {'x'.join(map(str, multiples))} supercell for min_length {own['min_length']} has {atoms} atoms, "
            f"more than max_atoms = {own['max_atoms']}",
            details={"multiples": multiples, "atoms": atoms},
        )
        return
    phonon = Phonopy(
        read_vasp(unitcell),
        supercell_matrix=[[n if i == j else 0 for j in range(3)] for i, n in enumerate(multiples)],
        primitive_matrix="P",
    )
    phonon.generate_displacements(distance=own["displacement"])
    phonon.save(a.workdir / "phonopy_disp.yaml")
    supercell = a.workdir / SUPERCELL
    supercell.parent.mkdir(exist_ok=True)
    write_vasp(supercell, phonon.supercell)
    cells = {f"p{index:03d}": cell for index, cell in enumerate(phonon.supercells_with_displacements)}
    written = {"supercell": supercell}
    for index, (label, cell) in enumerate(cells.items()):
        path = a.workdir / "displacements" / f"{index:03d}" / "POSCAR"
        path.parent.mkdir(parents=True, exist_ok=True)
        write_vasp(path, cell)
        written[label] = path
    for label, path in written.items():
        # The forces come back in POSCAR order and the POTCAR is per species: both need the unit cell's species line.
        if read_poscar_header(path).species != unit_header.species:
            a.fail(
                "vasp.input_mismatch",
                f"the {label} POSCAR lists its species as {' '.join(read_poscar_header(path).species)}, "
                f"not as the unit cell's {' '.join(unit_header.species)}",
            )
            return
    grid = automatic_kpoint_grid(own["phonon_kpoint_density"], poscar=supercell)
    mesh = " ".join(map(str, grid))
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    parameters = shared_parameters(a) | {
        "kpoint_grid": mesh,
        "centering": own["centering"],
        "excluded_remedies": EXCLUDED_REMEDIES,
    }
    for label in cells:
        files = {"POSCAR": written[label], "INCAR.static": TEMPLATES / "INCAR.phonon"}
        if potcar.is_file():
            files["POTCAR"] = potcar
        a.call("static", label=label, files=files, parameters=parameters)
    first_atoms = phonon.dataset["first_atoms"]
    write_json(
        a.workdir / RESULT,
        {
            "workflow": "vasp.phonons",
            "parameters": own,
            "relax": None if child is None else {"job_key": child.job_key, "job_id": child.job_id},
            "unit_cell": {
                "species": list(unit_header.species),
                "counts": list(unit_header.counts),
                "lattice": lattice,
                "positions": positions,
            },
            "band_path": {"bravais_lattice_extended": symbol, "spacegroup": result.spacegroup.it_number},
            "supercell_matrix": [[n if i == j else 0 for j in range(3)] for i, n in enumerate(multiples)],
            "supercell_atoms": atoms,
            "displacements": [
                {
                    "label": label,
                    "atom": int(entry["number"]),
                    "displacement": [float(x) for x in entry["displacement"]],
                }
                for label, entry in zip(cells, first_atoms, strict=True)
            ],
            "kpoint_grid": mesh,
            "centering": own["centering"],
            "excluded_remedies": EXCLUDED_REMEDIES,
            "nac": (a.workdir / DIELECTRIC).is_file(),
        },
    )
    a.state.merge({"supercell": "x".join(map(str, multiples)), "displacements": len(cells)})
    a.log.append(
        "note",
        f"standardize = {own['standardize']}: phonons of the {symbol} cell in a {'x'.join(map(str, multiples))} "
        f"supercell of {atoms} atoms; called {len(cells)} vasp.static displacements on the {mesh} {own['centering']} mesh"
        + ("; with the non-analytical correction from files/dielectric.json" if given.is_file() else ""),
    )
    a.gather("analyse", when="all_terminal")


def basis_change_of(path_lattice: list[list[float]], lattice: list[list[float]]) -> list[list[int]] | None:
    """Return the integer ``T`` with ``path_lattice = T lattice`` (rows as vectors), or ``None`` unless it is unimodular."""

    import numpy

    change = numpy.array(path_lattice) @ numpy.linalg.inv(numpy.array(lattice))
    rounded = numpy.rint(change)
    if numpy.abs(change - rounded).max() > 1e-6 or round(abs(numpy.linalg.det(rounded))) != 1:
        return None
    return [[int(x) for x in row] for row in rounded]


def to_unit_cell(point: list[Fraction], basis_change: list[list[int]]) -> list[Fraction]:
    """Return a k-point of the path cell's reciprocal basis in the unit cell's: ``q T^-T`` for ``L_path = T L_unit``.

    Reciprocal rows are ``B = L^-T``, so ``B_path = T^-T B_unit`` and ``k = q B_path = (q T^-T) B_unit``.
    """

    import numpy

    # T is unimodular, so its inverse is an integer matrix too.
    inverse_transposed = numpy.rint(numpy.linalg.inv(numpy.array(basis_change, dtype=float))).astype(int).T
    return [sum((point[i] * int(inverse_transposed[i][j]) for i in range(3)), Fraction(0)) for j in range(3)]


def _dot(u: list[float], v: list[float]) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _cross(u: list[float], v: list[float]) -> list[float]:
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


@run.step
def analyse(a: Attempt) -> None:
    """Record every displacement child, check that they agree, and compute the phonons from their forces."""

    own = settings(a)
    result = read_result(a)
    records = [child_record(child) for child in a.children]
    result["children"] = records
    write_json(a.workdir / RESULT, result)
    failed = [record for record in records if record["kind"] != "succeeded"]
    if failed:
        a.fail(
            "vasp.dependency_failed",
            f"displacements {', '.join(record['label'] for record in failed)} did not succeed",
            details={"failed": {record["label"]: record["failure"] for record in failed}},
        )
        return
    missing = [child.label for child in a.children if child.workdir is None]
    if missing:
        a.fail("vasp.child_workdir_missing", f"displacements {', '.join(missing)} have no workdir")
        return
    differences = inconsistencies(records)
    if differences:
        a.fail(
            "vasp.inconsistent_children",
            "the displacements ran with different settings: "
            + ", ".join((["KPOINTS"] if "kpoints" in differences else []) + list(differences.get("incar_tags", {}))),
            details=differences,
        )
        return
    atoms = result["supercell_atoms"]
    bad = [record["label"] for record in records if record["forces"] is None or len(record["forces"]) != atoms]
    if bad:
        a.fail(
            "vasp.forces_missing",
            f"displacements {', '.join(bad)} have no forces on all {atoms} supercell atoms in their OUTCAR",
        )
        return
    try:
        import numpy
        import phonopy
        from httk.analyse.integrations.phonopy import harmonic_from_phonopy
        from phonopy.file_IO import write_FORCE_SETS
        from phonopy.physical_units import get_physical_units
    except ImportError as exception:
        a.fail(
            "vasp.analysis_unavailable",
            f"the phonon analysis needs phonopy and httk-analyse: {exception} (pip install phonopy httk-analyse)",
        )
        return
    phonon = phonopy.load(a.workdir / "phonopy_disp.yaml", produce_fc=False, is_nac=False)
    labels = [entry["label"] for entry in result["displacements"]]
    if len(phonon.dataset["first_atoms"]) != len(labels):
        a.fail("vasp.inconsistent_children", "phonopy_disp.yaml does not list the displacements that were called")
        return
    by_label = {record["label"]: record["forces"] for record in records}
    phonon.forces = numpy.array([by_label[label] for label in labels])
    phonon.produce_force_constants()
    write_FORCE_SETS(phonon.dataset, a.workdir / "FORCE_SETS")
    if (a.workdir / DIELECTRIC).is_file():
        dielectric = json.loads((a.workdir / DIELECTRIC).read_text(encoding="utf-8"))
        units = get_physical_units()
        # Born rows are the field direction, as VASP prints them and phonopy reads them; the factor is e^2/(4 pi eps0)
        # in eV Angstrom, which phonopy needs explicitly for VASP units.
        phonon.nac_params = {
            "born": numpy.array(dielectric["born_charges"]),
            "dielectric": numpy.array(dielectric["electronic"]),
            "factor": units.Hartree * units.Bohr,
        }
        result["nac"] = {"used": True, "factor": units.Hartree * units.Bohr}
    else:
        result["nac"] = {"used": False}
    phonon.run_mesh(own["mesh_length"])
    frequencies = numpy.asarray(phonon.mesh.frequencies)
    imaginary = frequencies < -own["imaginary_tolerance"]
    result["mesh"] = {
        "length": own["mesh_length"],
        "numbers": [int(n) for n in phonon.mesh.mesh_numbers],
        "irreducible_qpoints": len(frequencies),
    }
    result["imaginary_modes"] = {
        "tolerance": own["imaginary_tolerance"],
        "count": int(imaginary.sum()),
        "lowest_frequency": float(frequencies.min()),
    }
    try:
        thermodynamics = harmonic_from_phonopy(
            phonon,
            own["temperatures"],
            zero_modes="omit",
            imaginary=own["imaginary"],
            cutoff_frequency=own["cutoff_frequency"],
        )
    except ValueError as exception:
        result["thermodynamics"] = None
        result["thermodynamics_omitted"] = str(exception)
        a.log.append("note", f"no thermodynamics (imaginary = {own['imaginary']}): {exception}")
    else:
        result["thermodynamics"] = {name: getattr(thermodynamics, name) for name in thermodynamics.__dataclass_fields__}
    if imaginary.any():
        a.log.append(
            "note",
            f"{int(imaginary.sum())} imaginary modes on the mesh below -{own['imaginary_tolerance']} THz, "
            f"the lowest {frequencies.min():.4f} THz",
        )
    # phonopy sizes the DOS grid by the spectrum's width, which a dispersionless (Einstein) spectrum lacks: a
    # 0.01 THz Gaussian then gives it one.
    phonon.run_total_dos(sigma=None if numpy.ptp(frequencies) > 1e-6 else 0.01)
    result["dos"] = {
        "frequencies": [float(x) for x in phonon.total_dos.frequency_points],
        "density": [float(x) for x in phonon.total_dos.dos],
    }
    path = json.loads((a.workdir / BAND_PATH).read_text(encoding="utf-8"))
    points = {label: [float(Fraction(x)) for x in point] for label, point in path["points"].items()}
    count = own["band_points"]
    segments = [
        [[s + (e - s) * i / (count - 1) for s, e in zip(points[start], points[end])] for i in range(count)]
        for start, end in path["path"]
    ]
    phonon.run_band_structure(segments)
    band = phonon.band_structure
    stacked = numpy.concatenate([numpy.asarray(f) for f in band.frequencies])
    result["band_structure"] = {
        "qpoints": [list(map(float, q)) for segment in band.qpoints for q in segment],
        "qpoint_label_indices": [i * count + j for i in range(len(segments)) for j in (0, count - 1)],
        "qpoint_labels": [label for segment in path["path"] for label in segment],
        # Sorted at each q-point, as atomate2.phonons' collect does: frequencies[mode][q].
        "frequencies": numpy.sort(stacked, axis=1).T.tolist(),
    }
    write_json(a.workdir / RESULT, result)
    headline = {
        "imaginary_modes": result["imaginary_modes"]["count"],
        "lowest_frequency": result["imaginary_modes"]["lowest_frequency"],
        "zero_point_energy": None
        if result["thermodynamics"] is None
        else result["thermodynamics"]["zero_point_energy"],
        "nac": result["nac"]["used"],
    }
    a.state.merge(headline)
    a.log.append("note", f"phonons: {headline}")
    if a.parameter("publish_data", False) is True:
        base = a.parameter("data_prefix", "") or ""
        for name in PUBLISHED:
            a.put(a.workdir / name, f"{base}/{name}" if base else name)
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
