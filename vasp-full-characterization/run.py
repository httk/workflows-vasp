#!/usr/bin/env python3
"""Characterize one material with the other packages of this repository as child jobs.

This runner runs no VASP itself. The schedule::

    start ── call relax ──► plan ── standardize the relaxed cell, decide insulator or metal, call
             (vasp.relax,           band_structure, elastic, eos, and dielectric (insulator) or
              tight INCAR.relax)    phonons (metal) ──► progress
    progress ── whenever a child ends: once band_structure has ended, call dos on its density;
                once dielectric has ended, call phonons with its Born charges; when no child is
                live, write characterization.json and finish

No VASP calculation runs twice. The tight relaxation is the only one: elastic, eos
and phonons start from the relaxed, standardized cell with ``relax = false``; dos
reads the band-structure child's self-consistent density (``files/CHGCAR`` and
``files/scf-INCAR``), so the band-structure scf is the only one; phonons reads the
dielectric child's ``dielectric.json`` for its non-analytical correction. Every
child gets the same cell (``standardized/POSCAR``, the HPKOT primitive cell of the
relaxed structure, as ``vasp.band-structure`` builds it), so the Born charges refer
to the cell phonons runs on.

The insulator decision comes from the relaxation, as atomate2's parent decides it:
the band edges of the final relax round's EIGENVAL occupations; an insulator has an
indirect gap above ``band_gap_tolerance``. Without a readable EIGENVAL with
occupations the job fails ``vasp.band_gap_unknown``: nothing is guessed.

The job fails ``vasp.relax_failed`` when the relaxation fails. Otherwise it waits
for every child, then succeeds when all succeeded and fails
``vasp.characterization_incomplete`` naming the failed ones. ``characterization.json``
(in the workdir and the job data) records every child and the decisions taken.
Every property is collected on the child runs. The job parameters are documented
in this repository's README.
"""

import json
import math
import shutil
from collections.abc import Mapping
from pathlib import Path

from httk.codes.vasp import read_poscar_header
from httk.workflow import Attempt, ChildResult, Runner

TEMPLATES = Path(__file__).resolve().parent / "templates"
RESULT = "characterization.json"
RELAXED = Path("relaxed") / "POSCAR"
STANDARDIZED = Path("standardized") / "POSCAR"
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
PARTS = ("band_structure", "dos", "elastic", "eos", "dielectric", "phonons")
# The parameters this job sets itself on a child; they win over <label>_parameters.
STRUCTURAL = {"elastic": {"relax": False}, "eos": {"relax": False}, "phonons": {"relax": False, "standardize": False}}
# Per-site, cell-total or cell-sized tags whose meaning changes with the standardized cell (band structure and
# phonons refuse them; magnetic materials are out of scope).
CELL_DEPENDENT_TAGS = ("MAGMOM", "NELECT", "NUPDOWN", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF")
TERMINAL = ("succeeded", "failed", "cancelled")

run = Runner("vasp.full-characterization")


# -- parameters -------------------------------------------------------------------------------------------------------


def flag(value: object, what: str) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line."""

    text = str(value).strip().lower()
    if text not in ("true", "1", "yes", "false", "0", "no"):
        raise ValueError(f"{what} must be true or false, not {value!r}")
    return text in ("true", "1", "yes")


def real(value: object, what: str, *, positive: bool) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{what} must be a number, not {value!r}")  # noqa: TRY004 (a bad parameter fails vasp.input_invalid)
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{what} must be a number, not {value!r}") from None
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        raise ValueError(f"{what} must be a finite {'positive' if positive else 'nonnegative'} number, not {value!r}")
    return number


def settings(a: Attempt) -> dict:
    """Return this job's own parameters, validated, with their defaults.

    :raises ValueError: Naming the first invalid parameter.
    """

    value = a.parameter("parts", " ".join(PARTS))
    words = value.split() if isinstance(value, str) else value
    if not isinstance(words, (list, tuple)) or not words or not all(isinstance(word, str) for word in words):
        raise ValueError(f"parts must name at least one of {' '.join(PARTS)}, not {value!r}")
    unknown = [word for word in words if word not in PARTS]
    if unknown:
        raise ValueError(f"parts names unknown parts {' '.join(unknown)}: choose from {' '.join(PARTS)}")
    for label in ("relax", *PARTS):
        given = a.parameter(f"{label}_parameters", None)
        if given is not None and not isinstance(given, Mapping):
            raise ValueError(f"{label}_parameters must be a dictionary of job parameters, not {given!r}")
    return {
        "parts": [part for part in PARTS if part in words],
        "standardize": flag(a.parameter("standardize", True), "standardize"),
        "band_gap_tolerance": real(a.parameter("band_gap_tolerance", 0.1), "band_gap_tolerance", positive=False),
        "symmetry_tolerance": real(a.parameter("symmetry_tolerance", 1e-3), "symmetry_tolerance", positive=True),
    }


def child_parameters(a: Attempt, label: str) -> dict:
    """The shared parameters, then ``<label>_parameters``, then the parameters this job sets itself."""

    parameters = dict(a.parameters)
    shared = {name: parameters[name] for name in SHARED if parameters.get(name) is not None}
    return shared | dict(a.parameter(f"{label}_parameters", None) or {}) | STRUCTURAL.get(label, {})


def potcar_files(a: Attempt) -> dict[str, Path]:
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    return {"POTCAR": potcar} if potcar.is_file() else {}


def call(a: Attempt, label: str, files: dict[str, Path] | None = None) -> dict:
    """Call one child on the standardized cell (unless *files* say otherwise) and return its bookkeeping entry."""

    files = files or {"POSCAR": a.workdir / STANDARDIZED, **potcar_files(a)}
    reference = a.call(label, label=label, files=files, parameters=child_parameters(a, label))
    return {"job_id": reference.job_id, "job_key": reference.job_key, "kind": "called", "failure": None}


def part(child: ChildResult) -> dict:
    """The bookkeeping entry of one observed child."""

    failure = None if child.failure is None else {"code": child.failure.code, "message": child.failure.message}
    return {"job_id": child.job_id, "job_key": child.job_key, "kind": child.kind, "failure": failure}


# -- decisions --------------------------------------------------------------------------------------------------------


def band_gap(eigenval_path: Path) -> dict:
    """Return the band edges of an EIGENVAL with occupations, both spin channels stacked along the band axis.

    :raises OSError: If the file cannot be read.
    Occupations outside [-1e-3, 1 + 1e-3] are taken as 0.5: VASP's Blöchl-corrected tetrahedron method
    (``ISMEAR = -5``, the relaxation's) gives such values only to partially occupied states at the Fermi level.

    :raises ValueError: If it has no occupations or no band edges.
    """

    from httk.analyse.matsci.electronic import band_edges
    from httk.atomistic.integrations.vasp.io import read_eigenval

    eigenval = read_eigenval(eigenval_path)
    occupations = eigenval.occupations_floats()
    if occupations is None:
        raise ValueError("it has no occupation columns")
    energies = eigenval.energies_floats()
    spins, kpoints, bands = range(eigenval.ispin), range(eigenval.nkpoints), range(eigenval.nbands)
    stacked = [[energies[s][k][b] for k in kpoints] for s in spins for b in bands]
    held = [[occupations[s][k][b] for k in kpoints] for s in spins for b in bands]
    repaired = sum(1 for row in held for value in row if not -1e-3 <= value <= 1 + 1e-3)
    held = [[value if -1e-3 <= value <= 1 + 1e-3 else 0.5 for value in row] for row in held]
    # The gap does not depend on the reference energy.
    edges = band_edges(stacked, held, maximum_occupation=1, occupation_tolerance=1e-3, energy_reference=0.0)
    return {
        "metallic": edges.metallic,
        "indirect_gap": edges.indirect_gap,
        "direct_gap": edges.direct_gap,
        "occupations_out_of_range": repaired,
    }


def density_files(child: ChildResult) -> dict[str, Path] | None:
    """The band-structure child's HPKOT cell, POTCAR and scf density, as ``vasp.dos`` takes them; ``None`` if absent."""

    workdir = child.workdir
    if not child.succeeded or workdir is None:
        return None
    files = {
        "POSCAR": workdir / "POSCAR",
        "POTCAR": workdir / "POTCAR",
        "CHGCAR": workdir / "CHGCAR",
        "scf-INCAR": workdir / "scf" / "INCAR",
        "scf-OUTCAR": workdir / "scf" / "OUTCAR",
    }
    if not all(path.is_file() for name, path in files.items() if name not in ("POTCAR", "scf-OUTCAR")):
        return None
    return {name: path for name, path in files.items() if path.is_file()}


def write_result(a: Attempt, document: dict) -> None:
    path = a.workdir / RESULT
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8")
    temporary.replace(path)


def finish(a: Attempt, parts: dict, decisions: dict) -> None:
    """Write and publish the result, and succeed when every called child succeeded."""

    nac = decisions.get("phonons_nac")
    write_result(a, {"workflow": "vasp.full-characterization", "parts": parts, "decisions": decisions})
    # The result and both cells go to data/, where collect.py reads them.
    for name in (RESULT, RELAXED.as_posix(), STANDARDIZED.as_posix()):
        a.put(a.workdir / name, name)
    failed = {label: entry["failure"] for label, entry in parts.items() if entry["kind"] != "succeeded"}
    # A polar insulator's phonons without the correction stay incomplete, even after the dielectric is repaired.
    if nac is False and decisions["insulator"] and "dielectric" in decisions["parts"]:
        failed.pop("phonons", None)
        failed["phonons (without the non-analytical correction: dielectric did not succeed)"] = parts["phonons"][
            "failure"
        ]
    if failed:
        a.fail(
            "vasp.characterization_incomplete",
            f"these parts did not succeed: {', '.join(failed)}",
            details={"failed": failed},
        )
        return
    a.succeed()


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def start(a: Attempt) -> None:
    """Validate the parameters, then call the tight relaxation."""

    try:
        settings(a)
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail("vasp.input_missing", f"the starting structure {poscar.name} is not in this payload")
        return
    # Every child runs on the standardized cell; the relaxed CONTCAR keeps the input POSCAR's comment.
    sources = {"incar_tags": a.parameter("incar_tags", None) or {}}
    for label in ("relax", *PARTS):
        sources[f"{label}_parameters.incar_tags"] = (a.parameter(f"{label}_parameters", None) or {}).get(
            "incar_tags"
        ) or {}
    sources[f"the comment of {poscar.name}"] = (
        {"MAGMOM": ""} if "[MAGMOM=" in read_poscar_header(poscar).comment else {}
    )
    for source, tags in sources.items():
        for tag in CELL_DEPENDENT_TAGS:
            # A null tag is removed (VASP's default), which depends on no cell.
            if tag in {str(name).strip().upper() for name, value in dict(tags).items() if value is not None}:
                a.fail(
                    "vasp.cell_dependent_tag_unsupported",
                    f"{tag} from {source} is not supported: the children run on the standardized cell, whose sites "
                    "and size differ from the input structure's (magnetic materials are out of scope)",
                    details={"tag": tag, "source": source},
                )
                return
    override = a.payload / "files" / "INCAR.relax"
    template = override if override.is_file() else TEMPLATES / "INCAR.relax"
    files = {"POSCAR": poscar, "INCAR.relax": template, **potcar_files(a)}
    a.call("relax", label="relax", files=files, parameters=child_parameters(a, "relax"))
    a.log.append("note", f"called the tight relaxation with {template.name} from {template.parent}")
    a.gather("plan", when="all_succeeded", on_impossible="relax_failed")


@run.step
def relax_failed(a: Attempt) -> None:
    """Fail by name with the relaxation child's failure: nothing can be characterized without it."""

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
def plan(a: Attempt) -> None:
    """Standardize the relaxed cell, decide insulator or metal, and call every child that needs only the cell."""

    own = settings(a)
    relax = a.children["relax"]
    if relax.workdir is None:
        a.fail("vasp.child_workdir_missing", f"the relaxation {relax.job_key} has no workdir")
        return
    if not (relax.workdir / "CONTCAR").is_file():
        a.fail("vasp.relax_failed", f"the relaxation {relax.job_key} left no CONTCAR")
        return
    relaxed, standardized = a.workdir / RELAXED, a.workdir / STANDARDIZED
    relaxed.parent.mkdir(exist_ok=True)
    standardized.parent.mkdir(exist_ok=True)
    shutil.copyfile(relax.workdir / "CONTCAR", relaxed)
    tolerance = own["symmetry_tolerance"]
    if own["standardize"]:
        import httk.atomistic
        import httk.core

        try:
            # As vasp.band-structure's prepare builds its HPKOT cell.
            structure = httk.core.load(str(relaxed), format="vasp-poscar", precision=tolerance)
            httk.core.save(httk.atomistic.band_path(structure, tolerance=tolerance).structure, str(standardized))
        except ImportError as exception:
            a.fail("vasp.band_path_failed", f"cannot recognize the symmetry of the relaxed structure: {exception}")
            return
        except ValueError as exception:
            a.fail("vasp.band_path_failed", f"cannot standardize the relaxed structure: {exception}")
            return
        # A supplied POTCAR is per species in the input's order.
        species, written = read_poscar_header(relaxed).species, read_poscar_header(standardized).species
        if written != species:
            a.fail(
                "vasp.input_mismatch",
                f"the HPKOT cell lists its species as {' '.join(written)}, not in the relaxed structure's order",
                details={"input": list(species), "written": list(written)},
            )
            return
    else:
        shutil.copyfile(relaxed, standardized)
    try:
        gap = band_gap(relax.workdir / "EIGENVAL")
    except ImportError as exception:
        a.fail("vasp.analysis_unavailable", f"the band gap needs httk-analyse: {exception} (pip install httk-analyse)")
        return
    except (OSError, ValueError) as exception:
        a.fail(
            "vasp.band_gap_unknown",
            f"cannot tell an insulator from a metal: the relaxation's EIGENVAL is unusable: {exception}",
        )
        return
    insulator = not gap["metallic"] and gap["indirect_gap"] > own["band_gap_tolerance"]
    requested = own["parts"]
    notes = []
    phonons_nac = None
    parts = {"relax": part(relax)}
    for label in ("band_structure", "elastic", "eos"):
        if label in requested:
            parts[label] = call(a, label)
    if "dos" in requested and "band_structure" not in requested:
        parts["dos"] = call(a, "dos")
        notes.append("dos without band_structure: dos runs its own scf")
    if insulator and "dielectric" in requested:
        parts["dielectric"] = call(a, "dielectric")
        if "phonons" in requested:
            notes.append("insulator: phonons waits for the dielectric child's Born charges")
    elif "phonons" in requested:
        parts["phonons"] = call(a, "phonons")
        phonons_nac = False
        notes.append(
            ("metal" if not insulator else "dielectric not requested") + ": phonons without non-analytical correction"
        )
    if not insulator and "dielectric" in requested:
        notes.append("metal: no dielectric child")
    decisions = {
        "parts": requested,
        "standardize": own["standardize"],
        "symmetry_tolerance": tolerance,
        "band_gap": 0.0 if gap["metallic"] else gap["indirect_gap"],
        "metallic": gap["metallic"],
        "occupations_out_of_range": gap["occupations_out_of_range"],
        "band_gap_tolerance": own["band_gap_tolerance"],
        "band_gap_source": "relax EIGENVAL",
        "insulator": insulator,
        "phonons_nac": phonons_nac,
        "notes": notes,
    }
    for note in notes:
        a.log.append("note", note)
    a.state.merge({"parts": parts, "decisions": decisions})
    if len(parts) == 1:
        finish(a, parts, decisions)
        return
    a.gather("progress", when="any_terminal")


@run.step
def progress(a: Attempt) -> None:
    """Record the children that ended, call the children waiting for them, and finish when none is live."""

    # Every live child is in this join (each wait rejoins them), so only terminal entries are kept from state.
    parts = {label: entry for label, entry in a.state["parts"].items() if entry["kind"] in TERMINAL}
    parts.update({child.label: part(child) for child in a.children})
    decisions = a.state["decisions"]
    requested, notes = decisions["parts"], decisions["notes"]
    band_structure = a.children.get("band_structure")
    if "dos" in requested and "dos" not in parts and band_structure is not None and band_structure.kind in TERMINAL:
        given = density_files(band_structure)
        if given is not None:
            parts["dos"] = call(a, "dos", files=given)
            notes.append("dos reads the band-structure child's scf density: no second scf")
        else:
            parts["dos"] = call(a, "dos")
            notes.append(f"band_structure {band_structure.kind} without its density: dos runs its own scf")
    dielectric = a.children.get("dielectric")
    if "phonons" in requested and "phonons" not in parts and dielectric is not None and dielectric.kind in TERMINAL:
        document = None if dielectric.workdir is None else dielectric.workdir / "dielectric.json"
        if dielectric.succeeded and document is not None and document.is_file():
            parts["phonons"] = call(
                a, "phonons", {"POSCAR": a.workdir / STANDARDIZED, **potcar_files(a), "dielectric.json": document}
            )
            decisions["phonons_nac"] = True
            notes.append("phonons with the dielectric child's Born charges (non-analytical correction)")
        else:
            parts["phonons"] = call(a, "phonons")
            decisions["phonons_nac"] = False
            notes.append(f"dielectric {dielectric.kind}: phonons without non-analytical correction")
    dos = a.children.get("dos")
    if dos is not None and dos.succeeded and dos.workdir is not None and (dos.workdir / "band_edges.json").is_file():
        edges = json.loads((dos.workdir / "band_edges.json").read_text(encoding="utf-8"))
        gap = 0.0 if edges["metallic"] else edges["indirect_gap"]
        decisions["dos_band_gap"] = {"band_gap": gap, "metallic": edges["metallic"], "source": "dos band_edges.json"}
        if (not edges["metallic"] and gap > decisions["band_gap_tolerance"]) != decisions["insulator"]:
            notes.append(f"the dos band gap {gap:.4f} eV disagrees with the relaxation's insulator decision")
    a.state.merge({"parts": parts, "decisions": decisions})
    live = [label for label, entry in parts.items() if entry["kind"] not in TERMINAL]
    if live:
        a.gather("progress", when="any_terminal", rejoin=[label for label in live if label in a.children.labels])
        return
    finish(a, parts, decisions)


if __name__ == "__main__":
    raise SystemExit(run.main())
