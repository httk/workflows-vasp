#!/usr/bin/env python3
"""The Birch–Murnaghan equation of state of one structure, by calling the leaf VASP workflows.

This runner runs no VASP itself. Through ``[workflow.calls]`` it composes child
jobs of the sibling packages, so every VASP run has its leaf's templates and
remedy ladder:

* ``start`` checks the parameters and, with ``relax`` (default true), calls
  ``vasp.relax`` on the input structure with this package's tight
  ``templates/INCAR.relax`` (child ``relax``); with ``relax = false`` the input
  structure is the reference as it is.
* ``deform`` copies the reference structure to ``relaxed/POSCAR``, computes one
  automatic k-point mesh on it, scales its lattice rows by ``1 + ε`` for each of
  ``number_of_frames`` linear strains ``ε`` spread over ``linear_strain``
  (volume ``V₀(1 + ε)³``), and calls ``vasp.relax-ions`` with ``ISIF = 4`` once
  per volume (children ``v00``, ``v01``, …): the cell shape and the atoms relax
  at fixed volume.
* ``statics`` (only with ``static = true``) calls ``vasp.static`` on each
  relaxed cell (children ``s00``, ``s01``, …).
* ``fit`` reads every child's INCAR, KPOINTS and remedy history, checks that
  they ran one numerical setting, takes each volume from its ``vasp.relax-ions``
  CONTCAR and each energy from the OUTCAR of the child that provides it (the
  static when ``static``), and fits httk-analyse's ``fit_birch_murnaghan``.

Every child gets the same ``kpoint_grid``, ``centering`` and ``excluded_remedies``
(no remedy may change the mesh, the smearing, ``LREAL`` or the cell), and the
shared parameters of this job. ``eos.json`` in the workdir records the plan,
every child's identity, settings and outcome, the volumes and energies and the
fitted values; ``collect.py`` reads only it and ``relaxed/POSCAR``. A failed
child fails this job with ``vasp.dependency_failed``; ``httk job continue`` reruns it.
"""

import dataclasses
import json
import math
import re
import shutil
from pathlib import Path

from httk.codes.vasp import automatic_kpoint_grid, job_remedy_history_path, read_incar
from httk.workflow import Attempt, Runner

TEMPLATES = Path(__file__).resolve().parent / "templates"
RESULT = "eos.json"
RELAXED = "relaxed/POSCAR"
# The job parameters every child gets when this job sets them.
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
# The ladder rungs that would change the numerical setting of one child against the others.
EXCLUDED_REMEDIES = "bump_kpoints equal_kpoints centering incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice"
# Tags a remedy may change in one child only: they change the electronic path, not the converged energy surface
# (EDIFF only tightens, by scale_ediff; NEDOS only samples the DOS, by the dentet rung).
ELECTRONIC_PATH_TAGS = {
    *("ALGO", "TIME", "AMIX", "BMIX", "AMIN", "ICHARG", "ISTART", "NELM", "NELMIN", "NBANDS", "NPAR", "NCORE"),
    *("KPAR", "ISYM", "SYMPREC", "EDIFF", "NEDOS"),
}

run = Runner("vasp.eos")


def flag(value: object, name: str) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line.

    :param value: The parameter value.
    :param name: The parameter name, for the error.
    :return: The boolean.
    :raises ValueError: If *value* is neither true nor false.
    """

    text = str(value).strip().lower()
    if text not in ("true", "1", "yes", "false", "0", "no"):
        raise ValueError(f"{name} must be true or false, not {value!r}")
    return text in ("true", "1", "yes")


def settings(a: Attempt) -> dict:
    """Return this job's own validated parameters.

    :param a: The running attempt.
    :return: The parameters by name.
    :raises ValueError: If a parameter is invalid.
    """

    strain = a.parameter("linear_strain", [-0.05, 0.05])
    words = strain.split() if isinstance(strain, str) else strain
    try:
        lo, hi = (float(word) for word in words)
    except (TypeError, ValueError):
        raise ValueError(f"linear_strain must be two numbers, not {strain!r}") from None
    if not (-0.3 < lo < 0 < hi < 0.3):
        raise ValueError(f"linear_strain must span zero within (-0.3, 0.3) for an interior minimum, not {strain!r}")
    frames = a.parameter("number_of_frames", 6)
    try:
        count = int(str(frames))
    except ValueError:
        raise ValueError(f"number_of_frames must be an integer, not {frames!r}") from None
    if count < 5:
        raise ValueError(f"number_of_frames must be at least 5 for a Birch–Murnaghan fit, not {count}")
    density = a.parameter("kpoint_density", 20.0)
    try:
        density = float(density)
    except (TypeError, ValueError):
        raise ValueError(f"kpoint_density must be a number, not {density!r}") from None
    if not (math.isfinite(density) and density > 0):
        raise ValueError(f"kpoint_density must be positive, not {density!r}")
    centering = a.parameter("centering", "Gamma")
    if centering not in ("Gamma", "Monkhorst-Pack"):
        raise ValueError(f"centering must be Gamma or Monkhorst-Pack, not {centering!r}")
    return {
        "relax": flag(a.parameter("relax", True), "relax"),
        "static": flag(a.parameter("static", False), "static"),
        "publish_data": flag(a.parameter("publish_data", False), "publish_data"),
        "strains": [lo + (hi - lo) * i / (count - 1) for i in range(count)],
        "kpoint_density": density,
        "centering": centering,
    }


def shared_parameters(a: Attempt) -> dict:
    parameters = dict(a.parameters)
    return {name: parameters[name] for name in SHARED if parameters.get(name) is not None}


def potcar_files(a: Attempt) -> dict[str, Path]:
    potcar = a.payload / str(a.parameter("potcar", None) or "files/POTCAR")
    return {"POTCAR": potcar} if potcar.is_file() else {}


def read_result(a: Attempt) -> dict:
    return json.loads((a.workdir / RESULT).read_text(encoding="utf-8"))


def write_result(a: Attempt, result: dict) -> None:
    path = a.workdir / RESULT
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(result, indent=1), encoding="utf-8")
    temporary.replace(path)


def volume(poscar: Path) -> float:
    """The cell volume of a POSCAR with a positive scale, in A^3: |det| times the scale cubed."""

    lines = poscar.read_text(encoding="utf-8").splitlines()
    m = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    determinant = (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )
    return abs(determinant) * float(lines[1].split()[0]) ** 3


def unsupported(lines: list[str]) -> str | None:
    """Why the lattice of these POSCAR lines cannot be edited as text, or ``None`` when it can."""

    if len(lines) < 8 or not lines[5].split() or any(re.fullmatch(r"\d+", word) for word in lines[5].split()):
        return "the reference POSCAR has no species line (a VASP-5 POSCAR is needed)"
    try:
        scale = float(lines[1].split()[0])
    except (IndexError, ValueError):
        return "the reference POSCAR has no scale line"
    if scale != 1.0:
        return f"the reference POSCAR has the scale {lines[1].strip()}, not 1"
    mode = lines[7].strip()[:1].lower()
    if mode == "s":
        return "the reference POSCAR has a Selective dynamics line"
    if mode != "d":
        return f"the reference POSCAR has {lines[7].strip()!r} coordinates, not Direct"
    return None


def child_record(child) -> dict:
    """What eos.json keeps of a finished child: identity, outcome, and the settings it ran."""

    record = {
        "label": child.label,
        "job_key": child.job_key,
        "job_id": child.job_id,
        "kind": child.kind,
        "failure": None if child.failure is None else {"code": child.failure.code, "message": child.failure.message},
    }
    if child.workdir is not None and (child.workdir / "INCAR").is_file():
        record["incar"] = read_incar(child.workdir / "INCAR")
        kpoints = child.workdir / "KPOINTS"
        record["kpoints"] = kpoints.read_text(encoding="utf-8") if kpoints.is_file() else None
    history = None if child.payload is None else job_remedy_history_path(child.payload)
    events = json.loads(history.read_text(encoding="utf-8")).get("events", []) if history and history.is_file() else []
    record["remedies"] = [change["operation"] for event in events for change in event.get("changes", [])]
    return record


def failed_children(a: Attempt, children, result: dict) -> bool:
    """Record *children* in *result*; fail the job and return true when one did not succeed or has no workdir."""

    for child in children:
        result["children"][child.label] = child_record(child)
    write_result(a, result)
    failed = [child for child in children if not child.succeeded]
    if failed:
        a.fail(
            "vasp.dependency_failed",
            f"{len(failed)} of {len(children)} children did not succeed: {', '.join(c.label for c in failed)}",
            details={"failed": [result["children"][c.label] for c in failed]},
        )
        return True
    missing = [child.label for child in children if child.workdir is None]
    if missing:
        a.fail("vasp.child_workdir_missing", f"no workdir found for {', '.join(missing)}", details={"labels": missing})
        return True
    return False


def inconsistency(records: list[dict]) -> dict | None:
    """The settings that differ among *records* (KPOINTS, or INCAR tags off the electronic path), or ``None``."""

    def normal(value: object) -> str | None:
        return None if value is None else str(value).strip().lower()

    differing: dict[str, dict[str, str | None]] = {}
    if len({normal(record.get("kpoints")) for record in records}) > 1:
        differing["KPOINTS"] = {record["label"]: record.get("kpoints") for record in records}
    tags = {tag for record in records for tag in record.get("incar", {})} - ELECTRONIC_PATH_TAGS
    for tag in sorted(tags):
        values = {record["label"]: record.get("incar", {}).get(tag) for record in records}
        if len({normal(value) for value in values.values()}) > 1:
            differing[tag] = values
    return differing or None


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def start(a: Attempt) -> None:
    """Check the parameters, and call the tight relaxation or take the input structure as the reference."""

    try:
        plan = settings(a)
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    poscar = a.payload / str(a.parameter("poscar", None) or "files/POSCAR")
    if not poscar.is_file():
        a.fail("vasp.input_missing", f"the starting structure {poscar.name} is not in this payload")
        return
    write_result(a, {"workflow": "vasp.eos", **plan, "children": {}})
    if not plan["relax"]:
        (a.workdir / RELAXED).parent.mkdir(exist_ok=True)
        shutil.copyfile(poscar, a.workdir / RELAXED)
        a.log.append("note", "relax = false: the input structure is the reference")
        a.advance("deform")
        return
    override = a.payload / "files" / "INCAR.relax"
    incar = override if override.is_file() else TEMPLATES / "INCAR.relax"
    files = {"POSCAR": poscar, "INCAR.relax": incar} | potcar_files(a)
    a.call("relax", label="relax", files=files, parameters=shared_parameters(a))
    a.log.append("note", f"called the tight relaxation with {incar.name} from {incar.parent.name}/")
    a.gather("deform", when="all_succeeded", on_impossible="relax_failed")


@run.step
def relax_failed(a: Attempt) -> None:
    """Fail by the relaxation that failed, with the child's failure."""

    child = a.children.get("relax")
    record = {} if child is None else child_record(child)
    a.fail("vasp.relax_failed", "the tight relaxation did not succeed", details=record)


@run.step
def deform(a: Attempt) -> None:
    """Scale the reference cell to every volume and call one fixed-volume relaxation per volume."""

    result = read_result(a)
    child = a.children.get("relax")
    if child is not None:
        result["children"]["relax"] = child_record(child)
        write_result(a, result)
        if child.workdir is None:
            a.fail("vasp.child_workdir_missing", "no workdir found for relax", details={"labels": ["relax"]})
            return
        contcar = child.workdir / "CONTCAR"
        if not contcar.is_file():
            a.fail("vasp.relax_failed", f"the relaxation {child.job_key} left no CONTCAR")
            return
        (a.workdir / RELAXED).parent.mkdir(exist_ok=True)
        shutil.copyfile(contcar, a.workdir / RELAXED)
    reference = a.workdir / RELAXED
    lines = reference.read_text(encoding="utf-8").splitlines(keepends=True)
    problem = unsupported(lines)
    if problem is not None:
        a.fail("vasp.input_unsupported", problem)
        return
    grid = " ".join(map(str, automatic_kpoint_grid(result["kpoint_density"], poscar=reference)))
    rows = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    parameters = shared_parameters(a) | {
        "isif": 4,
        "kpoint_grid": grid,
        "centering": result["centering"],
        "excluded_remedies": EXCLUDED_REMEDIES,
    }
    result.update(reference_volume=volume(reference), kpoint_grid=grid, excluded_remedies=EXCLUDED_REMEDIES)
    result["points"] = []
    for index, strain in enumerate(result["strains"]):
        label = f"v{index:02d}"
        path = a.workdir / "volumes" / label / "POSCAR"
        path.parent.mkdir(parents=True, exist_ok=True)
        edited = list(lines)
        for i, row in zip((2, 3, 4), rows, strict=True):
            edited[i] = "  " + "  ".join(f"{(1 + strain) * x:.16f}" for x in row) + "\n"
        path.write_text("".join(edited), encoding="utf-8")
        a.call("relax_ions", label=label, files={"POSCAR": path} | potcar_files(a), parameters=parameters)
        result["points"].append(
            {"label": label, "linear_strain": strain, "planned_volume": result["reference_volume"] * (1 + strain) ** 3}
        )
    write_result(a, result)
    a.log.append(
        "note",
        f"called {len(result['strains'])} fixed-volume relaxations (ISIF = 4) on the mesh {grid} "
        f"({result['centering']}) from a reference volume of {result['reference_volume']:.4f} A^3",
    )
    a.gather("statics" if result["static"] else "fit", when="all_terminal")


@run.step
def statics(a: Attempt) -> None:
    """Call one single point on each relaxed cell, and wait for them and the relaxations."""

    result = read_result(a)
    relaxations = [a.children[point["label"]] for point in result["points"]]
    if failed_children(a, relaxations, result):
        return
    parameters = shared_parameters(a) | {
        "kpoint_grid": result["kpoint_grid"],
        "centering": result["centering"],
        "excluded_remedies": EXCLUDED_REMEDIES,
    }
    for index, (point, child) in enumerate(zip(result["points"], relaxations, strict=True)):
        point["static_label"] = label = f"s{index:02d}"
        files = {"POSCAR": child.workdir / "CONTCAR"} | potcar_files(a)
        a.call("static", label=label, files=files, parameters=parameters)
    write_result(a, result)
    a.log.append("note", f"called {len(relaxations)} single points on the relaxed cells")
    a.gather("fit", when="all_terminal", rejoin=[point["label"] for point in result["points"]])


@run.step
def fit(a: Attempt) -> None:
    """Check the children ran one setting, read the volumes and energies, and fit the equation of state."""

    from httk.atomistic.integrations.vasp.io import OutcarFile

    result = read_result(a)
    if failed_children(a, list(a.children), result):
        return
    source = "static" if result["static"] else "relax_ions"
    providers = [point.get("static_label", point["label"]) for point in result["points"]]
    differing = inconsistency([result["children"][label] for label in providers])
    if differing is not None:
        a.fail(
            "vasp.inconsistent_children",
            f"the children did not run one numerical setting: {', '.join(differing)} differ",
            details=differing,
        )
        return
    for point, provider in zip(result["points"], providers, strict=True):
        point["volume"] = volume(a.children[point["label"]].workdir / "CONTCAR")
        outcar = a.children[provider].workdir / "OUTCAR"
        energy = OutcarFile(outcar).final_energies.energy_sigma0 if outcar.is_file() else None
        if energy is None:
            write_result(a, result)
            a.fail("vasp.energy_missing", f"the child {provider} left no energy(sigma->0) in its OUTCAR")
            return
        point["energy"] = float(energy)
    result["energy_source"] = source
    volumes = [point["volume"] for point in result["points"]]
    energies = [point["energy"] for point in result["points"]]
    try:
        from httk.analyse.matsci.eos import fit_birch_murnaghan
    except ImportError as exception:
        write_result(a, result)
        a.fail("vasp.analysis_unavailable", f"fitting the equation of state needs httk-analyse: {exception}")
        return
    try:
        fitted = fit_birch_murnaghan(volumes, energies)
    except ValueError as exception:
        write_result(a, result)
        a.fail("vasp.eos_fit_failed", f"the Birch–Murnaghan fit failed: {exception}")
        return
    # Every field, so collect.py rebuilds the fit without refitting.
    result["fit"] = dataclasses.asdict(fitted)
    write_result(a, result)
    if result["publish_data"]:
        prefix = a.parameter("data_prefix", "") or ""
        for name in (RESULT, RELAXED):
            a.put(a.workdir / name, f"{prefix}/{name}" if prefix else name)
    headline = ("equilibrium_volume", "equilibrium_energy", "bulk_modulus", "bulk_modulus_derivative", "rmse")
    a.state.merge({name: result["fit"][name] for name in headline} | {"energy_source": source})
    a.log.append(
        "note",
        f"Birch–Murnaghan fit of {len(volumes)} {source} energies: V0 = {fitted.equilibrium_volume:.4f} A^3, "
        f"E0 = {fitted.equilibrium_energy:.6f} eV, B0 = {fitted.bulk_modulus:.2f} GPa, "
        f"B0' = {fitted.bulk_modulus_derivative:.3f}, RMSE = {fitted.rmse:.2e} eV",
    )
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
