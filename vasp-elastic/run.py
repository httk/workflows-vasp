#!/usr/bin/env python3
"""The elastic stiffness tensor of a structure, by calling ``vasp.relax`` and ``vasp.relax-ions``.

This runner runs no VASP itself. It composes child jobs through httk's workflow
calls (``[workflow.calls]``, ``Attempt.call`` and ``Attempt.gather``):

* ``start`` validates the parameters and, with ``relax`` (default true), calls
  ``vasp.relax`` (child ``relax``) with this package's tight
  ``templates/INCAR.relax`` for its relax rounds; with ``relax = false`` the
  input structure is taken as relaxed.
* ``deform`` copies the relaxed structure to ``relaxed/POSCAR``, computes one
  automatic k-point mesh on it, and writes one strained cell per strain state
  and magnitude, ``deformations/<label>/POSCAR``: the Voigt strain ``e = m s``
  (``xx yy zz yz xz xy``, engineering shear strains), the symmetric small strain
  ``eps`` and the lattice rows ``a_i' = (I + eps) a_i``. It calls
  ``vasp.relax-ions`` (``ISIF = 2``) once per cell, all on that mesh and with the
  numerics-changing remedies excluded, and waits for all of them.
* ``fit`` records every child's INCAR, KPOINTS, remedies, strain and stress
  (OUTCAR, tensile-positive GPa) in ``elastic.json``, refuses children that ran
  with different settings, and fits the stiffness and a stress offset by least
  squares with httk-analyse's ``fit_stress_strain``.

``elastic.json`` holds every fit input and the fitted values, and ``collect.py``
reads only it and ``relaxed/POSCAR``. The job parameters are documented in this
repository's README.
"""

import json
import math
import shutil
from pathlib import Path

import numpy
from httk.codes.vasp import KPOINT_CENTERINGS, automatic_kpoint_grid, job_remedy_history_path, read_incar
from httk.workflow import Attempt, Runner

TEMPLATES = Path(__file__).resolve().parent / "templates"
RESULT = "elastic.json"
RELAXED = Path("relaxed") / "POSCAR"
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
# Remedies may change these between children: they do not change the converged energy surface
# (scale_ediff tightens EDIFF, the dentet rung sets NEDOS).
ELECTRONIC_TAGS = frozenset(
    {
        "EDIFF",
        "NEDOS",
        "ALGO",
        "TIME",
        "AMIX",
        "BMIX",
        "AMIN",
        "ICHARG",
        "ISTART",
        "NELM",
        "NELMIN",
        "NBANDS",
        "NPAR",
        "NCORE",
        "KPAR",
        "ISYM",
        "SYMPREC",
    }
)
# atomate2 generate_elastic_deformations: linspace(-0.01, 0.01, 5) without zero.
DEFAULT_MAGNITUDES = (-0.01, -0.005, 0.005, 0.01)
DEFAULT_STATES = tuple(tuple(1.0 if i == j else 0.0 for j in range(6)) for i in range(6))

run = Runner("vasp.elastic")


# -- parameters -------------------------------------------------------------------------------------------------------


def flag(value: object) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line."""

    text = str(value).strip().lower()
    if text not in ("true", "1", "yes", "false", "0", "no"):
        raise ValueError(f"relax must be true or false, not {value!r}")
    return text in ("true", "1", "yes")


def real(value: object, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{what} must be a number, not {value!r}")  # noqa: TRY004 (a bad parameter fails vasp.input_invalid)
    try:
        number = float(value)
    except ValueError:
        raise ValueError(f"{what} must be a number, not {value!r}") from None
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite, not {value!r}")
    return number


def strain_magnitudes(value: object) -> list[float]:
    """Return the strain magnitudes *value* lists (space-separated, or a list).

    :raises ValueError: Unless they are at least two distinct nonzero numbers.
    """

    words = value.split() if isinstance(value, str) else value
    if not isinstance(words, (list, tuple)):
        raise ValueError(f"strain_magnitudes must be a list of numbers, not {value!r}")  # noqa: TRY004 (as above)
    magnitudes = [real(word, "a strain magnitude") for word in words]
    if 0.0 in magnitudes:
        raise ValueError("strain_magnitudes may not contain zero: an unstrained cell adds nothing to the fit")
    if len(magnitudes) < 2 or len(set(magnitudes)) != len(magnitudes):
        raise ValueError(
            f"strain_magnitudes must be at least two distinct numbers, so the fit can separate the stress offset, "
            f"not {value!r}"
        )
    return magnitudes


def strain_states(value: object) -> list[list[float]]:
    """Return the strain states *value* lists (a list of six-number lists, or a JSON string of one).

    :raises ValueError: Unless every state has six numbers, none is zero, and together they span all six.
    """

    if isinstance(value, str):
        try:
            value = json.loads(value)
        except ValueError:
            raise ValueError(f"strain_states must be a list of six-number lists, not {value!r}") from None
    if not isinstance(value, (list, tuple)) or not all(
        isinstance(state, (list, tuple)) and len(state) == 6 for state in value
    ):
        raise ValueError(f"strain_states must be a list of six-number lists, not {value!r}")
    states = [[real(x, "a strain state component") for x in state] for state in value]
    if any(not any(state) for state in states):
        raise ValueError("strain_states may not contain an all-zero state")
    if numpy.linalg.matrix_rank(numpy.array(states)) < 6:
        raise ValueError(
            "strain_states must span all six Voigt components (rank 6), or the 21 stiffness components cannot be fitted"
        )
    return states


def shared_parameters(a: Attempt) -> dict:
    parameters = dict(a.parameters)
    return {name: parameters[name] for name in SHARED if parameters.get(name) is not None}


def write_result(a: Attempt, result: dict) -> None:
    path = a.workdir / RESULT
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(result, indent=1), encoding="utf-8")
    temporary.replace(path)


def read_result(a: Attempt) -> dict:
    return json.loads((a.workdir / RESULT).read_text(encoding="utf-8"))


# -- deformations -----------------------------------------------------------------------------------------------------


def deformation_gradient(strain: list[float]) -> list[list[float]]:
    """Return ``F = I + eps`` for a Voigt strain ``xx yy zz yz xz xy`` with engineering shear components."""

    xx, yy, zz, yz, xz, xy = strain
    return [[1 + xx, xy / 2, xz / 2], [xy / 2, 1 + yy, yz / 2], [xz / 2, yz / 2, 1 + zz]]


def poscar_lattice(text: str) -> list[list[float]]:
    """Return the lattice rows of a POSCAR this runner can strain.

    :raises ValueError: Unless it is a VASP-5 POSCAR with scale 1, direct coordinates and no selective dynamics.
    """

    lines = text.splitlines()
    if len(lines) < 8:
        raise ValueError("the POSCAR is too short")
    scale = lines[1].split()
    if len(scale) != 1 or float(scale[0]) != 1.0:
        raise ValueError(f"the POSCAR scale line must be 1, not {lines[1].strip()!r}")
    if all(word.isdigit() for word in lines[5].split()):
        raise ValueError("the POSCAR has no species line (a VASP-4 POSCAR)")
    mode = lines[7].strip()[:1].lower()
    if mode != "d":
        raise ValueError(f"the POSCAR coordinates must be Direct, without Selective dynamics, not {lines[7].strip()!r}")
    return [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]


def strained_poscar(text: str, gradient: list[list[float]]) -> str:
    """Return *text* with its lattice rows ``a_i`` replaced by ``F a_i`` and every other line kept byte-identical."""

    lattice = poscar_lattice(text)
    lines = text.splitlines(keepends=True)
    for i, row in enumerate(lattice):
        new = [sum(gradient[k][j] * row[j] for j in range(3)) for k in range(3)]
        ending = lines[i + 2][len(lines[i + 2].rstrip("\r\n")) :]
        lines[i + 2] = "  " + "  ".join(f"{x:.16f}" for x in new) + ending
    return "".join(lines)


# -- children ---------------------------------------------------------------------------------------------------------


def child_record(child, strain: list[float] | None) -> dict:
    """What elastic.json keeps of one deformation child: identity, outcome, settings, strain and stress."""

    record: dict = {
        "label": child.label,
        "job_key": child.job_key,
        "job_id": child.job_id,
        "kind": child.kind,
        "failure": None if child.failure is None else {"code": child.failure.code, "message": child.failure.message},
        "strain": strain,
        "incar": None,
        "kpoints": None,
        "remedies": [],
        "stress": None,
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
        stress = None if frame is None else frame.stress_gpa_voigt()
        record["stress"] = None if stress is None else list(stress)
    return record


def inconsistencies(records: list[dict]) -> dict:
    """Return how the children's settings differ: KPOINTS texts, or INCAR tags other than the electronic-path ones.

    :param records: The children's records, with their ``incar`` tags and ``kpoints`` text.
    :return: ``{"kpoints": {label: text}, "incar_tags": {tag: {label: value}}}`` members for what differs; empty when
        the children agree.
    """

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
    """Validate the parameters, then call the tight relaxation, or with ``relax = false`` go on to the deformations."""

    try:
        relax = flag(a.parameter("relax", True))
        strain_states(a.parameter("strain_states", None) or [list(state) for state in DEFAULT_STATES])
        strain_magnitudes(a.parameter("strain_magnitudes", None) or list(DEFAULT_MAGNITUDES))
        density = real(a.parameter("deformation_kpoint_density", 40.0), "deformation_kpoint_density")
        if density <= 0:
            raise ValueError(f"deformation_kpoint_density must be positive, not {density}")
        centering = a.parameter("centering", "Gamma")
        if centering not in KPOINT_CENTERINGS:
            raise ValueError(f"centering must be one of {', '.join(KPOINT_CENTERINGS)}, not {centering!r}")
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail("vasp.input_missing", f"the starting structure {poscar.name} is not in this payload")
        return
    if not relax:
        (a.workdir / RELAXED).parent.mkdir(exist_ok=True)
        shutil.copyfile(poscar, a.workdir / RELAXED)
        a.log.append("note", "relax = false: the input structure is strained as given")
        a.advance("deform")
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
    a.gather("deform", when="all_succeeded", on_impossible="relax_failed")


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
def deform(a: Attempt) -> None:
    """Write the strained cells of the relaxed structure and call one fixed-cell ionic relaxation per cell."""

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
    text = relaxed.read_text(encoding="utf-8")
    try:
        reference = poscar_lattice(text)
    except ValueError as exception:
        a.fail("vasp.input_unsupported", f"cannot strain the relaxed structure: {exception}")
        return
    states = strain_states(a.parameter("strain_states", None) or [list(state) for state in DEFAULT_STATES])
    magnitudes = strain_magnitudes(a.parameter("strain_magnitudes", None) or list(DEFAULT_MAGNITUDES))
    centering = a.parameter("centering", "Gamma")
    grid = automatic_kpoint_grid(float(a.parameter("deformation_kpoint_density", 40.0)), poscar=relaxed)
    mesh = " ".join(map(str, grid))
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    parameters = shared_parameters(a) | {
        "kpoint_grid": mesh,
        "centering": centering,
        "excluded_remedies": EXCLUDED_REMEDIES,
        "isif": 2,
    }
    deformations = {}
    for index, strain in enumerate([m * x for x in state] for state in states for m in magnitudes):
        label = f"d{index:02d}"
        cell = a.workdir / "deformations" / label / "POSCAR"
        cell.parent.mkdir(parents=True, exist_ok=True)
        cell.write_text(strained_poscar(text, deformation_gradient(strain)), encoding="utf-8")
        files = {"POSCAR": cell} | ({"POTCAR": potcar} if potcar.is_file() else {})
        a.call("relax_ions", label=label, files=files, parameters=parameters)
        deformations[label] = strain
    write_result(
        a,
        {
            "workflow": "vasp.elastic",
            "relax": None if child is None else {"job_key": child.job_key, "job_id": child.job_id},
            "reference_lattice": reference,
            "kpoint_grid": mesh,
            "centering": centering,
            "excluded_remedies": EXCLUDED_REMEDIES,
            "strain_convention": "Voigt xx yy zz yz xz xy, engineering shear strains; F = I + eps, a_i' = F a_i",
            "deformations": deformations,
        },
    )
    a.log.append("note", f"called {len(deformations)} vasp.relax-ions deformations on the {mesh} {centering} mesh")
    a.gather("fit", when="all_terminal")


@run.step
def fit(a: Attempt) -> None:
    """Record every deformation child, check that they agree, and fit the elastic tensor to their stresses."""

    result = read_result(a)
    strains = result["deformations"]
    records = [child_record(child, strains.get(child.label)) for child in a.children]
    result["children"] = records
    write_result(a, result)
    failed = [record for record in records if record["kind"] != "succeeded"]
    if failed:
        a.fail(
            "vasp.dependency_failed",
            f"deformations {', '.join(record['label'] for record in failed)} did not succeed",
            details={"failed": {record["label"]: record["failure"] for record in failed}},
        )
        return
    missing = [child.label for child in a.children if child.workdir is None]
    if missing:
        a.fail("vasp.child_workdir_missing", f"deformations {', '.join(missing)} have no workdir")
        return
    differences = inconsistencies(records)
    if differences:
        a.fail(
            "vasp.inconsistent_children",
            "the deformations ran with different settings: "
            + ", ".join((["KPOINTS"] if "kpoints" in differences else []) + list(differences.get("incar_tags", {}))),
            details=differences,
        )
        return
    no_stress = [record["label"] for record in records if record["stress"] is None]
    if no_stress:
        a.fail("vasp.stress_missing", f"deformations {', '.join(no_stress)} have no stress in their OUTCAR")
        return
    try:
        from httk.analyse.matsci.elasticity import fit_stress_strain
    except ImportError:
        a.fail("vasp.analysis_unavailable", "fitting the elastic tensor needs httk-analyse: pip install httk-analyse")
        return
    try:
        fitted = fit_stress_strain(
            [record["strain"] for record in records], [record["stress"] for record in records], fit_offset=True
        )
    except ValueError as exception:
        a.fail("vasp.elastic_fit_failed", f"the stress-strain fit failed: {exception}")
        return
    tensor = fitted.tensor
    stable = tensor.is_stable()
    try:
        compliance = [list(row) for row in tensor.compliance]
    except ValueError:
        compliance = None
    result["fit"] = {
        "stiffness": [list(row) for row in tensor.stiffness],
        "stress_offset": list(fitted.stress_offset),
        "residuals": [list(row) for row in fitted.residuals],
        "rmse": fitted.rmse,
        "condition_number": fitted.condition_number,
        "strain_range": [list(pair) for pair in fitted.strain_range],
        "stable": stable,
        "stability_eigenvalues": list(tensor.stability_eigenvalues),
        "compliance": compliance,
    }
    headline: dict[str, object] = {"stable": stable, "rmse": fitted.rmse}
    if stable:
        moduli = {
            name: float(getattr(tensor, name))
            for name in (
                "bulk_modulus_voigt",
                "bulk_modulus_reuss",
                "bulk_modulus_hill",
                "shear_modulus_voigt",
                "shear_modulus_reuss",
                "shear_modulus_hill",
                "universal_anisotropy",
            )
        }
        result["fit"] |= moduli
        headline |= {
            "bulk_modulus_hill": moduli["bulk_modulus_hill"],
            "shear_modulus_hill": moduli["shear_modulus_hill"],
        }
        a.log.append(
            "note",
            f"fitted a stable elastic tensor: K_H = {moduli['bulk_modulus_hill']:.2f} GPa, "
            f"G_H = {moduli['shear_modulus_hill']:.2f} GPa, rmse {fitted.rmse:.3g} GPa",
        )
    else:
        a.log.append(
            "note",
            f"the fitted elastic tensor is NOT stable (smallest eigenvalue {min(tensor.stability_eigenvalues):.3g} GPa): "
            "no moduli",
        )
    write_result(a, result)
    a.state.merge(headline)
    if a.parameter("publish_data", False) is True:
        base = a.parameter("data_prefix", "") or ""
        for name in (RESULT, RELAXED.as_posix()):
            a.put(a.workdir / name, f"{base}/{name}" if base else name)
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
