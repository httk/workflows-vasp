#!/usr/bin/env python3
"""One density-functional perturbation theory run (``IBRION = 8``, ``LEPSILON``) of a fixed structure.

``prepare`` stages the structure and the INCAR template ``templates/INCAR.dielectric``
(the electronic settings of ``vasp-static`` with the DFPT switches of atomate2's
``StaticSetGenerator(lepsilon=True)``), ``run`` executes VASP with the reviewed remedy
ladder, and ``publish`` publishes the result. Between them ``run`` reads the response from the
OUTCAR and writes ``dielectric.json`` (see this repository's README): a run that printed no
dielectric tensor fails ``vasp.dielectric_missing``. Insulators only.

``IBRION``, ``LEPSILON``, ``NCORE`` and ``NPAR`` are what this workflow is (VASP forbids the
last two with ``IBRION = 8``), so ``incar_tags`` may not set them
(``vasp.unsupported_incar_tag``) nor ``parallel_tag`` name the last two
(``vasp.unsupported_parallel_tag``); ``KPAR`` is fine. The structure never changes by a remedy:
``scale_lattice``, ``contcar_to_poscar`` and ``incar.NPAR`` are never applied, nor the operations
``excluded_remedies`` names. ``kpoint_grid = "n1 n2 n3"`` fixes the automatic mesh; the default
centering is Gamma, which the tetrahedron method (``ISMEAR = -5``) wants.

Nothing here imports anything but an installed *httk-workflow*, *httk-atomistic* and
*httk-workflow-vasp*; reference this package directory by git URI (e.g.
``git+https://github.com/httk/workflows-vasp#vasp-dielectric``) or with ``--workflow-dir``. The
runner file alone, without its ``templates/``, fails with ``vasp.template_missing``.
"""

import hashlib
import json
import os
import re
import shlex
import shutil
from pathlib import Path

from httk.atomistic.integrations.vasp.io import OutcarFile
from httk.codes.vasp import (
    REMEDY_OPERATIONS,
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    read_incar,
    read_poscar_header,
    run_vasp,
    validate_vasp_workdir,
    write_automatic_kpoints,
)
from httk.workflow import Attempt, Runner

COLLECT = "dielectric.json POSCAR INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# Kept across a remedied rerun: they make the rerun cheaper, and VASP overwrites
# them itself when it reuses them.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")
TEMPLATES = Path(__file__).resolve().parent / "templates"

DEFINING_TAGS = ("IBRION", "LEPSILON", "NCORE", "NPAR")
VOIGT_ORDER = ["XX", "YY", "ZZ", "XY", "YZ", "ZX"]

run = Runner("vasp.dielectric")


def kpoint_grid(value: object) -> tuple[int, ...] | None:
    """Return the mesh *value* names (``"n1 n2 n3"`` or three integers), or ``None`` when it names none."""

    words = value.split() if isinstance(value, str) else value
    if not isinstance(words, list) or len(words) != 3 or not all(type(word) in (int, str) for word in words):
        return None
    try:
        grid = tuple(int(word) for word in words)
    except ValueError:
        return None
    return grid if min(grid) > 0 else None


def excluded_remedies(value: object) -> set[str]:
    """Return the ladder operation names *value* lists (space-separated, or a list).

    :param value: The ``excluded_remedies`` parameter.
    :return: The operation names.
    :raises ValueError: If *value* lists anything but remedy operation names.
    """

    names = value.split() if isinstance(value, str) else value or []
    if not isinstance(names, list) or not all(
        isinstance(name, str) and (name in REMEDY_OPERATIONS or re.fullmatch(r"incar\.[A-Z][A-Z0-9_]*", name))
        for name in names
    ):
        raise ValueError(
            f"excluded_remedies must list remedy operations ({', '.join(REMEDY_OPERATIONS)} or incar.<TAG>), "
            f"not {value!r}"
        )
    return set(names)


def poscar_cell(path: Path) -> tuple[list[list[float]], list[list[float]]]:
    """Return the scaled lattice (Å) and the direct coordinates of the POSCAR at *path*.

    :raises ValueError: If the POSCAR has a non-positive scale or Cartesian coordinates.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    header = read_poscar_header(path)
    if header.scale <= 0:
        raise ValueError("a POSCAR with a negative (volume) scale is not supported")
    mode = lines[7].strip().lower()
    if mode.startswith("s"):
        mode = lines[8].strip().lower()
        first = 9
    else:
        first = 8
    if not mode.startswith("d"):
        raise ValueError("a Cartesian POSCAR is not supported: write direct coordinates")
    positions = [[float(x) for x in line.split()[:3]] for line in lines[first : first + sum(header.counts)]]
    return [[header.scale * x for x in row] for row in header.lattice], positions


def write_response(workdir: Path) -> str | None:
    """Write ``dielectric.json`` from the OUTCAR of *workdir*; return what is wrong, or ``None``."""

    outcar = OutcarFile(workdir / "OUTCAR")
    response = outcar.dielectric_response
    electronic = None if response is None else response.electronic_floats()
    if response is None or electronic is None:
        return (
            "the OUTCAR has no macroscopic static dielectric tensor: is the structure an insulator that VASP finished?"
        )
    ionic, born = response.ionic_floats(), response.born_charges_floats()
    lexeme = getattr(getattr(outcar, "final_energies", None), "energy_sigma0", None)
    header = read_poscar_header(workdir / "POSCAR")
    lattice, positions = poscar_cell(workdir / "POSCAR")
    if born is not None and len(born) != sum(header.counts):
        return f"the OUTCAR has {len(born)} Born charge tensors for the POSCAR's {sum(header.counts)} atoms"
    document = {
        "electronic": electronic,
        "ionic": ionic,
        "static": None if ionic is None else [[x + y for x, y in zip(r, q)] for r, q in zip(electronic, ionic)],
        "born_charges": born,
        "piezoelectric_electronic": response.piezoelectric_electronic_floats(),
        "piezoelectric_ionic": response.piezoelectric_ionic_floats(),
        "piezoelectric_column_order": VOIGT_ORDER,
        "lattice": lattice,
        "positions": positions,
        "species": list(header.species),
        "counts": list(header.counts),
        "poscar_sha256": hashlib.sha256((workdir / "POSCAR").read_bytes()).hexdigest(),
        "energy": None if lexeme is None else float(lexeme),
    }
    temporary = workdir / "dielectric.json.tmp"
    temporary.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, workdir / "dielectric.json")
    return None


@run.step
def prepare(a: Attempt) -> None:
    """Stage the structure and the INCAR template, derive the rest, and go on to run VASP."""

    validate_vasp_workdir(a.workdir)
    incar = a.workdir / "INCAR"
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail(
            "vasp.input_missing",
            f"the starting structure {poscar.name} is not in this payload",
            details={"expected": str(poscar)},
        )
        return
    for tag in dict(a.parameter("incar_tags", {}) or {}):
        if str(tag).strip().upper() in DEFINING_TAGS:
            a.fail(
                "vasp.unsupported_incar_tag",
                f"{tag} from incar_tags is not supported: this workflow is a DFPT run (IBRION = 8, LEPSILON), "
                "which VASP does not combine with NPAR or NCORE (use KPAR to parallelize)",
                details={"tag": str(tag).strip().upper(), "source": "incar_tags"},
            )
            return
    parallel_tag = str(a.parameter("parallel_tag", None) or "")
    if parallel_tag.strip().upper() in ("NPAR", "NCORE"):
        a.fail(
            "vasp.unsupported_parallel_tag",
            f"parallel_tag {parallel_tag} is not supported: VASP forbids NPAR and NCORE with IBRION = 8 (use KPAR)",
            details={"tag": parallel_tag.strip().upper()},
        )
        return
    grid_parameter = a.parameter("kpoint_grid", "") or ""
    grid = kpoint_grid(grid_parameter) if grid_parameter else ()
    if grid is None:
        a.fail(
            "vasp.input_invalid",
            f"kpoint_grid must be three positive integers, not {grid_parameter!r}",
            details={"kpoint_grid": grid_parameter},
        )
        return
    try:
        excluded_remedies(a.parameter("excluded_remedies", ""))
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    shutil.copyfile(poscar, a.workdir / "POSCAR")
    try:
        poscar_cell(a.workdir / "POSCAR")
    except (ValueError, IndexError) as exception:
        a.fail("vasp.input_unsupported", str(exception))
        return
    override = a.payload / "files" / "INCAR.dielectric"
    source = override if override.is_file() else TEMPLATES / "INCAR.dielectric"
    if not source.is_file():
        a.fail(
            "vasp.template_missing",
            f"the INCAR template {source} is missing: reference this workflow as its package directory "
            "(a git URI ending in #vasp-dielectric, or --workflow-dir), not as its runner file alone",
            details={"expected": str(source)},
        )
        return
    shutil.copyfile(source, incar)
    template = read_incar(incar)
    library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
        library = None
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(a.parameter("kpoint_density", 20.0) or 20.0),
        centering=a.parameter("centering", "Gamma"),
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        pseudopotential_library=library,
        parallel_tag=a.parameter("parallel_tag", None) or None,
        parallel_value=None if parallel_value is None else int(parallel_value),
        # The template's own EDIFF/EDIFFG are kept rather than derived; the job's tags win over both.
        incar_tags={
            **{tag: template[tag] for tag in ("EDIFF", "EDIFFG") if tag in template},
            **dict(a.parameter("incar_tags", {}) or {}),
        },
    )
    prepare_vasp_inputs(options, directory=a.workdir)
    if grid:
        write_automatic_kpoints(grid, a.workdir / "KPOINTS", centering=options.centering)
        a.log.append("note", f"fixed the k-point mesh at {' '.join(map(str, grid))} ({options.centering})")
    a.log.append("note", "prepared a vasp.dielectric calculation")
    a.advance("run")


@run.step(name="run")
def run_step(a: Attempt) -> None:
    """Run VASP, remedy a recognized failure, or fail with what was diagnosed."""

    # vasp.command resolves the job parameter, HTTK_VASP_COMMAND, then the workspace setting.
    argv = shlex.split(a.setting("vasp.command", a.parameter("vasp_command", None)) or "")
    if not argv:
        a.fail(
            "vasp.command_missing",
            "no VASP command is configured: set it with `httk workspace settings set --key vasp.command --value '...' WORKSPACE`, "
            "or set HTTK_VASP_COMMAND on the machine that runs this job, or give the job a vasp_command parameter",
        )
        return
    # A rerun must not read the previous run's outputs; CONTCAR and the run report
    # survive on purpose, since a remedy and a restart are derived from them.
    clean_vasp_outputs(a.workdir, keep=KEEP_BETWEEN_RUNS)
    try:
        report = run_vasp(argv, directory=a.workdir, timeout=a.parameter("timeout", 86400.0))
    except OSError as exception:
        a.fail("vasp.failed", f"could not run VASP at all: {exception}")
        return
    a.log.append("note", f"VASP {report.classification}")
    oszicar = a.workdir / "OSZICAR"
    energy = last_oszicar_energy(oszicar) if oszicar.is_file() else None
    state: dict[str, object] = {"classification": report.classification}
    if energy is not None:
        state["energy"] = energy
    if report.classification == "completed":
        message = write_response(a.workdir)
        if message:
            a.state.merge(state)
            a.fail("vasp.dielectric_missing", message)
            return
        a.advance("publish", state=state)
        return

    applied = int(a.state.get("remedies", 0))
    maximum = int(a.parameter("maximum_remedies", 8))
    history = job_remedy_history_path(a.payload)
    # Planned before the budget is consulted, so a job that stops here still says
    # which remedy it would have applied.
    try:
        decision = plan_vasp_remedy(
            report.diagnostics,
            directory=a.workdir,
            history_path=history,
            policy=a.parameter("remedy_policy", "reviewed-v1"),
            excluded_operations={"incar.NPAR", "scale_lattice", "contcar_to_poscar"}
            | excluded_remedies(a.parameter("excluded_remedies", "")),
        )
    except ValueError as exception:
        a.state.merge(state)
        a.fail("vasp.failed", f"planning a VASP remedy failed: {exception}")
        return
    if decision.give_up or applied >= maximum:
        a.state.merge(state)
        if applied >= maximum:
            message = f"VASP {report.classification} after {applied} remedies"
        else:
            message = f"VASP {report.classification} with no remaining remedy"
        a.fail("vasp.failed", message, details=decision.as_mapping())
        return
    apply_vasp_remedy(decision, directory=a.workdir, history_path=history)
    state["remedies"] = applied + 1
    a.state.merge(state)
    a.log.append("note", f"applied a remedy for {decision.problem}")
    a.retry(f"applied the {decision.policy} remedy for {decision.problem}")


@run.step
def publish(a: Attempt) -> None:
    """Publish the finished calculation and complete the job."""

    prefix = a.parameter("data_prefix", "vasp") or ""
    # Results stay in the persistent workdir; publish_data opts a job into a curated data/ copy too.
    to_data = a.parameter("publish_data", False) is True
    published = []
    for name in (a.parameter("collect", COLLECT) or "").split():
        if (a.workdir / name).is_file():
            if to_data:
                a.put(a.workdir / name, f"{prefix}/{name}")
            published.append(name)
    if to_data:
        a.log.append("note", f"published to data/{prefix}: {', '.join(published) or 'nothing'}")
    else:
        # Without publish_data the persistent workdir is the result.
        a.log.append("note", f"kept in the workdir: {', '.join(published) or 'nothing'}")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
