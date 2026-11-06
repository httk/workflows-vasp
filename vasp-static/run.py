#!/usr/bin/env python3
"""One single-point VASP calculation of a structure that is already chosen.

The workflow is spelled out on the ``httk.codes.vasp`` primitives so that this
package directory is a starting point to copy and edit: ``prepare`` stages the structure and
the INCAR template ``templates/INCAR.static`` — the electronic settings of
``vasp-relax``'s ``INCAR.relax`` with the ionic loop switched off
(``NSW = 0``, ``IBRION = -1``) — and derives everything else, ``run`` executes
VASP with the reviewed remedy ladder, and ``publish`` publishes the result. The
structure may be a POSCAR or the CONTCAR of an earlier relaxation, staged as an
input file of this job and named by the ``poscar`` input.

Edit the template of a copied package to change the settings, or replace it in
a single job with a payload file ``files/INCAR.static``. The template's own
``EDIFF`` and ``EDIFFG`` are kept; when absent they are derived from
``accuracy_per_atom``. The ``incar_tags`` parameter wins over the template and
every derived value. ``kpoint_grid = "n1 n2 n3"`` fixes the automatic mesh instead
of deriving it from ``kpoint_density``, and ``excluded_remedies`` names ladder
operations never to apply, so every child of a fan-out runs one numerical setting.

The job inputs and parameters are documented in this repository's README. Nothing
here imports anything but an installed *httk-workflow*, so this file and its
``templates/`` are the whole runner: it is the entry of this package directory,
reference it by git URI (e.g. ``git+https://github.com/httk/workflows-vasp#vasp-static``)
or with ``--workflow-dir`` (installed into the workspace first), or copy the package
directory and edit it. The runner file alone, without its ``templates/``, fails
with ``vasp.template_missing``.

A payload ``files/CHGCAR`` (and ``files/WAVECAR``; the ``chgcar`` and ``wavecar``
parameters) is copied into the workdir, and the template's ``ICHARG = 1`` then
reads it; without a non-empty CHGCAR it becomes ``ICHARG = 2``. An ``ICHARG``
from ``incar_tags`` is never changed, and one of 1 or 11 without a CHGCAR fails
``vasp.input_missing``. A payload ``files/scf-INCAR`` (``scf_incar``), the INCAR
the density was computed with, must agree on ``PREC``, ``ENCUT`` and ``ENAUG``
or the job fails ``vasp.charge_density_incompatible``.
"""

import re
import shlex
import shutil
from pathlib import Path

from httk.codes.vasp import (
    REMEDY_OPERATIONS,
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    parse_diagnostic_policy,
    parse_potcar_variants,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    rattle_poscar,
    read_checker_specs,
    read_incar,
    run_vasp,
    update_incar,
    validate_vasp_workdir,
    write_automatic_kpoints,
)
from httk.workflow import Attempt, Runner

COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# Kept across a remedied rerun: they make the rerun cheaper, and VASP overwrites
# them itself when it reuses them.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")
TEMPLATES = Path(__file__).resolve().parent / "templates"
# The INCAR tags that fix the FFT grid a CHGCAR is written on.
GRID_TAGS = ("PREC", "ENCUT", "ENAUG")

run = Runner("vasp.static")


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


def pinned_icharg(a: Attempt) -> set[str]:
    """Return ``{"incar.ICHARG"}`` when the job's ``incar_tags`` set ``ICHARG``, which no remedy may change."""

    tags = dict(a.parameter("incar_tags", {}) or {})
    return {"incar.ICHARG"} if any(str(tag).strip().upper() == "ICHARG" for tag in tags) else set()


def stage_charge_density(a: Attempt, explicit: dict[str, object]) -> bool:
    """Settle ``ICHARG`` against the staged ``CHGCAR`` and check the density, after input preparation.

    An ``ICHARG`` from *explicit* (the job's ``incar_tags``) is never changed; the
    template's ``ICHARG = 1`` becomes 2 without a staged ``CHGCAR`` to read.

    :param a: The attempt whose workdir holds the prepared inputs.
    :param explicit: The job's ``incar_tags``.
    :return: Whether the job may go on; when not, it has failed.
    """

    incar = a.workdir / "INCAR"
    chgcar = a.workdir / "CHGCAR"
    has_density = chgcar.is_file() and chgcar.stat().st_size > 0
    named = a.parameter("chgcar", "files/CHGCAR")
    expected = {"expected": str(a.payload / named)}
    requested = {str(tag).strip().upper(): value for tag, value in explicit.items()}
    if "ICHARG" in requested:
        value = requested["ICHARG"]
        if str(value).strip() in ("1", "11") and not has_density:
            a.fail(
                "vasp.input_missing",
                f"ICHARG = {str(value).strip()} from incar_tags reads a charge density, but this payload has no non-empty {named}",
                details={"ICHARG": value, **expected},
            )
            return False
    elif read_incar(incar).get("ICHARG") == "1" and not has_density:
        update_incar({"ICHARG": 2}, incar)
    scf_incar = a.payload / a.parameter("scf_incar", "files/scf-INCAR")
    if not scf_incar.is_file():
        return True
    if not has_density:
        a.fail(
            "vasp.input_missing",
            f"{scf_incar.name} describes a charge density, but this payload has no non-empty {named}",
            details=expected,
        )
        return False
    scf, staged = read_incar(scf_incar), read_incar(incar)

    def normalized(value: str | None) -> str | None:
        return None if value is None else value.strip().lower()

    differences = {
        tag: {"scf": scf.get(tag), "staged": staged.get(tag)}
        for tag in GRID_TAGS
        if normalized(scf.get(tag)) != normalized(staged.get(tag))
    }
    if differences:
        listed = "; ".join(
            f"{tag} {pair['scf']!r} in {scf_incar.name}, {pair['staged']!r} here" for tag, pair in differences.items()
        )
        a.fail(
            "vasp.charge_density_incompatible",
            f"this calculation cannot read the supplied charge density: {listed}",
            details={"differences": differences},
        )
        return False
    return True


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
        variants = parse_potcar_variants(
            a.parameter("potcar_variants", None) or a.setting("vasp.potcar_variants", None)
        )
        parse_diagnostic_policy(a.parameter("diagnostic_policy", None) or a.setting("vasp.diagnostic_policy", None))
        read_checker_specs(a.parameter("checkers", ""), base=a.payload)
    except ValueError as exception:
        a.fail("vasp.input_invalid", str(exception))
        return
    explicit = dict(a.parameter("incar_tags", {}) or {})
    shutil.copyfile(poscar, a.workdir / "POSCAR")
    override = a.payload / "files" / "INCAR.static"
    source = override if override.is_file() else TEMPLATES / "INCAR.static"
    if not source.is_file():
        a.fail(
            "vasp.template_missing",
            f"the INCAR template {source} is missing: reference this workflow as its package directory "
            "(a git URI ending in #vasp-static, or --workflow-dir), not as its runner file alone",
            details={"expected": str(source)},
        )
        return
    shutil.copyfile(source, incar)
    template = read_incar(incar)
    library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
        if variants:
            a.log.append("note", "potcar_variants is not used: this payload supplies its own POTCAR")
        library = None
    for name in ("CHGCAR", "WAVECAR"):
        supplied = a.payload / a.parameter(name.lower(), f"files/{name}")
        if supplied.is_file():
            shutil.copyfile(supplied, a.workdir / name)
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(a.parameter("kpoint_density", 20.0) or 20.0),
        centering=a.parameter("centering", VaspPreparationOptions.centering),
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        pseudopotential_library=library,
        potcar_variants=variants,
        parallel_tag=a.parameter("parallel_tag", None) or None,
        parallel_value=None if parallel_value is None else int(parallel_value),
        # The template's own EDIFF/EDIFFG are kept rather than derived; the job's tags win over both.
        incar_tags={
            **{tag: template[tag] for tag in ("EDIFF", "EDIFFG") if tag in template},
            **explicit,
        },
    )
    prepare_vasp_inputs(options, directory=a.workdir)
    if grid:
        write_automatic_kpoints(grid, a.workdir / "KPOINTS", centering=options.centering)
        a.log.append("note", f"fixed the k-point mesh at {' '.join(map(str, grid))} ({options.centering})")
    if not stage_charge_density(a, explicit):
        return
    a.log.append("note", "prepared a vasp.static calculation")
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
        report = run_vasp(
            argv,
            directory=a.workdir,
            timeout=a.parameter("timeout", 86400.0),
            policy=a.parameter("diagnostic_policy", None) or a.setting("vasp.diagnostic_policy", None),
            checkers=read_checker_specs(a.parameter("checkers", ""), base=a.payload),
        )
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
            excluded_operations=excluded_remedies(a.parameter("excluded_remedies", "")) | pinned_icharg(a),
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
    amplitude = float(a.parameter("rattle_amplitude", 0.0) or 0.0)
    if amplitude > 0:
        # The attempt is the entropy: reproducible, and different for every attempt.
        entropy = f"{a.context.job_key}:{a.context.attempt_ordinal}"
        rattle_poscar(a.workdir / "POSCAR", amplitude=amplitude, entropy=entropy)
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
