#!/usr/bin/env python3
"""One VASP relaxation: prepare inputs, run with remedies, publish the result.

The three steps are the whole workflow, spelled out on the ``httk.workflow.codes.vasp``
primitives so that this file is a starting point to copy and edit. ``prepare``
stages the structure and the INCAR of the job payload into the workdir and
derives everything else; ``run`` executes VASP under supervision and, when the
run fails in a way the reviewed remedy ladder recognizes, applies exactly one
remedy and asks for another attempt; ``publish`` publishes the files that
describe the finished calculation. ``vasp-relax-bash`` is the same workflow in
Bash, step for step.

The job inputs and parameters are documented in this repository's README. Nothing
here imports anything but an installed *httk-workflow*, so this one file is the
whole runner: it is the entry of this package directory, reference it by git URI
(e.g. ``git+https://github.com/httk/workflows-vasp#vasp-relax``) or with
``--workflow-dir``, publish it to a workspace runner store, or copy it and edit it.
"""

import shlex
import shutil

from httk.workflow import Attempt, Runner
from httk.workflow.codes.vasp import (
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    rattle_poscar,
    run_vasp,
    validate_vasp_workdir,
)

COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# Kept across a remedied rerun: they make the rerun cheaper, and VASP overwrites
# them itself when it reuses them.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")

run = Runner("vasp.relax")


@run.step
def prepare(a: Attempt) -> None:
    """Stage the payload inputs, derive the rest, and go on to run VASP."""

    validate_vasp_workdir(a.workdir)
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail(
            "vasp.input_missing",
            f"the starting structure {poscar.name} is not in this payload",
            details={"expected": str(poscar)},
        )
        return
    shutil.copyfile(poscar, a.workdir / "POSCAR")
    incar = a.payload / a.parameter("incar", "files/INCAR")
    if incar.is_file():
        shutil.copyfile(incar, a.workdir / "INCAR")
    else:
        # Everything an INCAR needs is derived below.
        (a.workdir / "INCAR").write_text("", encoding="utf-8")
    library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
        library = None
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(a.parameter("kpoint_density", 20.0) or 20.0),
        centering=a.parameter("centering", VaspPreparationOptions.centering),
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        pseudopotential_library=library,
        parallel_tag=a.parameter("parallel_tag", None) or None,
        parallel_value=None if parallel_value is None else int(parallel_value),
        incar_tags=dict(a.parameter("incar_tags", {}) or {}),
    )
    prepare_vasp_inputs(options, directory=a.workdir)
    a.log.append("note", "prepared a vasp.relax calculation")
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
    transactional = a.context.data_generation is not None
    published = []
    for name in (a.parameter("collect", COLLECT) or "").split():
        if (a.workdir / name).is_file():
            if transactional:
                a.put(a.workdir / name, f"{prefix}/{name}")
            published.append(name)
    if transactional:
        a.log.append("note", f"published to data/{prefix}: {', '.join(published) or 'nothing'}")
    else:
        # Without transactional data the persistent workdir is the result.
        a.log.append("note", f"kept in the workdir: {', '.join(published) or 'nothing'}")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
