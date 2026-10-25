#!/usr/bin/env python3
"""Any VASP calculation, from the input files the job supplies.

The job supplies every input file but the POTCAR — INCAR, POSCAR, and whatever
else the calculation reads (KPOINTS, WAVECAR, CHGCAR, ...) — in its payload's
``files/``, and this workflow runs VASP on them as given: nothing is derived,
except that the VASP launch wrapper may write planned ``KPAR``/``NCORE`` and a
rounded ``NBANDS`` into the INCAR (and logs them), unless the INCAR sets them
(``KPAR``, ``NCORE`` or ``NPAR`` for the planned pair) or the settings
``vasp.parallelization`` and ``vasp.nbands_rounding`` are ``off``. With both
``off`` (or on one rank) the INCAR reaches VASP byte for byte, as KPOINTS always
does. The steps, spelled out
on the ``httk.codes.vasp`` primitives so that this file is a starting point to
copy and edit: ``prepare`` copies every regular file directly under ``files/``
into the workdir — except ``vasp-run-report.json`` and ``POTCAR.provenance.json``,
which are this runner's own bookkeeping — and supplies a POTCAR when the job did not; ``run`` executes
VASP under supervision; and ``publish`` publishes the result.

The structure may come as the ``structure`` input or as ``files/POSCAR``; both
are the same payload file, ``files/POSCAR``, so a job gives one of them. A
POTCAR is taken from ``files/POTCAR``, or else assembled from the
pseudopotential library (the ``vasp.pseudo_library`` setting or the
``pseudopotential_library`` parameter). KPOINTS is optional, since VASP can
derive the k-points from ``KSPACING``.

Remedies are off by default: a run that does not complete fails the job at once,
with the diagnosis and the remedy the reviewed ladder would have applied in the
failure details. A job that sets ``maximum_remedies`` above zero opts into
applying up to that many remedies, each followed by another run.

The job inputs and parameters are documented in this repository's README. Nothing
here imports anything but an installed *httk-workflow*, so this one file is the
whole runner: it is the entry of this package directory, reference it by git URI
(e.g. ``git+https://github.com/httk/workflows-vasp#vasp-general``) or with
``--workflow-dir`` (installed into the workspace first), or copy it and edit it.
"""

import shlex
import shutil

from httk.codes.vasp import (
    apply_vasp_remedy,
    assemble_potcar,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    rattle_poscar,
    run_vasp,
    validate_vasp_workdir,
)
from httk.workflow import Attempt, Runner

COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# The runner's own bookkeeping, written by this calculation: never taken from the job's files.
BOOKKEEPING = ("vasp-run-report.json", "POTCAR.provenance.json")
# Kept across a remedied rerun: they make the rerun cheaper, and VASP overwrites
# them itself when it reuses them.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")

run = Runner("vasp.general")


@run.step
def prepare(a: Attempt) -> None:
    """Copy the job's input files into the workdir, supply a POTCAR, and go on to run VASP."""

    validate_vasp_workdir(a.workdir)
    files = a.payload / "files"
    skipped = []
    if files.is_dir():
        for path in sorted(files.iterdir()):
            if path.name in BOOKKEEPING:
                skipped.append(path.name)
            elif path.is_file():
                shutil.copyfile(path, a.workdir / path.name)
    if skipped:
        a.log.append("note", f"skipped the runner's own bookkeeping files in files/: {', '.join(skipped)}")
    missing = [name for name in ("INCAR", "POSCAR") if not (a.workdir / name).is_file()]
    if missing:
        a.fail(
            "vasp.input_missing",
            f"{' and '.join(missing)} {'is' if len(missing) == 1 else 'are'} not in this payload's files/",
            details={"missing": missing},
        )
        return
    if not (a.workdir / "POTCAR").is_file():
        library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
        if library is None:
            a.fail(
                "vasp.input_missing",
                "no files/POTCAR in this payload and no pseudopotential library to assemble one: "
                "set the vasp.pseudo_library workspace setting or the pseudopotential_library parameter",
                details={"missing": ["POTCAR"]},
            )
            return
        try:
            assemble_potcar(library, poscar=a.workdir / "POSCAR", output=a.workdir / "POTCAR")
        except FileNotFoundError as exception:
            a.fail("vasp.input_missing", f"could not assemble a POTCAR: {exception}", details={"missing": ["POTCAR"]})
            return
    a.log.append("note", "prepared a vasp.general calculation from the job's own input files")
    a.advance("run")


@run.step(name="run")
def run_step(a: Attempt) -> None:
    """Run VASP, remedy a recognized failure when the job allows it, or fail with what was diagnosed."""

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
    # Off by default: the job's inputs are the job's own, so nothing edits them unasked.
    maximum = int(a.parameter("maximum_remedies", 0))
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
        if maximum == 0:
            message = f"VASP {report.classification}; remedies are off (maximum_remedies = 0), so nothing was retried"
        elif applied >= maximum:
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
