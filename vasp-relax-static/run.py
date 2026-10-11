#!/usr/bin/env python3
"""The staged robust relaxation followed by one static calculation, in one job.

The relaxation is ``vasp-relax``'s, deliberately repeated here so that this file
and its package directory are a starting point to copy and edit on their own: the *httk* v1
``vasp-robust-relax-formenrg`` template (version 1.0 by Christopher Tholander),
spelled out on the ``httk.codes.vasp`` primitives. A structure far from
equilibrium is relaxed in stages, each from its own INCAR template in
``templates/``:

* ``preprerelax`` (``INCAR.preprerelax``) — a cheap, coarse relaxation that may
  also change the cell shape (``ISIF = 7``), at ``prerelax_kpoint_density``;
* ``prerelax`` (``INCAR.prerelax``) — a loose relaxation of positions and cell,
  at ``prerelax_kpoint_density``, which writes the charge density;
* ``relax1``, ``relax2``, … (``INCAR.relax``) — accurate relaxation rounds at
  ``kpoint_density``; from ``relax2`` on, each reads the previous round's charge
  density. From
  ``relax2`` on, a round whose final energy is within ten times ``EDIFFG`` of
  the previous round's ends the relaxation; otherwise, or when the round ran out
  of ``NSW`` ionic steps, another round follows. A round numbered above
  ``maximum_relax_rounds`` that still has not converged fails the job with
  ``vasp.relax_not_converged`` (v1 gave up at ``RELAXSTEP > 5``);
* ``static`` (``INCAR.static``) — one single point of the relaxed structure at
  ``kpoint_density``, with the electronic settings of ``INCAR.relax``, reading
  the final relaxation round's charge density.

The pre-stages and ``relax1`` also accept a run that only ran out of ionic
steps. Each finished relaxation stage is archived under ``stages/<stage>/`` and
its CONTCAR becomes the next stage's POSCAR. ``promote`` archives the final
relaxation round under ``relax/``, and the static stage runs at the top of the
workdir, which is where ``collect.py`` reads both.

This is the *httk* v1 formation-energy template without its database half: a
formation energy is assembled at job_records time from the total energies of several
such jobs and the elemental references they are compared against, which is an
analysis over finished jobs rather than a step of one job. What this runner
guarantees is the pair of numbers that assembly needs — one relaxed structure and
one total energy of it, both with the evidence of how they were produced.

Edit the templates of a copied package to change the settings, or override one
stage's template wholesale in a single job with a payload file
``files/INCAR.<stage>`` (e.g. ``files/INCAR.static``). A template's own ``EDIFF``
and ``EDIFFG`` are kept; when absent they are derived from ``accuracy_per_atom``
with the v1 margin (500 for a relaxation, 33 otherwise). The ``incar_tags``
parameter applies to every stage and wins over templates and derived values.
A remedy applied in one stage stays applied in every later stage, as in v1.

A template's ``ICHARG = 1`` (read the charge density) stays only when the
workdir holds a non-empty CHGCAR written by a run with the same ``PREC``,
``ENCUT`` and ``ENAUG``, and becomes ``ICHARG = 2`` otherwise: a density is not
read onto a different FFT grid. So ``relax1`` starts from atomic charges, since
``prerelax`` ran at ``PREC = Normal`` and the default cutoff, and later relax
rounds reuse the density. v1's own check started every stage from atomic
charges, whatever the template said. All
stages share one job: one remedy budget, one remedy ladder, one workdir.

The steps: ``prepare`` stages the structure and POTCAR once; ``stage`` prepares
the current stage's inputs from its template; ``run`` executes VASP under
supervision and either finishes the stage or, when the run fails in a way the
reviewed remedy ladder recognizes, applies exactly one remedy and retries;
``promote`` hands the relaxed structure to the static stage; and ``publish``
publishes both the relaxation and the single point.

The job inputs and parameters are documented in this repository's README. Nothing
here imports anything but an installed *httk-workflow*, so this file and its
``templates/`` are the whole runner: it is the entry of this package directory,
reference it by git URI (e.g. ``git+https://github.com/httk/workflows-vasp#vasp-relax-static``)
or with ``--workflow-dir`` (installed into the workspace first), or copy the package
directory and edit it. The runner file alone, without its ``templates/``, fails
with ``vasp.template_missing``.
"""

import math
import shlex
import shutil
from pathlib import Path

from httk.codes.vasp import (
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    contcar_to_poscar,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    rattle_poscar,
    read_incar,
    reapply_vasp_remedies,
    run_vasp,
    update_incar,
    validate_vasp_workdir,
)
from httk.workflow import Attempt, Runner

COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# Kept across a remedied rerun and across stages: they make the next run cheaper,
# a relax stage reads the CHGCAR of the stage before it, and VASP overwrites them
# itself when it reuses them.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")
# What stages/<stage>/ keeps of every finished relaxation stage.
STAGE_ARCHIVE = (
    "POSCAR",
    "INCAR",
    "KPOINTS",
    "OUTCAR",
    "OSZICAR",
    "XDATCAR",
    "CONTCAR",
    "vasprun.xml",
    "vasp-run-report.json",
    "vasp.out",
    "vasp.err",
)
NEXT_STAGE = {"preprerelax": "prerelax", "prerelax": "relax", "relax": "relax"}
TEMPLATES = Path(__file__).resolve().parent / "templates"
# The INCAR tags that fix the FFT grid a CHGCAR is written on.
GRID_TAGS = ("PREC", "ENCUT", "ENAUG")
# Where the final relaxation round is archived inside the workdir before the static stage overwrites it.
RELAX_ARCHIVE = "relax"

run = Runner("vasp.relax-static")


def chgcar_grid(tags: dict[str, str]) -> str:
    """Return the grid a CHGCAR of these INCAR tags is on: PREC, ENCUT and ENAUG, lower-cased, absent ones empty."""

    return ";".join(tags.get(tag, "").lower() for tag in GRID_TAGS)


@run.step
def prepare(a: Attempt) -> None:
    """Stage the structure and POTCAR of the payload, and go on to the first stage."""

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
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
    a.advance("stage", state={"stage": "preprerelax", "relax_round": 0})


@run.step
def stage(a: Attempt) -> None:
    """Prepare the current stage's inputs from its INCAR template, and go on to run VASP."""

    name = str(a.state["stage"])
    label = f"relax{a.state['relax_round']}" if name == "relax" else name
    incar = a.workdir / "INCAR"
    override = a.payload / "files" / f"INCAR.{name}"
    source = override if override.is_file() else TEMPLATES / f"INCAR.{name}"
    if not source.is_file():
        a.fail(
            "vasp.template_missing",
            f"the INCAR template {source} is missing: reference this workflow as its package directory "
            "(a git URI ending in #vasp-relax-static, or --workflow-dir), not as its runner file alone",
            details={"expected": str(source)},
        )
        return
    shutil.copyfile(source, incar)
    template = read_incar(incar)
    if name in ("preprerelax", "prerelax"):
        density = a.parameter("prerelax_kpoint_density", 10.0)
    else:
        density = a.parameter("kpoint_density", 20.0)
    # Only the first stage assembles a POTCAR from the library; later stages find it.
    library = None
    if not (a.workdir / "POTCAR").is_file():
        library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(density or 20.0),
        centering=a.parameter("centering", VaspPreparationOptions.centering),
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        # v1: a relaxation converges its electrons tighter, so its ionic steps are more accurate.
        ediff_margin=500.0 if float(template.get("NSW", "0")) > 1 else 33.0,
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
    # Remedies are sticky across stages, as in v1: what an earlier stage needed, this one gets too.
    reapply_vasp_remedies(directory=a.workdir, history_path=job_remedy_history_path(a.payload))
    # The last run's charge density is only read on the grid it was written on:
    # ICHARG = 1 stays only with a CHGCAR from an accepted run with the same PREC,
    # ENCUT and ENAUG, as recorded in job state (not the workdir INCAR, which a
    # retried stage step has already replaced with its template).
    chgcar = a.workdir / "CHGCAR"
    current = read_incar(incar)
    reusable = chgcar.is_file() and chgcar.stat().st_size and a.state.get("chgcar_grid") == chgcar_grid(current)
    if current.get("ICHARG") == "1" and not reusable:
        update_incar({"ICHARG": 2}, incar)
    a.log.append("note", f"prepared the vasp.relax-static stage {label}")
    a.advance("run")


@run.step(name="run")
def run_step(a: Attempt) -> None:
    """Run VASP, finish the stage, remedy a recognized failure, or fail with what was diagnosed."""

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
    # Running out of NSW ionic steps is progress, not a failure, of a relaxation stage.
    ionic = report.classification == "nonconverged" and {
        item.code for item in report.diagnostics if item.severity != "warning"
    } == {"ionic_nonconvergence"}
    if report.classification == "completed" or ionic:
        finish_stage(a, state, energy, ionic=ionic)
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
    # The stage keeps its remedied inputs: the retry runs VASP again without re-preparing.
    a.retry(f"applied the {decision.policy} remedy for {decision.problem}")


def finish_stage(a: Attempt, state: dict[str, object], energy: float | None, *, ionic: bool) -> None:
    """Publish after the static stage, promote a converged relaxation, or archive the stage and go on."""

    ran = read_incar(a.workdir / "INCAR")
    # The grid of the CHGCAR this accepted run wrote, committed with its outcome so a
    # retried stage step still knows it; a run with LCHARG = .FALSE. wrote none.
    writes_chgcar = not ran.get("LCHARG", ".TRUE.").lstrip(".").upper().startswith("F")
    state["chgcar_grid"] = chgcar_grid(ran) if writes_chgcar else None
    name = str(a.state["stage"])
    if name == "static":
        a.advance("publish", state=state)
        return
    rounds = int(a.state.get("relax_round", 0))
    label = f"relax{rounds}" if name == "relax" else name
    if energy is not None:
        state[f"{label}_energy"] = energy
    if name == "relax" and rounds >= 2:
        tags = ran
        # VASP's own default EDIFFG is ten times EDIFF.
        limit = 10 * abs(float(tags.get("EDIFFG", 10 * float(tags.get("EDIFF", "1e-4")))))
        previous = a.state.get(f"relax{rounds - 1}_energy")
        change = math.inf if energy is None or previous is None else abs(energy - float(str(previous)))
        if not ionic and change <= limit:
            a.advance("promote", state=state)
            return
        if rounds > int(a.parameter("maximum_relax_rounds", 5)):
            a.state.merge(state)
            reason = "ran out of NSW ionic steps" if ionic else f"changed the energy by {change:g} eV (limit {limit:g})"
            a.fail(
                "vasp.relax_not_converged",
                f"the relaxation did not converge in rounds relax1 to relax{rounds}: relax{rounds} {reason}",
                details={"rounds": rounds, "energy_change": None if ionic else change, "limit": limit},
            )
            return
    contcar = a.workdir / "CONTCAR"
    if not contcar.is_file() or not contcar.read_text(encoding="utf-8", errors="replace").strip():
        a.state.merge(state)
        a.fail("vasp.no_relaxed_structure", f"the stage {label} left no CONTCAR to continue from")
        return
    archive = a.workdir / "stages" / label
    archive.mkdir(parents=True, exist_ok=True)
    for file in STAGE_ARCHIVE:
        if (a.workdir / file).is_file():
            shutil.copyfile(a.workdir / file, archive / file)
    # The reference POSCAR only lends its comment line, which carries the MAGMOM override.
    contcar_to_poscar(contcar, reference=a.workdir / "POSCAR", output=a.workdir / "POSCAR")
    state["stage"] = NEXT_STAGE[name]
    state["relax_round"] = rounds + 1 if state["stage"] == "relax" else 0
    a.log.append("note", f"finished the stage {label}" + (" out of ionic steps" if ionic else ""))
    a.advance("stage", state=state)


@run.step
def promote(a: Attempt) -> None:
    """Archive the converged relaxation, adopt its CONTCAR, and go on to the static stage."""

    contcar = a.workdir / "CONTCAR"
    if not contcar.is_file() or not contcar.read_text(encoding="utf-8", errors="replace").strip():
        a.fail("vasp.no_relaxed_structure", "the relaxation completed without leaving a CONTCAR to run statically")
        return
    archive = a.workdir / RELAX_ARCHIVE
    archive.mkdir(exist_ok=True)
    for name in (a.parameter("collect", COLLECT) or "").split():
        if (a.workdir / name).is_file():
            shutil.copyfile(a.workdir / name, archive / name)
    # The relaxed cell is the structure of the single point; the reference POSCAR
    # only lends it its comment line, which carries the MAGMOM override.
    contcar_to_poscar(contcar, reference=a.workdir / "POSCAR", output=a.workdir / "POSCAR")
    a.log.append("note", "promoted the relaxed structure to a static calculation")
    a.advance(
        "stage",
        state={
            "stage": "static",
            "relax_energy": a.state.get("energy"),
            "relax_classification": a.state.get("classification"),
        },
    )


@run.step
def publish(a: Attempt) -> None:
    """Publish both the relaxation and the single point, and complete the job."""

    base = a.parameter("data_prefix", "") or ""
    # Results stay in the persistent workdir; publish_data opts a job into a curated data/ copy too.
    to_data = a.parameter("publish_data", False) is True
    for stage_name, directory in (("static", a.workdir), (RELAX_ARCHIVE, a.workdir / RELAX_ARCHIVE)):
        prefix = f"{base}/{stage_name}" if base else stage_name
        published = []
        for name in (a.parameter("collect", COLLECT) or "").split():
            if (directory / name).is_file():
                if to_data:
                    a.put(directory / name, f"{prefix}/{name}")
                published.append(name)
        if to_data:
            a.log.append("note", f"published to data/{prefix}: {', '.join(published) or 'nothing'}")
        else:
            # Without publish_data the persistent workdir is the result.
            a.log.append("note", f"kept in the workdir: {', '.join(published) or 'nothing'}")
    stages = sorted(path.name for path in (a.workdir / "stages").iterdir()) if (a.workdir / "stages").is_dir() else []
    a.log.append("note", f"earlier relaxation stages kept in the workdir under stages/: {', '.join(stages) or 'none'}")
    a.state["static_energy"] = a.state.get("energy")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
