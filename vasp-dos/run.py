#!/usr/bin/env python3
"""The electronic density of states and band gaps of a structure that is already chosen (normally relaxed).

The workflow is spelled out on the ``httk.codes.vasp`` primitives so that this
package directory is a starting point to copy and edit. A self-consistent run on
a uniform mesh writes the charge density, and a non-self-consistent run
(``ICHARG = 11``) reads it on a dense uniform mesh with the tetrahedron method,
writing the DOSCAR and the EIGENVAL the band gaps come from. Both runs compute
the input cell as given.

The steps, each with its stage fixed by the step itself, so that rerunning any
step after a crash repeats only that step's own work:

* ``prepare`` stages the structure and checks the inputs; when the payload gives
  a charge density (``files/CHGCAR`` with the ``files/scf-INCAR`` it was computed
  with, optionally ``files/scf-OUTCAR``), it goes straight to the nscf stage;
* ``stage_scf`` and ``run_scf`` prepare and run the self-consistent stage from
  ``templates/INCAR.scf``, on the automatic mesh of ``kpoint_density``, and
  archive it under ``scf/``, leaving its CHGCAR at the top of the workdir;
* ``stage_nscf`` and ``run_nscf`` prepare and run the non-self-consistent stage
  from ``templates/INCAR.nscf`` on the automatic mesh of ``dos_kpoint_density``;
* ``analyse`` finds the band edges from the EIGENVAL occupations, relative to
  the DOSCAR Fermi energy, and writes ``band_edges.json``;
* ``publish`` publishes both stages.

Edit the templates of a copied package to change the settings, or replace one
in a single job with a payload file ``files/INCAR.scf`` or ``files/INCAR.nscf``.
The template's own ``EDIFF`` and ``EDIFFG`` are kept; when absent they are
derived from ``accuracy_per_atom``. The ``incar_tags`` parameter applies to both
stages and wins over the templates and every derived value. A remedy applied in
the scf stage stays applied in the nscf stage, as in the other staged runners.

The guards, each failing the job by name. Before any VASP run, ``prepare``
checks ``incar_tags``, each stage's payload INCAR or template, and a given
``files/scf-INCAR``: any ``ICHARG`` in ``incar_tags`` or an nscf INCAR's
``ICHARG`` other than 11 is ``vasp.invalid_dos_inputs``, and a ``files/CHGCAR``
without its ``files/scf-INCAR`` is ``vasp.input_missing``. The nscf stage checks
again what VASP is about to run: one whose ``ICHARG`` is not 11 is
``vasp.invalid_dos_inputs``; one without a CHGCAR is
``vasp.charge_density_missing``; and one whose density-defining tags differ from
the scf run's is ``vasp.charge_density_incompatible``. Remedies that move the
structure (``scale_lattice``, ``contcar_to_poscar``) are never applied, nor, in
the nscf stage, an ``ICHARG`` change. The nscf stage always runs with
``ICHARG = 11`` and ``LCHARG = .FALSE.``, also when an scf remedy set
``ICHARG = 2``.

The scope is collinear LDA/GGA(+U): ``LHFCALC``, ``METAGGA``, ``LNONCOLLINEAR``
and ``LSORBIT`` fail with ``vasp.unsupported_calculation``, in ``prepare`` and
again on each stage's prepared INCAR. Hybrid functionals need zero-weight
explicit k-points in a self-consistent run, and meta-GGAs the kinetic-energy
density, which a non-self-consistent ``ICHARG = 11`` run does not read.

The job inputs and parameters are documented in this repository's README. The
``analyse`` step imports an installed *httk-atomistic* and *httk-analyse*
(``vasp.analysis_unavailable`` without the latter). This file and its
``templates/`` are the whole runner: it is the entry of this package directory,
reference it by git URI (e.g.
``git+https://github.com/httk/workflows-vasp#vasp-dos``) or with
``--workflow-dir`` (installed into the workspace first), or copy the package
directory and edit it. The runner file alone, without its ``templates/``, fails
with ``vasp.template_missing``.
"""

import dataclasses
import json
import os
import shlex
import shutil
from pathlib import Path

from httk.codes.vasp import (
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    read_incar,
    reapply_vasp_remedies,
    run_vasp,
    update_incar,
    validate_vasp_workdir,
)
from httk.workflow import Attempt, Runner

COLLECT = (
    "INCAR KPOINTS POSCAR OUTCAR OSZICAR EIGENVAL DOSCAR vasprun.xml band_edges.json vasp-run-report.json "
    "POTCAR.provenance.json"
)
# Kept across a remedied rerun and between the stages: the nscf stage reads the scf CHGCAR.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")
# What scf/ always keeps of the self-consistent run, besides the collect files.
SCF_ARCHIVE = ("INCAR", "KPOINTS", "POSCAR", "OUTCAR", "OSZICAR", "vasprun.xml", "vasp-run-report.json")
# What publish_data always publishes per stage, besides the collect files: what collect.py reads.
REQUIRED = {"scf": ("OUTCAR",), "nscf": ("EIGENVAL", "DOSCAR", "band_edges.json")}
# Tags outside the collinear LDA/GGA(+U) scope; METAGGA in any value, the others when true.
SCOPE_TAGS = ("LHFCALC", "LNONCOLLINEAR", "LSORBIT", "METAGGA")
# The INCAR tags that define the charge density the nscf stage reads.
DENSITY_TAGS = (
    *("PREC", "ENCUT", "ENAUG", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF", "ISPIN", "MAGMOM", "LMAXMIX"),
    *("GGA", "METAGGA", "LASPH", "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "NELECT", "LNONCOLLINEAR", "LSORBIT"),
    *("GGA_COMPAT", "LUSE_VDW", "AGGAC", "AGGAX", "ALDAC", "PARAM1", "PARAM2", "LMAXPAW", "XC"),
)
# Remedies never applied: the structure is the input's, and the nscf stage must stay non-self-consistent.
EXCLUDED = {
    "scf": {"scale_lattice", "contcar_to_poscar"},
    "nscf": {"scale_lattice", "contcar_to_poscar", "incar.ICHARG"},
}
# The k-point density of each stage: its parameter and default.
DENSITY = {"scf": ("kpoint_density", 20.0), "nscf": ("dos_kpoint_density", 30.0)}
TEMPLATES = Path(__file__).resolve().parent / "templates"

run = Runner("vasp.dos")


def is_true(value: str | None) -> bool:
    """Return whether an INCAR boolean reads true (``.TRUE.``, ``T``, ``true``, ...)."""

    return value is not None and value.strip().lstrip(".").upper().startswith("T")


def out_of_scope(tags: dict[str, str]) -> list[str]:
    """Return the tags of *tags* outside the collinear LDA/GGA(+U) scope."""

    return [tag for tag in SCOPE_TAGS if tag in tags and (tag == "METAGGA" or is_true(tags[tag]))]


@run.step
def prepare(a: Attempt) -> None:
    """Stage the structure, check the inputs, and go on to the scf stage, or to the nscf stage on a given density."""

    validate_vasp_workdir(a.workdir)
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail(
            "vasp.input_missing",
            f"the starting structure {poscar.name} is not in this payload",
            details={"expected": str(poscar)},
        )
        return
    files = a.payload / "files"
    given = (files / "CHGCAR").is_file()
    if given and not (files / "scf-INCAR").is_file():
        a.fail(
            "vasp.input_missing",
            "files/CHGCAR is given without files/scf-INCAR, the INCAR it was computed with: "
            "its density-defining tags must be known to check that the nscf stage can read it",
            details={"expected": str(files / "scf-INCAR")},
        )
        return
    # Checked here, before any VASP run, on every source of either stage's INCAR: the job's
    # incar_tags (applied to both stages), each stage's payload override or template, and a
    # given density's INCAR.
    explicit = {
        str(tag).strip().upper(): str(value) for tag, value in dict(a.parameter("incar_tags", {}) or {}).items()
    }
    sources = {"incar_tags": explicit}
    for name in ("scf", "nscf"):
        override = files / f"INCAR.{name}"
        source = override if override.is_file() else TEMPLATES / f"INCAR.{name}"
        # A missing template fails its stage with vasp.template_missing.
        if source.is_file():
            sources[f"{source.parent.name}/INCAR.{name}"] = read_incar(source)
    if given:
        sources["files/scf-INCAR"] = read_incar(files / "scf-INCAR")
    # incar_tags apply to both stages, so any ICHARG there is wrong for one of them.
    for source, tags in sources.items():
        icharg = tags.get("ICHARG")
        if icharg is not None and (source == "incar_tags" or source.endswith("INCAR.nscf") and icharg.strip() != "11"):
            a.fail(
                "vasp.invalid_dos_inputs",
                f"ICHARG = {icharg} from {source}: the scf stage starts from atomic charges and the nscf stage "
                "needs ICHARG = 11",
                details={"source": source, "ICHARG": icharg},
            )
            return
    for source, tags in sources.items():
        unsupported = out_of_scope(tags)
        if unsupported:
            a.fail(
                "vasp.unsupported_calculation",
                f"{source} sets {', '.join(unsupported)}: this workflow computes collinear LDA/GGA(+U) "
                "densities of states only",
                details={"source": source, "tags": {tag: tags[tag] for tag in unsupported}},
            )
            return
    shutil.copyfile(poscar, a.workdir / "POSCAR")
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
    if not given:
        a.advance("stage_scf", state={"given_density": False})
        return
    shutil.copyfile(files / "CHGCAR", a.workdir / "CHGCAR")
    if (files / "scf-OUTCAR").is_file():
        # For provenance only: the energy belongs to the job that computed the density.
        (a.workdir / "scf").mkdir(exist_ok=True)
        shutil.copyfile(files / "scf-OUTCAR", a.workdir / "scf" / "OUTCAR")
    a.log.append("note", "reading the given charge density files/CHGCAR: no scf stage runs")
    a.advance("stage_nscf", state={"given_density": True, "scf_incar": sources["files/scf-INCAR"]})


def stage_inputs(a: Attempt, name: str) -> bool:
    """Prepare one stage's inputs from its INCAR template; return whether they are supported."""

    incar = a.workdir / "INCAR"
    override = a.payload / "files" / f"INCAR.{name}"
    source = override if override.is_file() else TEMPLATES / f"INCAR.{name}"
    if not source.is_file():
        a.fail(
            "vasp.template_missing",
            f"the INCAR template {source} is missing: reference this workflow as its package directory "
            "(a git URI ending in #vasp-dos, or --workflow-dir), not as its runner file alone",
            details={"expected": str(source)},
        )
        return False
    shutil.copyfile(source, incar)
    template = read_incar(incar)
    # Only the first stage that runs assembles a POTCAR from the library; a later one finds it.
    library = None
    if not (a.workdir / "POTCAR").is_file():
        library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    parallel_value = a.parameter("parallel_value", None)
    parameter, default = DENSITY[name]
    options = VaspPreparationOptions(
        kpoint_density=float(a.parameter(parameter, default) or default),
        # Gamma-centred meshes for the tetrahedron method, as in vasp-relax-static-coarse.
        centering=a.parameter("centering", "Gamma"),
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        ediff_margin=33.0,
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
    reapply_vasp_remedies(directory=a.workdir, history_path=job_remedy_history_path(a.payload))
    tags = read_incar(incar)
    # Authoritative over prepare's early check: this is the INCAR VASP will run with.
    unsupported = out_of_scope(tags)
    if unsupported:
        a.fail(
            "vasp.unsupported_calculation",
            f"the {name} stage sets {', '.join(unsupported)}: this workflow computes collinear LDA/GGA(+U) "
            "densities of states only",
            details={"stage": name, "tags": {tag: tags[tag] for tag in unsupported}},
        )
        return False
    return True


@run.step
def stage_scf(a: Attempt) -> None:
    """Prepare the self-consistent stage on the uniform mesh, and go on to run it."""

    if stage_inputs(a, "scf"):
        a.log.append("note", "prepared the vasp.dos stage scf")
        a.advance("run_scf")


@run.step
def stage_nscf(a: Attempt) -> None:
    """Prepare the non-self-consistent stage on the dense uniform mesh, and go on to run it."""

    if not stage_inputs(a, "nscf"):
        return
    # A replayed scf remedy may have set ICHARG = 2; the scf CHGCAR must not be overwritten.
    update_incar({"ICHARG": 11, "LCHARG": ".FALSE."}, a.workdir / "INCAR")
    a.log.append("note", "prepared the vasp.dos stage nscf")
    a.advance("run_nscf")


def nscf_inputs_problem(a: Attempt) -> tuple[str, str, dict[str, object]] | None:
    """Return the failure that forbids launching the nscf stage, or ``None`` when it may run."""

    tags = read_incar(a.workdir / "INCAR")
    if tags.get("ICHARG", "").strip() != "11":
        return (
            "vasp.invalid_dos_inputs",
            f"the nscf stage INCAR has ICHARG = {tags.get('ICHARG')}, not 11",
            {"ICHARG": tags.get("ICHARG")},
        )
    chgcar = a.workdir / "CHGCAR"
    if not chgcar.is_file() or not chgcar.stat().st_size:
        return ("vasp.charge_density_missing", "there is no CHGCAR for the nscf stage to read", {})
    scf = dict(a.state.get("scf_incar") or {})

    def normalized(value: object) -> str | None:
        return None if value is None else str(value).strip().lower()

    differences = {
        tag: {"scf": scf.get(tag), "nscf": tags.get(tag)}
        for tag in DENSITY_TAGS
        if normalized(scf.get(tag)) != normalized(tags.get(tag))
    }
    if differences:
        listed = "; ".join(
            f"{tag} {pair['scf']!r} in scf, {pair['nscf']!r} in nscf" for tag, pair in differences.items()
        )
        return (
            "vasp.charge_density_incompatible",
            f"the nscf stage cannot read the scf charge density: {listed}",
            {"differences": differences},
        )
    return None


def run_stage(a: Attempt, name: str) -> None:
    """Run VASP for one stage, finish it, remedy a recognized failure, or fail with what was diagnosed."""

    # vasp.command resolves the job parameter, HTTK_VASP_COMMAND, then the workspace setting.
    argv = shlex.split(a.setting("vasp.command", a.parameter("vasp_command", None)) or "")
    if not argv:
        a.fail(
            "vasp.command_missing",
            "no VASP command is configured: set it with `httk workspace settings set --key vasp.command --value '...' WORKSPACE`, "
            "or set HTTK_VASP_COMMAND on the machine that runs this job, or give the job a vasp_command parameter",
        )
        return
    # A rerun must not read the previous run's outputs; the CHGCAR the nscf stage reads is kept.
    clean_vasp_outputs(a.workdir, keep=KEEP_BETWEEN_RUNS)
    if name == "nscf":
        # Checked before every launch, retries included.
        problem = nscf_inputs_problem(a)
        if problem is not None:
            code, message, details = problem
            a.fail(code, message, details=details)
            return
    try:
        report = run_vasp(argv, directory=a.workdir, timeout=a.parameter("timeout", 86400.0))
    except OSError as exception:
        a.fail("vasp.failed", f"could not run VASP at all: {exception}")
        return
    a.log.append("note", f"VASP {name} {report.classification}")
    if report.classification == "completed":
        if name == "nscf":
            a.advance("analyse", state={"nscf_classification": report.classification})
            return
        oszicar = a.workdir / "OSZICAR"
        state = {
            "scf_energy": last_oszicar_energy(oszicar) if oszicar.is_file() else None,
            "scf_incar": read_incar(a.workdir / "INCAR"),
            "scf_classification": report.classification,
        }
        archive = a.workdir / "scf"
        archive.mkdir(exist_ok=True)
        for file in dict.fromkeys((*SCF_ARCHIVE, *(a.parameter("collect", COLLECT) or "").split())):
            if (a.workdir / file).is_file():
                shutil.copyfile(a.workdir / file, archive / file)
        a.advance("stage_nscf", state=state)
        return

    state: dict[str, object] = {"classification": report.classification}
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
            excluded_operations=EXCLUDED[name],
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
    # The stage keeps its remedied inputs: the retry runs VASP again without re-preparing.
    a.retry(f"applied the {decision.policy} remedy for {decision.problem}")


@run.step
def run_scf(a: Attempt) -> None:
    """Run the self-consistent stage, archive it under ``scf/``, and go on to the nscf stage."""

    run_stage(a, "scf")


@run.step
def run_nscf(a: Attempt) -> None:
    """Run the non-self-consistent stage, and go on to analyse it."""

    run_stage(a, "nscf")


@run.step
def analyse(a: Attempt) -> None:
    """Find the band edges of the nscf run, write ``band_edges.json``, and go on to publish."""

    try:
        from httk.analyse.matsci.electronic import band_edges
    except ImportError as exception:
        a.fail(
            "vasp.analysis_unavailable",
            f"the band edges need httk-analyse: {exception} (pip install httk-analyse)",
        )
        return
    from httk.atomistic.integrations.vasp.io import read_eigenval
    from httk.atomistic.integrations.vasp.io.doscar import read_doscar

    eigenval = read_eigenval(a.workdir / "EIGENVAL")
    occupations = eigenval.occupations_floats()
    if occupations is None:
        a.fail("vasp.eigenvalues_without_occupations", "the nscf EIGENVAL has no occupation columns")
        return
    tolerance = float(a.parameter("occupation_tolerance", 1e-3))
    energies = eigenval.energies_floats()
    # The Bloechl-corrected tetrahedron method (ISMEAR = -5) gives states at the Fermi level of a
    # metal occupations slightly outside [0, 1] (1.0012, -0.0007); only a partially occupied state
    # gets one, so it counts as partial (0.5). Insulators get exact 0 and 1.
    out_of_range = sum(not -tolerance <= o <= 1 + tolerance for spin in occupations for k in spin for o in k)
    occupations = [[[o if -tolerance <= o <= 1 + tolerance else 0.5 for o in k] for k in spin] for spin in occupations]

    def edges(spins: range) -> dict[str, object]:
        # band_edges takes (bands, kpoints); spin channels stack along the band axis, up first.
        stacked = [[energies[s][k][b] for k in range(eigenval.nkpoints)] for s in spins for b in range(eigenval.nbands)]
        held = [[occupations[s][k][b] for k in range(eigenval.nkpoints)] for s in spins for b in range(eigenval.nbands)]
        result = band_edges(stacked, held, maximum_occupation=1, occupation_tolerance=tolerance, energy_reference=fermi)
        return dataclasses.asdict(result)

    def channel(spin: int) -> dict[str, object] | None:
        # A channel without both occupied and empty states (say, an empty minority channel) has no edges.
        try:
            return edges(range(spin, spin + 1))
        except ValueError:
            return None

    try:
        doscar = read_doscar(a.workdir / "DOSCAR", spin_mode="collinear" if eigenval.ispin == 2 else "nonmagnetic")
        fermi = float(doscar.fermi_energy)
        document = edges(range(eigenval.ispin))
        per_spin = [channel(spin) for spin in range(eigenval.ispin)] if eigenval.ispin == 2 else []
        direct = [edge["direct_gap"] for edge in per_spin if edge is not None]
        if per_spin and not document["metallic"]:
            if not direct:
                raise ValueError("neither spin channel has both occupied and empty states")
            # Spin-conserving transitions only: the smallest same-k gap within one spin channel.
            document["direct_gap"] = min(direct)
    except ValueError as exception:
        a.fail("vasp.band_edges_failed", f"cannot find the band edges of the nscf run: {exception}")
        return
    document |= {
        "ispin": eigenval.ispin,
        "fermi_energy": fermi,
        "occupation_tolerance": tolerance,
        "occupations_out_of_range": out_of_range,
        "nkpoints": eigenval.nkpoints,
        "nbands": eigenval.nbands,
        "given_density": bool(a.state.get("given_density", False)),
    }
    if per_spin:
        document["per_spin"] = per_spin
    path = a.workdir / "band_edges.json"
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, path)
    metallic = document["metallic"]
    band_gap = 0.0 if metallic else document["indirect_gap"]
    direct_band_gap = 0.0 if metallic else document["direct_gap"]
    a.state.merge(
        {"metallic": metallic, "band_gap": band_gap, "direct_band_gap": direct_band_gap, "fermi_energy": fermi}
    )
    found = "metallic" if metallic else f"a band gap of {band_gap:.4f} eV (direct {direct_band_gap:.4f} eV)"
    a.log.append("note", f"found {found}, Fermi energy {fermi:.4f} eV")
    a.advance("publish")


@run.step
def publish(a: Attempt) -> None:
    """Publish both stages, and complete the job."""

    base = a.parameter("data_prefix", "") or ""
    # Results stay in the persistent workdir; publish_data opts a job into a curated data/ copy too.
    to_data = a.parameter("publish_data", False) is True
    collect = (a.parameter("collect", COLLECT) or "").split()
    for stage_name, directory in (("scf", a.workdir / "scf"), ("nscf", a.workdir)):
        prefix = f"{base}/{stage_name}" if base else stage_name
        published = []
        for name in dict.fromkeys((*REQUIRED[stage_name], *collect)):
            if (directory / name).is_file():
                if to_data:
                    a.put(directory / name, f"{prefix}/{name}")
                published.append(name)
        if to_data:
            a.log.append("note", f"published to data/{prefix}: {', '.join(published) or 'nothing'}")
        else:
            # Without publish_data the persistent workdir is the result.
            a.log.append("note", f"kept in the workdir: {', '.join(published) or 'nothing'}")
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
