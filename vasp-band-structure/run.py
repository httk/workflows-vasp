#!/usr/bin/env python3
"""The electronic band structure of a structure that is already chosen (normally relaxed).

The workflow is spelled out on the ``httk.codes.vasp`` primitives so that this
package directory is a starting point to copy and edit. The structure is
standardized to its HPKOT primitive cell (Hinuma et al., the convention of
SeeK-path) with ``httk.atomistic.band_path``; the band k-points refer to the
reciprocal basis of that cell, so it, and not the input cell, is the structure
both VASP runs compute. Then a self-consistent run on a uniform mesh writes the
charge density, and a non-self-consistent run (``ICHARG = 11``) reads it along
the line-mode path.

The steps, each with its stage fixed by the step itself, so that rerunning any
step after a crash repeats only that step's own work:

* ``prepare`` recognizes the symmetry of the input structure within
  ``symmetry_tolerance`` (Å; recognition symmetrizes the sites within it) and
  writes the HPKOT cell as the workdir POSCAR, in the input POSCAR's species
  order, the line-mode ``KPOINTS.bands`` and ``band_path.json``, the exact labels,
  coordinates and path (its ``transform`` maps the basis of
  ``httk.atomistic.primitive_cell`` of the input to the band cell's basis, not
  the input POSCAR's basis);
* ``stage_scf`` and ``run_scf`` prepare and run the self-consistent stage from
  ``templates/INCAR.scf``, on the automatic mesh of ``kpoint_density``, and
  archive it under ``scf/``, leaving its CHGCAR at the top of the workdir;
* ``stage_bands`` and ``run_bands`` prepare and run the band stage from
  ``templates/INCAR.bands`` on ``KPOINTS.bands``;
* ``publish`` publishes both stages.

Edit the templates of a copied package to change the settings, or replace one
in a single job with a payload file ``files/INCAR.scf`` or ``files/INCAR.bands``.
The template's own ``EDIFF`` and ``EDIFFG`` are kept; when absent they are
derived from ``accuracy_per_atom``. The ``incar_tags`` parameter applies to both
stages and wins over the templates and every derived value. A remedy applied in
the scf stage stays applied in the band stage, as in the other staged runners.

The guards, each failing the job by name. Before any VASP run, ``prepare``
checks ``incar_tags`` and each stage's payload INCAR or template:
``MAGMOM``, ``NELECT``, ``NUPDOWN`` and the absolute FFT grids ``NG[XYZ]``,
``NG[XYZ]F`` (or the POSCAR comment's ``[MAGMOM=...]``) are refused with
``vasp.cell_dependent_tag_unsupported``, since the HPKOT cell's sites and size
differ from the input's (a ``null`` in ``incar_tags``, which leaves the tag to
VASP's default, is no such override); any ``ICHARG`` in ``incar_tags``, a band INCAR's
``ICHARG`` other than 11, or a ``band_divisions`` that is not an integer of at
least 2 is ``vasp.invalid_band_inputs``; and a recognition failure is
``vasp.band_path_failed``. The stages check again what VASP is about to run:
a band stage whose ``ICHARG`` is not 11 is
``vasp.invalid_band_inputs``; one without the scf CHGCAR is
``vasp.charge_density_missing``; and one whose density-defining tags differ from
the scf run's is ``vasp.charge_density_incompatible``. Remedies that move the
structure (``scale_lattice``, ``contcar_to_poscar``) are never applied, nor, in
the band stage, an ``ICHARG`` change. The band stage always runs with
``ICHARG = 11`` and ``LCHARG = .FALSE.``, also when an scf remedy set
``ICHARG = 2``.

The scope is collinear LDA/GGA(+U): ``LHFCALC``, ``METAGGA``, ``LNONCOLLINEAR``
and ``LSORBIT`` fail with ``vasp.unsupported_calculation``, in ``prepare`` and
again on each stage's prepared INCAR. Hybrid functionals
need zero-weight explicit k-points in a self-consistent run, and meta-GGAs the
kinetic-energy density, which a non-self-consistent ``ICHARG = 11`` run does not
read.

The job inputs and parameters are documented in this repository's README.
Unlike the other runners, this one imports an installed *httk-atomistic* (with
spglib) besides *httk-workflow*; this file and its ``templates/`` are the whole
runner: it is the entry of this package directory, reference it by git URI
(e.g. ``git+https://github.com/httk/workflows-vasp#vasp-band-structure``) or with
``--workflow-dir`` (installed into the workspace first), or copy the package
directory and edit it. The runner file alone, without its ``templates/``, fails
with ``vasp.template_missing``.
"""

import json
import shlex
import shutil
from pathlib import Path

import httk.atomistic
import httk.core
from httk.atomistic.integrations.vasp.io import write_line_mode_kpoints
from httk.codes.vasp import (
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    read_incar,
    read_poscar_header,
    reapply_vasp_remedies,
    run_vasp,
    update_incar,
    validate_vasp_workdir,
)
from httk.workflow import Attempt, Runner

COLLECT = (
    "INCAR KPOINTS POSCAR OUTCAR OSZICAR EIGENVAL vasprun.xml band_path.json vasp-run-report.json "
    "POTCAR.provenance.json"
)
# Kept across a remedied rerun and between the stages: the band stage reads the scf CHGCAR.
KEEP_BETWEEN_RUNS = ("WAVECAR", "CHGCAR", "CHG")
# What scf/ always keeps of the self-consistent run, besides the collect files.
SCF_ARCHIVE = ("INCAR", "KPOINTS", "POSCAR", "OUTCAR", "OSZICAR", "vasprun.xml", "vasp-run-report.json")
# What publish_data always publishes per stage, besides the collect files: what collect.py reads.
REQUIRED = {"scf": ("OUTCAR",), "bands": ("POSCAR", "EIGENVAL", "band_path.json")}
# Per-site, cell-total or cell-sized (absolute FFT grid) tags whose meaning changes with the
# reduction to the HPKOT primitive cell.
CELL_DEPENDENT_TAGS = ("MAGMOM", "NELECT", "NUPDOWN", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF")
# Tags outside the collinear LDA/GGA(+U) scope; METAGGA in any value, the others when true.
SCOPE_TAGS = ("LHFCALC", "LNONCOLLINEAR", "LSORBIT", "METAGGA")
# The INCAR tags that define the charge density the band stage reads.
DENSITY_TAGS = (
    *("PREC", "ENCUT", "ENAUG", "NGX", "NGY", "NGZ", "NGXF", "NGYF", "NGZF", "ISPIN", "MAGMOM", "LMAXMIX"),
    *("GGA", "METAGGA", "LASPH", "LDAU", "LDAUTYPE", "LDAUL", "LDAUU", "LDAUJ", "NELECT", "LNONCOLLINEAR", "LSORBIT"),
    *("GGA_COMPAT", "LUSE_VDW", "AGGAC", "AGGAX", "ALDAC", "PARAM1", "PARAM2", "LMAXPAW", "XC"),
)
# Remedies never applied: the structure must stay the HPKOT cell, and the band stage must stay non-self-consistent.
EXCLUDED = {
    "scf": {"scale_lattice", "contcar_to_poscar"},
    "bands": {"scale_lattice", "contcar_to_poscar", "incar.ICHARG"},
}
TEMPLATES = Path(__file__).resolve().parent / "templates"

run = Runner("vasp.band-structure")


def is_true(value: str | None) -> bool:
    """Return whether an INCAR boolean reads true (``.TRUE.``, ``T``, ``true``, ...)."""

    return value is not None and value.strip().lstrip(".").upper().startswith("T")


def out_of_scope(tags: dict[str, str]) -> list[str]:
    """Return the tags of *tags* outside the collinear LDA/GGA(+U) scope."""

    return [tag for tag in SCOPE_TAGS if tag in tags and (tag == "METAGGA" or is_true(tags[tag]))]


@run.step
def prepare(a: Attempt) -> None:
    """Standardize the structure to its HPKOT primitive cell, write the band path, and go on to the scf stage."""

    validate_vasp_workdir(a.workdir)
    poscar = a.payload / a.parameter("poscar", "files/POSCAR")
    if not poscar.is_file():
        a.fail(
            "vasp.input_missing",
            f"the starting structure {poscar.name} is not in this payload",
            details={"expected": str(poscar)},
        )
        return
    header = read_poscar_header(poscar)
    # Checked here, before any VASP run, on every source of either stage's INCAR: the job's
    # incar_tags (applied to both stages) and each stage's payload override or template.
    explicit = {
        str(tag).strip().upper(): str(value)
        for tag, value in dict(a.parameter("incar_tags", {}) or {}).items()
        # A null tag is removed (VASP's default), which depends on no cell.
        if value is not None
    }
    sources = {"incar_tags": explicit}
    for name in ("scf", "bands"):
        override = a.payload / "files" / f"INCAR.{name}"
        source = override if override.is_file() else TEMPLATES / f"INCAR.{name}"
        # A missing template fails its stage with vasp.template_missing.
        if source.is_file():
            sources[f"{source.parent.name}/INCAR.{name}"] = read_incar(source)
    comment = {"MAGMOM": ""} if "[MAGMOM=" in header.comment else {}
    for source, tags in {**sources, f"the comment of {poscar.name}": comment}.items():
        for tag in CELL_DEPENDENT_TAGS:
            if tag in tags:
                a.fail(
                    "vasp.cell_dependent_tag_unsupported",
                    f"{tag} from {source} is not supported: the band structure runs on the HPKOT primitive cell, "
                    "whose sites and size differ from the input structure's",
                    details={"tag": tag, "source": source},
                )
                return
    # incar_tags apply to both stages, so any ICHARG there is wrong for one of them.
    for source, tags in sources.items():
        icharg = tags.get("ICHARG")
        if icharg is not None and (source == "incar_tags" or source.endswith("bands") and icharg.strip() != "11"):
            a.fail(
                "vasp.invalid_band_inputs",
                f"ICHARG = {icharg} from {source}: the scf stage starts from atomic charges and the band stage "
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
                "band structures only",
                details={"source": source, "tags": {tag: tags[tag] for tag in unsupported}},
            )
            return
    divisions = a.parameter("band_divisions", 40)
    if type(divisions) is not int or divisions < 2:
        a.fail(
            "vasp.invalid_band_inputs",
            f"band_divisions = {divisions!r}: the k-points per band-path segment must be an integer of at least 2",
            details={"band_divisions": divisions},
        )
        return
    tolerance = float(a.parameter("symmetry_tolerance", 1e-3))
    try:
        # By format, so any file name works, as in the other runners.
        structure = httk.core.load(str(poscar), format="vasp-poscar", precision=tolerance)
        result = httk.atomistic.band_path(structure, tolerance=tolerance)
    except ImportError as exception:
        a.fail(
            "vasp.band_path_failed",
            f"cannot recognize the symmetry of {poscar.name}: {exception} (pip install httk-atomistic[default])",
        )
        return
    except ValueError as exception:
        a.fail("vasp.band_path_failed", f"cannot find the band path of {poscar.name}: {exception}")
        return
    written = a.workdir / "POSCAR"
    httk.core.save(result.structure, str(written))
    # A supplied POTCAR, and LDAU* arrays, are per species in the input's order.
    species = read_poscar_header(written).species
    if species != header.species:
        a.fail(
            "vasp.input_mismatch",
            f"the HPKOT cell lists its species as {' '.join(species)}, not in the order of {poscar.name}",
            details={"input": list(header.species), "written": list(species)},
        )
        return
    potcar = a.payload / a.parameter("potcar", "files/POTCAR")
    if potcar.is_file():
        shutil.copyfile(potcar, a.workdir / "POTCAR")
    symbol = result.bravais_lattice_extended
    write_line_mode_kpoints(
        a.workdir / "KPOINTS.bands",
        result.points,
        result.path,
        divisions=divisions,
        comment=f"HPKOT {symbol} band path (httk)",
    )
    document = {
        "convention": "HPKOT",
        "bravais_lattice": result.bravais_lattice,
        "bravais_lattice_extended": symbol,
        "spacegroup": result.spacegroup.it_number,
        "symmetry_tolerance": tolerance,
        "points": {label: [str(value) for value in point.to_fractions()] for label, point in result.points.items()},
        "path": [list(segment) for segment in result.path],
        "transform": [[int(value) for value in row] for row in result.transform.to_fractions()],
    }
    (a.workdir / "band_path.json").write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    spacegroup = result.spacegroup.it_number
    a.log.append("note", f"standardized to the HPKOT {symbol} cell of space group {spacegroup}")
    a.advance("stage_scf", state={"bravais_lattice_extended": symbol, "spacegroup": spacegroup})


def stage_inputs(a: Attempt, name: str) -> bool:
    """Prepare one stage's inputs from its INCAR template; return whether they are supported."""

    incar = a.workdir / "INCAR"
    override = a.payload / "files" / f"INCAR.{name}"
    source = override if override.is_file() else TEMPLATES / f"INCAR.{name}"
    if not source.is_file():
        a.fail(
            "vasp.template_missing",
            f"the INCAR template {source} is missing: reference this workflow as its package directory "
            "(a git URI ending in #vasp-band-structure, or --workflow-dir), not as its runner file alone",
            details={"expected": str(source)},
        )
        return False
    shutil.copyfile(source, incar)
    template = read_incar(incar)
    # Only the first stage assembles a POTCAR from the library; the band stage finds it.
    library = None
    if not (a.workdir / "POTCAR").is_file():
        library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(a.parameter("kpoint_density", 20.0) or 20.0),
        centering=a.parameter("centering", VaspPreparationOptions.centering),
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
            "band structures only",
            details={"stage": name, "tags": {tag: tags[tag] for tag in unsupported}},
        )
        return False
    return True


@run.step
def stage_scf(a: Attempt) -> None:
    """Prepare the self-consistent stage on the uniform mesh, and go on to run it."""

    if stage_inputs(a, "scf"):
        a.log.append("note", "prepared the vasp.band-structure stage scf")
        a.advance("run_scf")


@run.step
def stage_bands(a: Attempt) -> None:
    """Prepare the non-self-consistent stage on the band path, and go on to run it."""

    if not stage_inputs(a, "bands"):
        return
    override = a.payload / "files" / "INCAR.bands"
    explicit = {
        str(tag).strip().upper(): str(value)
        for tag, value in dict(a.parameter("incar_tags", {}) or {}).items()
        # A null tag is removed (VASP's default), which depends on no cell.
        if value is not None
    }
    for source, tags in (
        ("incar_tags", explicit),
        ("files/INCAR.bands", read_incar(override) if override.is_file() else {}),
    ):
        if "ICHARG" in tags and tags["ICHARG"].strip() != "11":
            a.fail(
                "vasp.invalid_band_inputs",
                f"ICHARG = {tags['ICHARG']} from {source}: the band stage is non-self-consistent and needs ICHARG = 11",
                details={"source": source, "ICHARG": tags["ICHARG"]},
            )
            return
    incar = a.workdir / "INCAR"
    # A replayed scf remedy may have set ICHARG = 2; the scf CHGCAR must not be overwritten.
    update_incar({"ICHARG": 11, "LCHARG": ".FALSE."}, incar)
    # After the replay, which may have edited the automatic mesh.
    shutil.copyfile(a.workdir / "KPOINTS.bands", a.workdir / "KPOINTS")
    a.log.append("note", "prepared the vasp.band-structure stage bands")
    a.advance("run_bands")


def band_inputs_problem(a: Attempt) -> tuple[str, str, dict[str, object]] | None:
    """Return the failure that forbids launching the band stage, or ``None`` when it may run."""

    tags = read_incar(a.workdir / "INCAR")
    if tags.get("ICHARG", "").strip() != "11":
        return (
            "vasp.invalid_band_inputs",
            f"the band stage INCAR has ICHARG = {tags.get('ICHARG')}, not 11",
            {"ICHARG": tags.get("ICHARG")},
        )
    chgcar = a.workdir / "CHGCAR"
    if not chgcar.is_file() or not chgcar.stat().st_size:
        return ("vasp.charge_density_missing", "the scf stage left no CHGCAR for the band stage to read", {})
    scf = dict(a.state.get("scf_incar") or {})

    def normalized(value: object) -> str | None:
        return None if value is None else str(value).strip().lower()

    differences = {
        tag: {"scf": scf.get(tag), "bands": tags.get(tag)}
        for tag in DENSITY_TAGS
        if normalized(scf.get(tag)) != normalized(tags.get(tag))
    }
    if differences:
        listed = "; ".join(
            f"{tag} {pair['scf']!r} in scf, {pair['bands']!r} in bands" for tag, pair in differences.items()
        )
        return (
            "vasp.charge_density_incompatible",
            f"the band stage cannot read the scf charge density: {listed}",
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
    # A rerun must not read the previous run's outputs; the CHGCAR the band stage reads is kept.
    clean_vasp_outputs(a.workdir, keep=KEEP_BETWEEN_RUNS)
    if name == "bands":
        # Checked before every launch, retries included.
        problem = band_inputs_problem(a)
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
        if name == "bands":
            a.advance("publish", state={"bands_classification": report.classification})
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
        a.advance("stage_bands", state=state)
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
    """Run the self-consistent stage, archive it under ``scf/``, and go on to the band stage."""

    run_stage(a, "scf")


@run.step
def run_bands(a: Attempt) -> None:
    """Run the non-self-consistent band stage, and go on to publish."""

    run_stage(a, "bands")


@run.step
def publish(a: Attempt) -> None:
    """Publish both stages, and complete the job."""

    base = a.parameter("data_prefix", "") or ""
    # Results stay in the persistent workdir; publish_data opts a job into a curated data/ copy too.
    to_data = a.parameter("publish_data", False) is True
    collect = (a.parameter("collect", COLLECT) or "").split()
    for stage_name, directory in (("scf", a.workdir / "scf"), ("bands", a.workdir)):
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
