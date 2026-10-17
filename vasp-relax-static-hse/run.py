#!/usr/bin/env python3
"""Relax with HSE from a PBE-relaxed start, by calling ``vasp.relax-static-coarse`` once per functional.

This runner runs no VASP itself. It composes two jobs of the coarse-mesh driver
through httk's workflow calls (``[workflow.calls]``, ``Attempt.call`` and
``Attempt.gather``), so the optimizer, the offset correction, the time guard and
the remedy ladder are that package's, unchanged:

* ``mode = "no_pbe"`` (default): ``start`` calls the driver on the input
  structure with its own PBE templates (child ``pbe``); ``after_pbe`` reads the
  PBE result and calls the driver again with HSE (child ``hse``).
* ``mode = "start_from_pbe"``: the input structure is already PBE-relaxed;
  ``start`` calls the HSE stage on it directly. ``pbe_provenance`` (a record id,
  job key or other pointer to where the structure came from) is required, and
  ``files/WAVECAR.pbe``, ``files/WAVECAR.pbe.final`` and ``files/hessian.pbe.json``
  in the payload are used as seeds when present.

What the HSE stage gets from the PBE stage:

* its starting structure: the PBE stage's ``relax/CONTCAR`` (the geometry of
  its final static);
* its k-point meshes: the PBE stage's coarse and final meshes, unless the job
  gives ``coarse_mesh`` / ``final_mesh`` itself (or in ``hse_parameters``), so
  the WAVECARs fit;
* with ``seed_wavecar`` (default), the PBE stage's last coarse-mesh and
  final-mesh WAVECARs as the HSE stage's first restarts, each only when the
  HSE mesh of that kind is the PBE one;
* with ``hessian_from_pbe`` (default off), the PBE stage's last BFGS Hessian
  as the HSE stage's ``hessian_seed``, projected by the driver onto the
  symmetry of the HSE start (``hessian_seed_symmetrize``). Otherwise it is used
  as it is: the coordinates of the two stages differ only by the PBE strain of
  the start (about 1 %), which changes the matrix by about as much.

The INCARs: the HSE stage runs ``files/INCAR.hse.coarse`` and
``files/INCAR.hse.final`` of the payload, else this package's templates of those
names; the PBE stage runs ``files/INCAR.pbe.coarse`` / ``files/INCAR.pbe.final``
when given, else the driver's own templates. Every other job parameter is passed
to both stages, then ``pbe_parameters`` to the PBE stage and ``hse_parameters``
to the HSE stage on top.

``stages.json`` in the workdir records per stage the child job key and id, the
called workflow's installed id, the input and seed sha256, the meshes and the
outcome, and the PBE-to-HSE geometry shift. ``after_hse`` copies the HSE stage's
results into the workdir in the driver's layout (``relax/`` and the static at
the top) and ``collect.py`` reads them there. A stage that fails fails this job
with ``vasp.pbe_stage_failed`` or ``vasp.hse_stage_failed`` and the child's
failure in the details; a PBE stage that ends outside its limits is still used
as the HSE start, and its ``within_limits`` is recorded.
"""

import hashlib
import json
import math
import shutil
from pathlib import Path

from httk.workflow import Attempt, Runner

TEMPLATES = Path(__file__).resolve().parent / "templates"
CALLEE = "coarse"
STAGES = "stages.json"
# The driver's files of the final static that the workdir top keeps, as vasp-relax-static-coarse publishes them.
COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json POSCAR"
# This runner's parameters; every other job parameter goes to both stages.
OWN = {
    "mode": "no_pbe",
    "pbe_parameters": {},
    "hse_parameters": {},
    "pbe_provenance": None,
    "seed_wavecar": True,
    "hessian_from_pbe": False,
    "poscar": "files/POSCAR",
    "potcar": "files/POTCAR",
}

run = Runner("vasp.relax-static-hse")


def parameter(a: Attempt, name: str):
    value = a.parameter(name, OWN[name])
    return OWN[name] if value is None else value


def flag(value: object) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line."""

    return value is True or str(value).strip().lower() in ("true", "1", "yes")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def shared_parameters(a: Attempt) -> dict:
    return {name: value for name, value in dict(a.parameters).items() if name not in OWN}


def read_stages(a: Attempt) -> dict:
    path = a.workdir / STAGES
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {"stages": {}}


def write_stages(a: Attempt, stages: dict) -> None:
    path = a.workdir / STAGES
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(stages, indent=1), encoding="utf-8")
    temporary.replace(path)


def lattice(poscar: Path) -> dict:
    """Lattice lengths (A) and volume (A^3) of a VASP-5 POSCAR with a positive scale."""

    lines = poscar.read_text(encoding="utf-8").splitlines()
    scale = float(lines[1].split()[0])
    cell = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    volume = abs(
        cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
        - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
        + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0])
    )
    return {"lengths": [math.sqrt(sum(x * x for x in row)) for row in cell], "volume": volume}


def child_record(child) -> dict:
    """What stages.json keeps of a finished child: identity, called workflow and outcome."""

    job = json.loads((child.payload / "job.json").read_text(encoding="utf-8")) if child.payload else {}
    report_path = None if child.workdir is None else child.workdir / "relax-driver" / "report.json"
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path and report_path.is_file() else None
    record = {
        "label": child.label,
        "job_key": child.job_key,
        "job_id": child.job_id,
        "kind": child.kind,
        "workflow": job.get("workflow", {}),
        "failure": None if child.failure is None else {"code": child.failure.code, "message": child.failure.message},
    }
    if report is not None:
        finals = [r for r in report["evaluations"] if r["evaluation"] == report["final_evaluation"]]
        record.update(
            within_limits=report["within_limits"],
            coarse_evaluations=report["coarse_evaluations"],
            corrections=report["corrections"],
            meshes=report["meshes"],
            final_energy=finals[0]["energy"] if finals else None,
            final_max_force=finals[0]["max_force"] if finals else None,
            final_max_stress=finals[0]["max_stress"] if finals else None,
            evaluations=len(report["evaluations"]),
            seconds=sum(r["seconds"] for r in report["evaluations"]),
        )
    return record


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def start(a: Attempt) -> None:
    """Call the PBE stage, or with a PBE-relaxed input the HSE stage directly."""

    mode = str(parameter(a, "mode"))
    if mode not in ("no_pbe", "start_from_pbe"):
        a.fail("vasp.input_invalid", f"mode must be no_pbe or start_from_pbe, not {mode!r}")
        return
    poscar = a.payload / str(parameter(a, "poscar"))
    if not poscar.is_file():
        a.fail("vasp.input_missing", f"the starting structure {poscar.name} is not in this payload")
        return
    potcar = a.payload / str(parameter(a, "potcar"))
    files = {"POSCAR": poscar} | ({"POTCAR": potcar} if potcar.is_file() else {})
    stages = {"workflow": "vasp.relax-static-hse", "mode": mode, "input_sha256": sha256(poscar), "stages": {}}

    if mode == "start_from_pbe":
        provenance = a.parameter("pbe_provenance", None)
        if not provenance:
            a.fail(
                "vasp.pbe_provenance_missing",
                "start_from_pbe needs pbe_provenance: where the PBE-relaxed structure came from",
            )
            return
        given = {
            "coarse": a.payload / "files" / "WAVECAR.pbe",
            "final": a.payload / "files" / "WAVECAR.pbe.final",
            "hessian": a.payload / "files" / "hessian.pbe.json",
        }
        stages["stages"]["pbe"] = {"source": "given", "provenance": provenance, "structure_sha256": sha256(poscar)}
        call_hse(a, stages, poscar, files.get("POTCAR"), {}, {k: p for k, p in given.items() if p.is_file()})
        return

    for kind in ("coarse", "final"):
        incar = a.payload / "files" / f"INCAR.pbe.{kind}"
        if incar.is_file():
            files[f"INCAR.{kind}"] = incar
    parameters = shared_parameters(a) | dict(parameter(a, "pbe_parameters") or {})
    a.call(CALLEE, label="pbe", files=files, parameters=parameters)
    stages["stages"]["pbe"] = {"source": "child", "parameters": parameters}
    write_stages(a, stages)
    a.gather("after_pbe", when="all_succeeded", on_impossible="triage")


@run.step
def after_pbe(a: Attempt) -> None:
    """Read the PBE stage's result and seeds, and call the HSE stage from it."""

    child = a.children["pbe"]
    stages = read_stages(a)
    record = child_record(child)
    relaxed = child.workdir / "relax" / "CONTCAR"
    if not relaxed.is_file():
        a.fail("vasp.pbe_stage_failed", f"the PBE stage {child.job_key} left no relax/CONTCAR")
        return
    record["structure_sha256"] = sha256(relaxed)
    stages["stages"]["pbe"].update(record)
    driver = child.workdir / "relax-driver"
    seeds = {}
    snapshot = driver / f"driver-{record['evaluations']:03d}.json"
    if snapshot.is_file():
        state = json.loads(snapshot.read_text(encoding="utf-8"))
        # The WAVECAR each mesh kind would restart from next, as the PBE stage's last snapshot names it.
        for kind, cache in (state.get("wavecars") or {}).items():
            if (driver / kind / cache["file"]).is_file():
                seeds[kind] = driver / kind / cache["file"]
        hessian = a.workdir / "hessian.pbe.json"
        matrix = state["optimizer"]["B"]
        hessian.write_text(
            json.dumps({"hessian": matrix, "source": f"{child.job_key}:{snapshot.name}"}), encoding="utf-8"
        )
        seeds["hessian"] = hessian
    potcar = a.payload / str(parameter(a, "potcar"))
    call_hse(a, stages, relaxed, potcar if potcar.is_file() else None, record.get("meshes") or {}, seeds)


def call_hse(a: Attempt, stages: dict, structure: Path, potcar: Path | None, meshes: dict, seeds: dict) -> None:
    """Call the HSE stage on *structure* with the PBE meshes and the seeds that fit, and wait for it."""

    files: dict[str, Path] = {"POSCAR": structure} | ({"POTCAR": potcar} if potcar is not None else {})
    for kind in ("coarse", "final"):
        override = a.payload / "files" / f"INCAR.hse.{kind}"
        files[f"INCAR.{kind}"] = override if override.is_file() else TEMPLATES / f"INCAR.hse.{kind}"
    parameters = shared_parameters(a) | dict(parameter(a, "hse_parameters") or {})
    used: dict[str, dict] = {}
    for kind in ("coarse", "final"):
        pbe_mesh = meshes.get(kind)
        if pbe_mesh and not str(parameters.get(f"{kind}_mesh") or "").split():
            parameters[f"{kind}_mesh"] = " ".join(pbe_mesh)
        hse_mesh = str(parameters.get(f"{kind}_mesh") or "").split()
        seed = seeds.get(kind)
        # A WAVECAR restarts only on its own mesh; without a known PBE mesh (a given seed) the job's mesh must be set.
        fits = seed is not None and bool(hse_mesh) and (pbe_mesh is None or hse_mesh == list(pbe_mesh))
        if fits and flag(parameter(a, "seed_wavecar")):
            name = "WAVECAR.pbe" if kind == "coarse" else "WAVECAR.pbe.final"
            files[name] = seed
            parameters[f"{kind}_wavecar"] = f"files/{name}"
            used[f"{kind}_wavecar"] = {"sha256": sha256(seed), "mesh": " ".join(hse_mesh)}
        else:
            # Never let the driver's default seed name pick up a stray file.
            parameters[f"{kind}_wavecar"] = f"files/WAVECAR.pbe.{kind}.unused"
    if flag(parameter(a, "hessian_from_pbe")) and "hessian" in seeds:
        files["hessian.pbe.json"] = seeds["hessian"]
        parameters["hessian_seed"] = "files/hessian.pbe.json"
        # Learned along the PBE path, which may have started at lower symmetry than the PBE result the HSE stage
        # starts from: the driver projects it onto that symmetry instead of refusing it.
        parameters["hessian_seed_symmetrize"] = True
        used["hessian_seed"] = {"sha256": sha256(seeds["hessian"])}
    a.call(CALLEE, label="hse", files=files, parameters=parameters)
    stages["stages"]["hse"] = {
        "source": "child",
        "parameters": parameters,
        "structure_sha256": sha256(structure),
        "incar_sha256": {kind: sha256(files[f"INCAR.{kind}"]) for kind in ("coarse", "final")},
        "seeds": used,
        "start_lattice": lattice(structure),
        "shift_from_start": None,
    }
    write_stages(a, stages)
    a.gather("after_hse", when="all_succeeded", on_impossible="triage")


@run.step
def after_hse(a: Attempt) -> None:
    """Copy the HSE stage's results into the workdir, record the geometry shift, and complete the job."""

    child = a.children["hse"]
    stages = read_stages(a)
    record = child_record(child)
    hse = stages["stages"]["hse"]
    hse.update(record)
    archive = a.workdir / "relax"
    if archive.is_dir():
        shutil.rmtree(archive)
    shutil.copytree(child.workdir / "relax", archive)
    for name in COLLECT.split():
        if (child.workdir / name).is_file():
            shutil.copyfile(child.workdir / name, a.workdir / name)
    shutil.copyfile(child.workdir / "relax-driver" / "report.json", a.workdir / "hse-report.json")
    end = lattice(archive / "CONTCAR")
    begin = hse["start_lattice"]
    hse["end_lattice"] = end
    hse["shift_from_start"] = {
        "lengths_percent": [100 * (e / b - 1) for e, b in zip(end["lengths"], begin["lengths"])],
        "volume_percent": 100 * (end["volume"] / begin["volume"] - 1),
    }
    write_stages(a, stages)
    base = a.parameter("data_prefix", "") or ""
    if a.parameter("publish_data", False) is True:
        for stage_name, directory in (("static", a.workdir), ("relax", archive)):
            prefix = f"{base}/{stage_name}" if base else stage_name
            for name in COLLECT.split() + (["record.json"] if stage_name == "relax" else []):
                if (directory / name).is_file():
                    a.put(directory / name, f"{prefix}/{name}")
        for name in (STAGES, "hse-report.json"):
            a.put(a.workdir / name, f"{base}/{name}" if base else name)
    a.state.merge(
        {
            "within_limits": record.get("within_limits"),
            "static_energy": record.get("final_energy"),
            "hse_coarse_evaluations": record.get("coarse_evaluations"),
            "pbe_coarse_evaluations": stages["stages"].get("pbe", {}).get("coarse_evaluations"),
            "volume_shift_percent": hse["shift_from_start"]["volume_percent"],
        }
    )
    a.log.append(
        "note",
        f"HSE stage {child.job_key} {'within' if record.get('within_limits') else 'OUTSIDE'} the limits after "
        f"{record.get('coarse_evaluations')} coarse evaluations; volume {hse['shift_from_start']['volume_percent']:+.2f} % "
        "from the PBE start",
    )
    a.succeed()


@run.step
def triage(a: Attempt) -> None:
    """Fail by the stage that failed, with the child's failure."""

    stages = read_stages(a)
    failed = a.children.failed
    child = failed[0] if failed else None
    label = child.label if child is not None else "unknown"
    if child is not None:
        stages["stages"].setdefault(label, {}).update(child_record(child))
        write_stages(a, stages)
    a.fail(
        f"vasp.{label}_stage_failed",
        f"the {label.upper()} stage did not succeed",
        details={} if child is None else stages["stages"][label],
    )


if __name__ == "__main__":
    raise SystemExit(run.main())
