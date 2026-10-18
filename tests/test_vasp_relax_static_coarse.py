"""The ``vasp-relax-static-coarse`` driver, end to end through the real manager, against a model VASP.

The fake VASP here is not the scripted one of ``test_vasp_runners.py``: the
driver needs forces and stress that belong to an energy surface, so this one
evaluates a model with analytic forces and stress at whatever geometry the
POSCAR holds and writes them to ``vasprun.xml`` as VASP does (stress in kbar,
positive = compressive). The model is a harmonic cell term in the metric tensor
(``E = kc/4 |C C^T - G*|^2``) plus a harmonic term per atom around its target
site, optionally a double well along x for the first atom (the soft mode). The
final k-point mesh moves the targets slightly, which is the coarse-to-final
offset the driver must measure and correct. Forces and stress carry noise
proportional to ``EDIFF`` (seeded from the inputs, so every run is reproducible),
and a tight EDIFF can shift energy, forces and stress (``tight_shift``). The
energy can carry a term linear in the volume that the stress does not
(``pulay``, kbar). SCF
counts are lower when a WAVECAR is read; a WAVECAR is read only if it was
written on the same mesh, and VASP's stdout says whether it was, as VASP 6.6
does. A run can be made to fail (ZPOTRF), to reach NELM, to ignore its WAVECAR,
to say nothing about it, to sleep, or to exit with a status and nothing else.

The tests marked as fault injections call the runner's steps directly with a
stand-in for :class:`httk.workflow.Attempt`, so an attempt can be stopped at an
exact point and repeated.
"""

import importlib.util
import itertools
import json
import math
import time
from pathlib import Path
from typing import Any

import pytest
from httk.workflow import TaskManager, Workspace
from httk.workflow.introspection import resolve_job
from httk.workflow.scaffold import new_job

REPO_ROOT = Path(__file__).resolve().parent.parent
PACKAGE = REPO_ROOT / "vasp-relax-static-coarse"

_POSCAR = """model
1.0
4.1000000000 0.0000000000 0.0000000000
0.0500000000 4.0500000000 0.0000000000
0.0000000000 0.0000000000 4.1500000000
Si
2
Direct
0.0200000000 0.0100000000 0.0000000000
0.5000000000 0.5000000000 0.4800000000
"""
# A body-centred cubic start, 2.5 % too large: 48 rotations times 2 translations.
_CUBIC_POSCAR = """cubic
1.0
4.1000000000 0.0000000000 0.0000000000
0.0000000000 4.1000000000 0.0000000000
0.0000000000 0.0000000000 4.1000000000
Si
2
Direct
0.0000000000 0.0000000000 0.0000000000
0.5000000000 0.5000000000 0.5000000000
"""
_CUBIC_MODEL = {"gstar": [[16.0, 0.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 16.0]], "final_site_shift": [0.0, 0.0, 0.0]}
_POTCAR = "  TITEL  = PAW_PBE Si 05Jan2001\n   ZVAL   =    4.000    mass and valenz\n"

_MODEL = {
    "gstar": [[16.0, 0.0, 0.0], [0.0, 16.0, 0.0], [0.0, 0.0, 17.64]],
    "sstar": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    "kc": 0.05,
    "ka": 10.0,
    "soft": None,
    "final_mesh": ["6", "6", "6"],
    "final_metric_scale": 1.002,
    "final_site_shift": [0.0005, 0.0, 0.0],
    "noise_force": 1000.0,
    "noise_stress": 1.0e5,
    "tight_shift": None,
    "pulay": 0.0,
    "fail_on": [],
    "nelm_on": [],
    "ignore_wavecar": [],
    "quiet": False,
    "sleep_on": {},
    "exit_on": {},
}

_FAKE_VASP = r'''#!/usr/bin/env python3
"""A model VASP: analytic energy, forces and stress of the model at the POSCAR geometry."""

import hashlib
import json
import random
import sys
import time
from pathlib import Path

MODEL = json.loads(MODEL_JSON)
here = Path(__file__).resolve().parent
if (here / "fake-sleep").is_file():
    time.sleep(float((here / "fake-sleep").read_text()))
attempts = Path("fake-vasp-attempts")
count = int(attempts.read_text()) if attempts.is_file() else 0
attempts.write_text(str(count + 1))
if str(count) in MODEL["sleep_on"]:
    time.sleep(MODEL["sleep_on"][str(count)])

tags = {}
for line in Path("INCAR").read_text().splitlines():
    for part in line.split("#")[0].split("!")[0].split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            tags[key.strip().upper()] = value.strip()
poscar = Path("POSCAR").read_text()
kpoints = Path("KPOINTS").read_text().splitlines()
mesh = kpoints[3].split()
wavecar = Path("WAVECAR").read_text() if Path("WAVECAR").is_file() else None
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(dict(tags=tags, POSCAR=poscar, mesh=mesh, centering=kpoints[2], WAVECAR=wavecar)) + "\n")
if str(count) in MODEL["exit_on"]:
    # An ordinary failure: no outputs, no diagnosis, a nonzero status.
    raise SystemExit(MODEL["exit_on"][str(count)])

if count in MODEL["fail_on"]:
    if Path("WAVECAR").is_file():
        Path("WAVECAR").write_text("")
    print("LAPACK: Routine ZPOTRF failed!")
    Path("OUTCAR").write_text(" fake vasp 6.4.1\n   NELM   =     60\n   NSW    =     0\n")
    Path("OSZICAR").write_text("DAV:   1    -0.100000000000E+02\n")
    raise SystemExit(1)
# VASP reads a WAVECAR only when it fits; otherwise it starts afresh and says so.
read = tags.get("ISTART") == "1" and wavecar == "wave " + " ".join(mesh) and count not in MODEL["ignore_wavecar"]
if not MODEL["quiet"]:
    print(" the WAVECAR file was read successfully" if read else " WAVECAR not read")

lines = poscar.splitlines()
scale = float(lines[1])
cell = [[scale * float(x) for x in lines[i].split()] for i in (2, 3, 4)]
n = sum(int(x) for x in lines[6].split())
frac = [[float(x) for x in lines[8 + i].split()[:3]] for i in range(n)]
final = mesh == MODEL["final_mesh"]
gstar = [[g * (MODEL["final_metric_scale"] if final else 1.0) for g in row] for row in MODEL["gstar"]]
sstar = [list(row) for row in MODEL["sstar"]]
if final:
    sstar[0] = [s + d for s, d in zip(sstar[0], MODEL["final_site_shift"])]


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def tr(a):
    return [list(row) for row in zip(*a)]


volume = (
    cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
    - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
    + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0])
)
delta = [[g - s for g, s in zip(grow, srow)] for grow, srow in zip(mm(cell, tr(cell)), gstar)]
# A negative baseline, since a positive final energy is diagnosed as a failure.
energy = -10.0 + 0.25 * MODEL["kc"] * sum(d * d for row in delta for d in row)
# dE/deps at fixed fractional coordinates.
dede = [[MODEL["kc"] * x for x in row] for row in mm(mm(tr(cell), delta), cell)]
# The two atoms interact through their relative displacement d from the target, so the forces sum to zero as
# VASP's do: harmonic, or a double well along x for the soft mode.
assert n == 2
ds = [((f0 - s0) - (f1 - s1) + 0.5) % 1.0 - 0.5 for f0, s0, f1, s1 in zip(frac[0], sstar[0], frac[1], sstar[1])]
d = mm([ds], cell)[0]
grad = [MODEL["ka"] * x for x in d]
weights = [0.5 * MODEL["ka"] * x * x for x in d]
if MODEL["soft"] is not None:
    depth, x0 = MODEL["soft"]
    weights[0] = depth * (d[0] ** 2 - x0 ** 2) ** 2
    grad[0] = 4 * depth * d[0] * (d[0] ** 2 - x0 ** 2)
energy += sum(weights)
forces = [[-g for g in grad], list(grad)]
for i in range(3):
    for j in range(3):
        dede[i][j] += 0.5 * (d[i] * grad[j] + d[j] * grad[i])
stress = [[-x / volume * 1602.1766208 for x in row] for row in dede]
# An energy the stress does not belong to: -pulay (kbar) times the volume, as J-0210's energy and stress disagreed
# by an isotropic 2.15 kbar.
energy -= MODEL["pulay"] * volume / 1602.1766208

ediff = float(tags.get("EDIFF", "1e-6"))
if MODEL["tight_shift"] is not None and ediff < 1e-7:
    # What a tighter EDIFF changes at one geometry: the energy, a force pair along x, the diagonal stress.
    shift_energy, shift_force, shift_stress = MODEL["tight_shift"]
    energy += shift_energy
    forces[0][0] += shift_force
    forces[1][0] -= shift_force
    for i in range(3):
        stress[i][i] += shift_stress
rng = random.Random(hashlib.sha256((poscar + tags.get("EDIFF", "") + " ".join(mesh)).encode()).hexdigest())
forces = [[f + MODEL["noise_force"] * ediff * rng.gauss(0, 1) for f in row] for row in forces]
noise = [[MODEL["noise_stress"] * ediff * rng.gauss(0, 1) for _ in range(3)] for _ in range(3)]
stress = [[stress[i][j] + 0.5 * (noise[i][j] + noise[j][i]) for j in range(3)] for i in range(3)]
energy += ediff * rng.gauss(0, 1) + (0.01 if final else 0.0)

scf = 60 if count in MODEL["nelm_on"] else (6 if read else 16)
Path("OUTCAR").write_text(
    " fake vasp 6.4.1\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\n"
    "   NSW    =     0    number of steps for IOM\n"
    "   maximum number of plane-waves:    1234\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\n"
    "   free  energy   TOTEN  =       %.8f eV\n"
    "   energy  without entropy=      %.8f  energy(sigma->0) =      %.8f\n"
    " General timing and accounting informations for this job:\n" % (energy, energy, energy)
)
oszicar = ["       N       E                     dE             d eps       ncg     rms"]
for step in range(1, scf + 1):
    oszicar.append("DAV:  %d    %.12E   -0.10000E-06   -0.10000E-06   128   0.100E-03" % (step, energy))
oszicar.append("   1 F= %.8E E0= %.8E  d E =-.105000E+02" % (energy, energy))
Path("OSZICAR").write_text("\n".join(oszicar) + "\n")
Path("CONTCAR").write_text(poscar)


def varray(name, rows):
    return '<varray name="%s">%s</varray>' % (name, "".join("<v>%s</v>" % " ".join("%.10f" % x for x in row) for row in rows))


Path("vasprun.xml").write_text(
    "<modeling><calculation>"
    + "<scstep/>" * scf
    + '<structure><crystal>%s</crystal>%s</structure>' % (varray("basis", cell), varray("positions", [[x % 1.0 for x in row] for row in frac]))
    + varray("forces", forces)
    + varray("stress", stress)
    + '<energy><i name="e_fr_energy">%.10f</i></energy>' % energy
    + "</calculation></modeling>\n"
)
if count in MODEL["nelm_on"]:
    # An unconverged run leaves a spoiled WAVECAR behind.
    Path("WAVECAR").write_text("")
elif tags.get("LWAVE", "").strip(".").upper().startswith("T"):
    Path("WAVECAR").write_text("wave " + " ".join(mesh))
'''

_PARAMETERS = {"coarse_mesh": "3 3 3", "final_mesh": "6 6 6", "centering": "Gamma"}


def _fake(root: Path, model: dict[str, Any] | None) -> Path:
    executable = root / "fake-vasp"
    executable.write_text(
        _FAKE_VASP.replace("MODEL_JSON", repr(json.dumps({**_MODEL, **(model or {})}))), encoding="utf-8"
    )
    executable.chmod(0o755)
    return executable


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    model: dict[str, Any] | None = None,
    parameters: dict[str, object] | None = None,
    sleep: float | None = None,
    end_time_after: float | None = None,
    poscar: str = _POSCAR,
    files: dict[str, str] | None = None,
) -> tuple[Workspace, str, Path]:
    """Submit one job of the package against the model VASP, run it, and return the workspace, job and fake's dir.

    ``end_time_after`` gives the manager an allocation that ends that many seconds from now (no margin), and
    ``sleep`` makes every model VASP run take that long. ``files`` adds payload files by name and content.
    """

    root.mkdir(parents=True)
    executable = _fake(root, model)
    if sleep is not None:
        (root / "fake-sleep").write_text(str(sleep), encoding="utf-8")
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(root / "workspace")
    staged = root / "inputs"
    staged.mkdir()
    staged_files = {}
    for name, content in (("POSCAR", poscar), ("POTCAR", _POTCAR), *(files or {}).items()):
        (staged / name).write_text(content, encoding="utf-8")
        staged_files[name] = staged / name
    job = new_job(
        workspace,
        PACKAGE,
        files=staged_files,
        parameters={"publish_data": True, **_PARAMETERS, **(parameters or {})},
        name="packaged vasp-relax-static-coarse",
        install=True,
    )
    _run(workspace, end_time_after)
    return workspace, job.job_id, root


def _run(workspace: Workspace, end_time_after: float | None = None) -> None:
    end_time = None if end_time_after is None else time.time() + end_time_after
    with TaskManager(workspace, heartbeat_interval=0.01, end_time=end_time, deadline_margin=0.0) as manager:
        manager.run_until_idle(timeout=300.0)


def _job(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _state(payload: Path) -> dict[str, Any]:
    return json.loads((payload / ".httk-job" / "state.json").read_text(encoding="utf-8"))


def _failure(payload: Path) -> dict[str, Any]:
    return json.loads((payload / "state.json").read_text(encoding="utf-8"))["failure"]


def _runs(workdir: Path) -> list[dict[str, Any]]:
    path = workdir / "fake-vasp-runs.jsonl"
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def _report(workdir: Path) -> dict[str, Any]:
    return json.loads((workdir / "relax-driver" / "report.json").read_text(encoding="utf-8"))


def _records(workdir: Path) -> list[dict[str, Any]]:
    """The records of the evaluations so far, before the report exists."""

    paths = sorted((workdir / "relax-driver" / "evaluations").glob("*/record.json"))
    return [json.loads(path.read_text(encoding="utf-8")) for path in paths]


# -- a stand-in Attempt, for calling the runner's steps directly ------------------------------------------------------


def _module():
    """A fresh import of the package's run.py."""

    spec = importlib.util.spec_from_file_location("vasp_relax_static_coarse_run", PACKAGE / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _JobState(dict):
    def merge(self, values: dict[str, Any]) -> None:
        self.update(values)


class _Log:
    def __init__(self) -> None:
        self.entries: list[tuple[str, str]] = []

    def append(self, kind: str, message: str) -> None:
        self.entries.append((kind, message))


class _Context:
    def __init__(self, resources: dict[str, Any]) -> None:
        self.resources = resources


class _Attempt:
    """What the runner's steps use of an Attempt; the outcome is kept instead of published."""

    def __init__(self, payload: Path, parameters: dict[str, Any], resources: dict[str, Any], deadline: float | None):
        self.payload = payload
        self.workdir = payload / "run"
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.parameters = parameters
        self.state = _JobState()
        self.log = _Log()
        self.context = _Context(resources)
        self.deadline = deadline
        self.outcome: tuple | None = None

    def parameter(self, name: str, default: Any = None) -> Any:
        return self.parameters.get(name, default)

    def setting(self, name: str, default: Any = None) -> Any:
        return default

    def advance(self, step: str, *, state=None, resources=None, priority=None) -> None:
        assert self.outcome is None
        if state:
            self.state.merge(state)
        self.outcome = ("advance", step, resources)

    def fail(self, code: str, message: str, *, details=None, retryable=False, priority=None) -> None:
        assert self.outcome is None
        self.outcome = ("fail", code, message, details)

    def retry(self, reason: str) -> None:
        self.outcome = ("retry", reason)

    def put(self, source, destination) -> None:
        pass

    def succeed(self) -> None:
        self.outcome = ("succeed",)


def _attempt(
    root: Path,
    *,
    model: dict[str, Any] | None = None,
    parameters: dict[str, Any] | None = None,
    poscar: str = _POSCAR,
    resources: dict[str, Any] | None = None,
    deadline: float | None = None,
    base: dict[str, Any] | None = None,
) -> _Attempt:
    root.mkdir(parents=True)
    executable = _fake(root, model)
    payload = root / "payload"
    (payload / "files").mkdir(parents=True)
    (payload / ".httk-job").mkdir()
    (payload / "files" / "POSCAR").write_text(poscar, encoding="utf-8")
    (payload / "files" / "POTCAR").write_text(_POTCAR, encoding="utf-8")
    values = {"vasp_command": str(executable), **(_PARAMETERS if base is None else base), **(parameters or {})}
    return _Attempt(payload, values, resources or {}, deadline)


def _step(module, a: _Attempt, step: str) -> tuple:
    a.outcome = None
    getattr(module, step)(a)
    assert a.outcome is not None
    return a.outcome


def _drive(module, a: _Attempt, step: str = "prepare", limit: int = 200) -> tuple:
    """Run steps as a manager would until the job ends or stops; return the last outcome."""

    for _ in range(limit):
        outcome = _step(module, a, step)
        if outcome[0] == "advance":
            step = outcome[1]
        elif outcome[0] != "retry":
            return outcome
    raise AssertionError("the job did not end")


# -- the relaxation ---------------------------------------------------------------------------------------------------


def test_the_driver_relaxes_on_the_coarse_mesh_and_corrects_the_offset_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "offset", monkeypatch)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    report = _report(workdir)
    records = report["evaluations"]
    runs = _runs(workdir)
    assert len(runs) == len(records)

    # Every evaluation is a fresh single point; the coarse ones on the coarse mesh at PRECFOCK Fast, the final ones
    # on the final mesh at PRECFOCK Normal and the final EDIFF.
    for run, record in zip(runs, records):
        assert (run["tags"]["NSW"], run["tags"]["IBRION"], run["tags"]["LWAVE"]) == ("0", "-1", ".TRUE.")
        assert run["mesh"] == (["3", "3", "3"] if record["kind"] == "coarse" else ["6", "6", "6"])
        assert run["tags"]["PRECFOCK"] == ("Fast" if record["kind"] == "coarse" else "Normal")
        assert record["mesh"] == " ".join(run["mesh"])
        assert float(run["tags"]["EDIFF"]) == record["EDIFF"]
    # Two-level EDIFF: loose until the accepted point is within 0.02 eV/A and 1 kbar, then tight, and the switch
    # point is evaluated again at the tight EDIFF, which becomes the accepted point.
    switch = report["switched_at"]
    assert switch is not None
    coarse = [r for r in records if r["kind"] == "coarse"]
    assert all(r["EDIFF"] == 1e-6 for r in coarse if r["evaluation"] <= switch)
    assert all(r["EDIFF"] == 1e-8 for r in coarse if r["evaluation"] > switch)
    tighten = records[switch]
    assert tighten["purpose"] == "tighten" and tighten["y"] == records[switch - 1]["y"]
    assert tighten["rebased"] is True and tighten["accepted"] is True
    assert tighten["ediff_shift"]["from_evaluation"] == switch
    assert records[switch - 1]["max_force"] <= 0.02 and records[switch - 1]["max_stress"] <= 1.0

    # The first final static missed the limits, its offset was measured against the tight coarse evaluation at
    # the same geometry, the coarse surface was corrected once, and the second final static is within the limits.
    finals = [r for r in records if r["kind"] == "final"]
    assert [r["purpose"] for r in finals] == ["final", "final"]
    assert finals[0]["within_limits"] is False and finals[1]["within_limits"] is True
    assert report["within_limits"] is True and report["corrections"] == 1
    assert len(report["offsets"]) == 2
    for offset, final in zip(report["offsets"], finals):
        partner = records[offset["coarse_evaluation"] - 1]
        assert partner["kind"] == "coarse" and partner["y"] == final["y"]
        assert partner["EDIFF"] == 1e-8 and final["EDIFF"] == 1e-8
        assert offset["final_evaluation"] == final["evaluation"]
    # The model's offset: the final mesh stretches the target metric by 0.2 % and moves one site by 0.002 A.
    assert report["offsets"][0]["max_stress"] > 0.2 and report["offsets"][0]["max_force"] > 0.003
    corrected = [r for r in coarse if r["offset_applied"]]
    assert corrected and all(r["evaluation"] > finals[0]["evaluation"] for r in corrected)
    assert corrected[-1]["corrected_max_force"] <= 0.003 and corrected[-1]["corrected_max_stress"] <= 0.2
    assert finals[1]["max_force"] <= 0.003 and finals[1]["max_stress"] <= 0.2
    # The corrected-surface model and the physical prediction stay within a tenth of the limits along this path.
    assert all("correction_model_discrepancy" in r for r in corrected)
    bound = report["correction_model_bound"]
    assert bound["within"] is True and 0 < bound["force"] <= 0.0003 and bound["stress"] <= 0.02

    # WAVECAR restarts stay on their mesh: every evaluation after the first of its kind asks for its kind's WAVECAR
    # from the evaluation before it, and VASP's output shows it was read.
    first = {}
    previous = {}
    for run, record in zip(runs, records):
        first.setdefault(record["kind"], record["evaluation"])
        restarted = record["evaluation"] != first[record["kind"]]
        assert record["restart_requested"] is restarted and record["restart_observed"] is restarted
        assert record["restart_source"] == (previous[record["kind"]] if restarted else None)
        assert run["tags"]["ISTART"] == ("1" if restarted else "0")
        assert run["tags"]["ICHARG"] == ("0" if restarted else "2")
        assert record["scf_steps"] == (6 if restarted else 16)
        previous[record["kind"]] = record["evaluation"]
    assert report["restarts"] == {"requested": len(records) - 2, "observed": len(records) - 2, "unknown": 0}
    # Only the WAVECAR each kind restarts from is kept.
    assert sorted(p.name for p in (workdir / "relax-driver" / "coarse").iterdir()) == [
        f"WAVECAR-{previous['coarse']:03d}"
    ]

    # The layout collect.py reads: relax/ holds the coarse evaluation at the final geometry, the final static is at
    # the top, and the published data carries both and the report.
    assert (
        json.loads((workdir / "relax" / "record.json").read_text())["evaluation"]
        == report["offsets"][-1]["coarse_evaluation"]
    )
    assert (workdir / "relax" / "CONTCAR").read_text() == (workdir / "POSCAR").read_text()
    assert (payload / "data" / "relax-driver" / "report.json").is_file()
    assert (payload / "data" / "static" / "OUTCAR").is_file() and (payload / "data" / "relax" / "CONTCAR").is_file()
    state = _state(payload)
    assert state["within_limits"] is True
    assert float(state["static_energy"]) == pytest.approx(finals[1]["energy"])


def test_a_soft_mode_relaxes_into_one_well_without_oscillating(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # A double well along x for the first atom, started on its negative-curvature shoulder (x = 0.05 A, wells at
    # +-0.15 A, barrier 5.6 meV), with an initial Hessian model 35 times too stiff for that mode.
    model = {"soft": [11.1, 0.15], "sstar": [[-0.05 / 4.1 + 0.02, 0.01, 0.0], [0.5, 0.5, 0.5]]}
    workspace, job_id, _ = _campaign(tmp_path / "soft", monkeypatch, model=model)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    report = _report(payload / "run")
    coarse = [r for r in report["evaluations"] if r["kind"] == "coarse"]
    accepted = [r for r in coarse if r["accepted"]]
    # The accepted energies only go down (within the rise tolerance) at each EDIFF level, and the first atom never
    # crosses the barrier.
    for level in (1e-6, 1e-8):
        energies = [r["energy"] for r in accepted if r["EDIFF"] == level and not r.get("rebased")]
        assert all(b <= a + 2e-4 for a, b in itertools.pairwise(energies))
    sstar = model["sstar"]
    displacement = [
        sum(((r["frac"][0][k] - sstar[0][k]) - (r["frac"][1][k] - sstar[1][k])) * r["cell"][k][0] for k in range(3))
        for r in accepted
    ]
    assert all(x > 0 for x in displacement)
    assert abs(displacement[-1] - 0.15) < 0.01
    assert report["within_limits"] is True
    assert len(coarse) < 60


def test_an_energy_rise_is_rejected_and_the_trust_radius_halved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # An initial Hessian model 20 times too soft for the atoms: the Newton step overshoots once the trust radius no
    # longer caps it, the energy rises by more than the 0.2 meV tolerance, and the optimizer must back off.
    workspace, job_id, _ = _campaign(tmp_path / "rise", monkeypatch, parameters={"hessian_atom": 0.5})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    coarse = [r for r in _report(payload / "run")["evaluations"] if r["kind"] == "coarse"]
    rejected = [i for i, r in enumerate(coarse) if r["accepted"] is False]
    assert rejected
    for i in rejected:
        last_accepted = next(r for r in reversed(coarse[:i]) if r["accepted"])
        # Only ever compared within one EDIFF level.
        assert coarse[i]["EDIFF"] == last_accepted["EDIFF"]
        assert coarse[i]["energy"] > last_accepted["energy"] + 2e-4
        assert coarse[i]["trust_factor"] == pytest.approx(0.5 * coarse[i - 1]["trust_factor"])
    assert _report(payload / "run")["within_limits"] is True


def test_more_coarse_evaluations_than_allowed_fail_by_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "cap", monkeypatch, parameters={"maximum_evaluations": 4})
    kind, payload = _job(workspace, job_id)
    assert kind == "failed"
    assert _failure(payload)["code"] == "vasp.relax_not_converged"
    assert int(_state(payload)["evaluations"]) == 3
    assert len(_runs(payload / "run")) == 4


def test_an_offset_probe_measures_the_offset_at_the_switch_point_too(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "probe", monkeypatch, parameters={"offset_probe": True})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    report = _report(payload / "run")
    records = report["evaluations"]
    switch = report["switched_at"]
    probe = records[switch + 1]
    assert (probe["kind"], probe["purpose"], probe["y"]) == ("final", "probe", records[switch]["y"])
    assert report["offsets"][0]["purpose"] == "probe"
    assert report["offsets"][0]["coarse_evaluation"] == switch + 1
    assert [o["purpose"] for o in report["offsets"]][1:] == ["final", "final"]


# -- R1: a tighter EDIFF replaces the baseline ------------------------------------------------------------------------


@pytest.mark.parametrize("sign", (1.0, -1.0))
def test_a_tighter_ediff_replaces_the_accepted_point_instead_of_competing_with_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, sign: float
) -> None:
    # At the tight EDIFF the model's energy moves by 1 meV (five times the rise tolerance), a force pair by
    # 0.004 eV/A and the diagonal stress by 0.3 kbar, either way.
    shift = [sign * 0.001, sign * 0.004, sign * 0.3]
    workspace, job_id, _ = _campaign(
        tmp_path / "shift", monkeypatch, model={"tight_shift": shift}, parameters={"offset_probe": True}
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    report = _report(payload / "run")
    records = report["evaluations"]
    switch = report["switched_at"]
    tighten = records[switch]
    # The tight repeat of the switch point is the new accepted point whatever its energy, and records the shift.
    assert tighten["purpose"] == "tighten" and tighten["accepted"] is True and tighten["rebased"] is True
    assert tighten["ediff_shift"]["energy"] == pytest.approx(sign * 0.001, abs=1e-5)
    # Up to the loose evaluation's noise: 0.001 eV/A and 0.1 kbar per component.
    assert tighten["ediff_shift"]["max_force"] == pytest.approx(0.004, abs=0.003)
    assert tighten["ediff_shift"]["max_stress"] == pytest.approx(0.3, abs=0.3)
    # The probe follows the tight baseline and pairs with it.
    probe_offset = report["offsets"][0]
    assert probe_offset["purpose"] == "probe" and probe_offset["coarse_evaluation"] == switch + 1
    assert probe_offset["coarse_EDIFF"] == 1e-8 == probe_offset["final_EDIFF"]
    # Every later coarse evaluation is tight, and the relaxation converges instead of rejecting every step.
    later = [r for r in records if r["kind"] == "coarse" and r["evaluation"] > switch + 1]
    assert later and all(r["EDIFF"] == 1e-8 for r in later)
    assert sum(1 for r in later if r["accepted"] is False) <= 2
    assert report["within_limits"] is True


# -- R2: the time guard -----------------------------------------------------------------------------------------------


def test_the_time_guard_defers_an_evaluation_and_a_later_manager_resumes_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The same relaxation run without interruption, for the reference trajectory.
    reference_ws, reference_id, _ = _campaign(tmp_path / "reference", monkeypatch)
    reference = _report(_job(reference_ws, reference_id)[1] / "run")["evaluations"]

    # Every model run takes 2 s and the manager's allocation ends 9 s from now: evaluations stop once the time left
    # is below what the next one needs, before the allocation ends, and the job waits with that mintime.
    workspace, job_id, root = _campaign(tmp_path / "guard", monkeypatch, sleep=2.0, end_time_after=9.0)
    kind, payload = _job(workspace, job_id)
    assert kind == "ready"
    stopped = int(_state(payload)["evaluations"])
    assert 1 <= stopped < len(reference)
    # The next evaluation waits with a mintime of 1.4 x its expected time plus the measured overhead.
    driver = json.loads((payload / "run" / "relax-driver" / f"driver-{stopped:03d}.json").read_text())
    assert min(r["seconds"] for r in _records(payload / "run")) >= 2.0
    overhead = driver["overhead"]
    key = "coarse:loose" if driver["next"]["ediff"] > 1e-7 else "coarse:tight"
    # The estimate itself (and the first tight one's factor) is tested with fixed timings below.
    expected = _module().expected_seconds(_attempt(tmp_path / "estimate", parameters={}), driver, key)
    full = 1.4 * expected + overhead["archive"] + overhead["staging"]
    mintime = json.loads((payload / "state.json").read_text(encoding="utf-8"))["resources"]["mintime"]
    assert mintime in (math.ceil(full), math.ceil(full - overhead["staging"]))
    # No evaluation was cut off: every VASP run so far is a finished evaluation.
    assert len(_runs(payload / "run")) == stopped

    # A manager without an allocation end resumes the relaxation where it stopped, and it ends exactly as the
    # uninterrupted one did.
    (root / "fake-sleep").unlink()
    _run(workspace)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    records = _report(payload / "run")["evaluations"]
    assert [r["y"] for r in records] == [r["y"] for r in reference]
    assert [r["kind"] for r in records] == [r["kind"] for r in reference]


def test_an_attempt_killed_during_vasp_is_resumed_from_its_last_committed_evaluation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference_ws, reference_id, _ = _campaign(tmp_path / "reference", monkeypatch)
    reference = _report(_job(reference_ws, reference_id)[1] / "run")["evaluations"]

    # The fifth run takes far longer than any estimate: the allocation ends while VASP runs, and the manager stops
    # the attempt. 15 s leave the first four (well under a second each) room on a loaded machine: with 6 s Altair saw
    # the guard defer the fifth before it started (O1 review, 2026-10-18).
    workspace, job_id, _ = _campaign(
        tmp_path / "killed", monkeypatch, model={"sleep_on": {"4": 60.0}}, end_time_after=15.0
    )
    kind, payload = _job(workspace, job_id)
    # The stopped run is not a VASP failure: no remedy, and the evaluation waits for a manager with at least
    # 1.4 times the time the stopped run took. Four evaluations are committed; the fifth run never finished.
    assert kind == "ready", _failure(payload) if kind == "failed" else kind
    assert int(_state(payload)["evaluations"]) == 4 and "remedies" not in _state(payload)
    assert (payload / "run" / "fake-vasp-attempts").read_text() == "5"
    assert len(_runs(payload / "run")) == 4
    stopped = json.loads((payload / "run" / "relax-driver" / "driver-004.json").read_text())["slowest"]["coarse:loose"]
    assert stopped > 1.0
    assert json.loads((payload / "state.json").read_text())["resources"]["mintime"] >= 1.4 * stopped
    cached = (payload / "run" / "relax-driver" / "coarse" / "WAVECAR-004").read_text()
    assert cached == "wave 3 3 3"
    # The manager's SIGTERM is what the stopped run's process report says, and its files are kept.
    kept = payload / "run" / "relax-driver" / "incomplete" / "005-01"
    incomplete = json.loads((kept / "record.json").read_text())
    assert incomplete["termination"].startswith("signal_15") and incomplete["handling"] == "deferred"
    assert (kept / "INCAR").is_file()

    # The next manager resumes at evaluation 5, and the relaxation ends as the uninterrupted one did.
    _run(workspace)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded", _failure(payload) if kind == "failed" else kind
    records = _report(payload / "run")["evaluations"]
    assert len(_runs(payload / "run")) == len(records)
    assert _runs(payload / "run")[4]["WAVECAR"] == cached
    assert [r["y"] for r in records] == [r["y"] for r in reference]


def test_an_evaluation_longer_than_the_attempt_maxtime_fails_before_it_starts(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "maxtime", parameters={"expected_seconds": 80.0}, resources={"maxtime": 100})
    outcome = _drive(module, a)
    assert outcome[:2] == ("fail", "vasp.insufficient_time")
    assert outcome[3]["needed_seconds"] == pytest.approx(112.0) and outcome[3]["maxtime_seconds"] == 100
    assert _runs(a.workdir) == [] and not (a.workdir / "INCAR").exists()


def test_an_evaluation_whose_successor_cannot_fit_maxtime_stops_after_committing(tmp_path: Path) -> None:
    # The first tight evaluation is estimated at 1e5 times the loose ones: no attempt of 100 s could run it.
    module = _module()
    a = _attempt(tmp_path / "successor", parameters={"tight_cost_factor": 1e5}, resources={"maxtime": 100})
    outcome = _drive(module, a)
    assert outcome[:2] == ("fail", "vasp.insufficient_time")
    assert outcome[3]["purpose"] == "tighten"
    records = _records(a.workdir)
    assert int(a.state["evaluations"]) == len(records) == len(_runs(a.workdir))
    driver = json.loads((a.workdir / "relax-driver" / f"driver-{len(records):03d}.json").read_text())
    assert driver["next"]["purpose"] == "tighten"


def test_the_first_evaluation_is_guarded_by_expected_seconds(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "first", parameters={"expected_seconds": 10.0}, deadline=time.time() + 5)
    _step(module, a, "prepare")
    outcome = _step(module, a, "evaluate")
    assert outcome[:2] == ("advance", "evaluate") and outcome[2]["mintime"] == "0:14"
    assert _runs(a.workdir) == []


def test_the_time_guard_is_checked_again_after_staging(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # Enough time before staging (3 s for 1.4 s), not after a staging that takes 2 s.
    module = _module()
    a = _attempt(tmp_path / "staging", parameters={"expected_seconds": 1.0}, deadline=time.time() + 3)
    original = module.stage_inputs

    def slow(*arguments):
        staged = original(*arguments)
        time.sleep(2.0)
        return staged

    monkeypatch.setattr(module, "stage_inputs", slow)
    _step(module, a, "prepare")
    outcome = _step(module, a, "evaluate")
    assert outcome[:2] == ("advance", "evaluate") and outcome[2]["mintime"] == "0:02"
    assert any("after staging" in message for _, message in a.log.entries)
    assert _runs(a.workdir) == []


# -- R3: remedies, meshes and restart files ---------------------------------------------------------------------------


def test_an_electronic_failure_is_remedied_with_its_charge_reset_and_the_cached_wavecar(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "remedy", monkeypatch, model={"nelm_on": [3]})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    runs = _runs(workdir)
    records = _report(workdir)["evaluations"]
    # Run 3 reached NELM and spoiled its WAVECAR; the retry is evaluation 4 at the same geometry, from the cached
    # WAVECAR of evaluation 3, with the remedy's ICHARG = 2 and mixing, which stay in every later evaluation.
    assert len(runs) == len(records) + 1
    assert runs[3]["POSCAR"] == runs[4]["POSCAR"]
    assert runs[3]["WAVECAR"] == runs[4]["WAVECAR"] == "wave 3 3 3"
    assert runs[3]["tags"]["ICHARG"] == "0"
    for run in runs[4:]:
        assert run["tags"]["ICHARG"] == "2" and float(run["tags"]["AMIX"]) == 0.1
    assert records[3]["ICHARG"] == "2" and ["incar.ICHARG", 2] in records[3]["remedies"]
    assert records[3]["restart_requested"] is True and records[3]["restart_observed"] is True
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"electronic_nonconvergence": 1}


def test_a_remedy_that_would_change_the_surface_stops_the_job(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ZPOTRF's first remedy scales the lattice, which the driver would override: no remedy is applied or retried.
    workspace, job_id, _ = _campaign(tmp_path / "zpotrf", monkeypatch, model={"fail_on": [3]})
    kind, payload = _job(workspace, job_id)
    assert kind == "failed"
    failure = _failure(payload)
    assert failure["code"] == "vasp.unsupported_remedy" and "scale_lattice" in failure["message"]
    assert len(_runs(payload / "run")) == 4
    assert not (payload / ".httk-job" / "vasp-remedies.json").exists()
    # The failed run truncated its WAVECAR; the cached one of evaluation 3 is intact.
    assert (payload / "run" / "relax-driver" / "coarse" / "WAVECAR-003").read_text() == "wave 3 3 3"


def test_a_remedy_history_that_changes_the_mesh_stops_before_vasp(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "mesh")
    history = {
        "format": "httk-vasp-remedy-history",
        "format_version": 2,
        "attempts": {"kpoints_class": 1},
        "events": [{"problem": "kpoints_class", "changes": [{"operation": "bump_kpoints", "value": 1}]}],
    }
    (a.payload / ".httk-job" / "vasp-remedies.json").write_text(json.dumps(history), encoding="utf-8")
    _step(module, a, "prepare")
    outcome = _step(module, a, "evaluate")
    # Refused from the history itself, before any input is written.
    assert outcome[:2] == ("fail", "vasp.unsupported_remedy")
    assert "bump_kpoints" in outcome[2]
    assert _runs(a.workdir) == [] and not (a.workdir / "KPOINTS").exists()


def test_a_wavecar_is_staged_only_for_the_inputs_it_was_written_for(tmp_path: Path) -> None:
    # No centering given: the package's default is Gamma, as tetrahedron integration needs.
    module = _module()
    a = _attempt(tmp_path / "signature", base={"coarse_mesh": "3 3 3", "final_mesh": "6 6 6"})
    _step(module, a, "prepare")
    assert _step(module, a, "evaluate")[:2] == ("advance", "evaluate")
    assert _step(module, a, "evaluate")[:2] == ("advance", "evaluate")
    # The band count changes: the cached WAVECAR is not staged, and the evaluation starts afresh.
    a.parameters["coarse_incar_tags"] = {"PRECFOCK": "Fast", "NBANDS": 24}
    assert _step(module, a, "evaluate")[:2] == ("advance", "evaluate")
    runs = _runs(a.workdir)
    records = _records(a.workdir)
    assert all(run["centering"].startswith("G") for run in runs)
    assert [run["tags"]["ISTART"] for run in runs] == ["0", "1", "0"]
    assert [r["restart_requested"] for r in records] == [False, True, False]
    assert runs[2]["WAVECAR"] is None
    driver = json.loads((a.workdir / "relax-driver" / "driver-003.json").read_text())
    assert driver["wavecars"]["coarse"]["evaluation"] == 3
    assert driver["wavecars"]["coarse"]["signature"]["NBANDS"] == "24"


# -- R4: restart evidence ---------------------------------------------------------------------------------------------


def test_a_restart_counts_as_observed_only_when_vasp_says_it_read_the_wavecar(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A payload seed written on another mesh, and a third run that ignores a fitting WAVECAR.
    workspace, job_id, _ = _campaign(
        tmp_path / "evidence", monkeypatch, model={"ignore_wavecar": [2]}, files={"WAVECAR": "wave 9 9 9"}
    )
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    records = _report(payload / "run")["evaluations"]
    assert (records[0]["restart_requested"], records[0]["restart_source"]) == (True, "seed")
    assert records[0]["restart_observed"] is False and records[0]["scf_steps"] == 16
    assert (records[1]["restart_requested"], records[1]["restart_observed"]) == (True, True)
    assert (records[2]["restart_requested"], records[2]["restart_observed"]) == (True, False)
    assert records[2]["scf_steps"] == 16

    # A VASP that says nothing either way leaves every requested restart unknown.
    workspace, job_id, _ = _campaign(tmp_path / "quiet", monkeypatch, model={"quiet": True})
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    report = _report(payload / "run")
    assert report["restarts"]["unknown"] == report["restarts"]["requested"] > 0
    assert report["restarts"]["observed"] == 0


# -- R5: the finished job never runs its final static again -----------------------------------------------------------


class _Crash(Exception):
    pass


def _reference_trajectory(tmp_path: Path) -> list[list[float]]:
    module = _module()
    a = _attempt(tmp_path / "reference")
    assert _drive(module, a) == ("succeed",)
    return [r["y"] for r in _records(a.workdir)]


def _incomplete(workdir: Path) -> list[dict[str, Any]]:
    paths = sorted((workdir / "relax-driver" / "incomplete").glob("*/record.json"))
    return [json.loads(path.read_text(encoding="utf-8")) for path in paths]


def test_a_retry_after_the_terminal_state_was_merged_only_publishes(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "terminal")
    step = "prepare"
    # Run until the attempt that finished the last final static has merged its state and is about to publish.
    while True:
        outcome = _step(module, a, step)
        if outcome[:2] == ("advance", "publish"):
            break
        step = outcome[1]
    assert a.state["stage"] == "done"
    runs = len(_runs(a.workdir))
    report = _report(a.workdir)
    # The outcome was never published, and the report, relax/ and the static's files at the top were lost.
    for name in ("OUTCAR", "OSZICAR", "vasprun.xml"):
        (a.workdir / name).unlink()
    (a.workdir / "relax-driver" / "report.json").unlink()
    (a.workdir / "relax" / "CONTCAR").write_text("partial", encoding="utf-8")
    # The retried evaluate runs no VASP, restores the publication from the archive and publishes the same result.
    assert _step(module, a, "evaluate")[:2] == ("advance", "publish")
    assert len(_runs(a.workdir)) == runs
    assert _report(a.workdir) == report
    last = report["final_evaluation"]
    archived = a.workdir / "relax-driver" / "evaluations" / f"{last:03d}" / "OUTCAR"
    assert (a.workdir / "OUTCAR").read_text() == archived.read_text()
    assert (a.workdir / "relax" / "CONTCAR").read_text() != "partial"
    assert _step(module, a, "publish") == ("succeed",)


@pytest.mark.parametrize("when", ("after", "before"))
def test_an_attempt_stopped_around_its_snapshot_resumes_without_changing_the_trajectory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, when: str
) -> None:
    reference = _reference_trajectory(tmp_path)
    module = _module()
    a = _attempt(tmp_path / when)
    original = module.save_driver

    def save(attempt, driver, count):
        if count == 4 and not (attempt.workdir / "crashed").exists():
            (attempt.workdir / "crashed").touch()
            if when == "after":
                original(attempt, driver, count)
            raise _Crash
        original(attempt, driver, count)

    monkeypatch.setattr(module, "save_driver", save)
    with pytest.raises(_Crash):
        _drive(module, a)
    assert int(a.state["evaluations"]) == 3
    assert len(_runs(a.workdir)) == 4
    # Stopped after the snapshot: the retry continues after evaluation 4 without running it again. Stopped before:
    # the retry runs evaluation 4 again, from the WAVECAR the committed snapshot names.
    assert _drive(module, a, "evaluate") == ("succeed",)
    records = _records(a.workdir)
    runs = _runs(a.workdir)
    assert len(runs) == len(records) + (0 if when == "after" else 1)
    if when == "before":
        assert runs[3]["POSCAR"] == runs[4]["POSCAR"] and runs[3]["WAVECAR"] == runs[4]["WAVECAR"]
    assert [r["y"] for r in records] == reference


# -- R6: one corrected model ------------------------------------------------------------------------------------------


def _model_energy(cell: list[list[float]], frac: list[list[float]]) -> tuple[float, list, list]:
    """The fake VASP's model without noise: energy (eV), forces (eV/A) and stress (kbar)."""

    def mm(a, b):
        return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]

    tr = [list(row) for row in zip(*cell)]
    kc, ka = _MODEL["kc"], _MODEL["ka"]
    volume = abs(
        cell[0][0] * (cell[1][1] * cell[2][2] - cell[1][2] * cell[2][1])
        - cell[0][1] * (cell[1][0] * cell[2][2] - cell[1][2] * cell[2][0])
        + cell[0][2] * (cell[1][0] * cell[2][1] - cell[1][1] * cell[2][0])
    )
    delta = [[g - s for g, s in zip(gr, sr)] for gr, sr in zip(mm(cell, tr), _MODEL["gstar"])]
    energy = 0.25 * kc * sum(x * x for row in delta for x in row)
    dede = [[kc * x for x in row] for row in mm(mm(tr, delta), cell)]
    sstar = _MODEL["sstar"]
    ds = [((f0 - s0) - (f1 - s1) + 0.5) % 1.0 - 0.5 for f0, s0, f1, s1 in zip(frac[0], sstar[0], frac[1], sstar[1])]
    d = mm([ds], cell)[0]
    energy += sum(0.5 * ka * x * x for x in d)
    grad = [ka * x for x in d]
    for i in range(3):
        for j in range(3):
            dede[i][j] += 0.5 * (d[i] * grad[j] + d[j] * grad[i])
    stress = [[-x / volume * 1602.1766208 for x in row] for row in dede]
    return energy, [[-g for g in grad], grad], stress


def test_the_corrected_surface_has_one_energy_and_its_gradient_away_from_the_anchor(tmp_path: Path) -> None:
    module = _module()
    reference = module.read_poscar(_write(tmp_path / "POSCAR", _POSCAR))
    reference["V0"] = abs(module.det3(reference["cell"]))
    reference["L"] = reference["V0"] ** (1 / 3)
    reference["symmetry"] = module.symmetry_operations(reference["cell"], reference["frac"], ["Si", "Si"], 1e-3)
    # Altair's probe: paired force offsets of 0.1 eV/A along x and a 1 kbar diagonal stress offset, anchored at
    # the reference; evaluated at 1 % strain together with an atomic displacement.
    anchor = [0.0] * 12
    offset = {"forces": [[0.1, 0, 0], [-0.1, 0, 0]], "stress": [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]], "y": anchor}
    offset["gradient"] = module.symmetrize(
        reference, module.gradient(reference, anchor, offset["forces"], offset["stress"])
    )
    driver = {"reference": reference, "offset": offset}
    y = [0.01, -0.02, 0.005, -0.01, 0.02, -0.005, 0.01 * reference["L"], 0, 0, 0, 0, 0.002 * reference["L"]]

    def surface(point):
        cell, frac = module.geometry(reference, point)
        energy, forces, stress = _model_energy(cell, frac)
        record = {"kind": "coarse", "y": point, "energy": energy, "forces": forces, "stress": stress}
        return module.model(driver, record), record

    (_, g), record = surface(y)
    h = 1e-5
    for i in range(12):
        plus, minus = list(y), list(y)
        plus[i] += h
        minus[i] -= h
        derivative = (surface(plus)[0][0] - surface(minus)[0][0]) / (2 * h)
        assert derivative == pytest.approx(g[i], abs=1e-6)
    # The physical prediction differs there by about 0.1 eV/A x 1 % and is recorded.
    discrepancy = record["correction_model_discrepancy"]
    assert 0.0005 < discrepancy["force"] < 0.003 and discrepancy["stress"] < 0.05


# -- R7: symmetry -----------------------------------------------------------------------------------------------------


def test_a_symmetric_structure_keeps_its_symmetry_under_noisy_forces(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id, _ = _campaign(tmp_path / "cubic", monkeypatch, model=_CUBIC_MODEL, poscar=_CUBIC_POSCAR)
    kind, payload = _job(workspace, job_id)
    assert kind == "succeeded"
    report = _report(payload / "run")
    assert report["symmetry_operations"] == {"start": 96, "final": 96}
    final = report["evaluations"][-1]
    lengths = [math.sqrt(sum(x * x for x in row)) for row in final["cell"]]
    assert max(lengths) - min(lengths) < 1e-9 and lengths[0] == pytest.approx(4.0, abs=0.01)


def _seed(a: _Attempt, matrix: list[list[float]]) -> None:
    (a.payload / "files" / "hessian.json").write_text(json.dumps({"hessian": matrix}), encoding="utf-8")
    a.parameters["hessian_seed"] = "files/hessian.json"


@pytest.mark.parametrize(
    ("atom_x", "atom_y", "code"),
    ((70.0, 70.0, None), (7.0, 70.0, "vasp.hessian_seed_invalid"), (float("nan"), 70.0, "vasp.hessian_seed_invalid")),
)
def test_a_hessian_seed_must_be_finite_and_invariant_under_the_symmetry(
    tmp_path: Path, atom_x: float, atom_y: float, code: str | None
) -> None:
    module = _module()
    a = _attempt(tmp_path / "seed", poscar=_CUBIC_POSCAR, model=_CUBIC_MODEL)
    diagonal = [atom_x, atom_y, 70.0] * 2 + [5.0] * 6
    _seed(a, [[diagonal[i] if i == j else 0.0 for j in range(12)] for i in range(12)])
    outcome = _step(module, a, "prepare")
    if code is None:
        assert outcome[:2] == ("advance", "evaluate")
    else:
        assert outcome[:2] == ("fail", code)


def test_a_hessian_seed_learned_at_lower_symmetry_is_projected_when_asked(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "project", poscar=_CUBIC_POSCAR, model=_CUBIC_MODEL)
    diagonal = [7.0, 70.0, 70.0] * 2 + [5.0] * 6
    _seed(a, [[diagonal[i] if i == j else 0.0 for j in range(12)] for i in range(12)])
    a.parameters["hessian_seed_symmetrize"] = True
    assert _step(module, a, "prepare")[:2] == ("advance", "evaluate")
    hessian = json.loads((a.workdir / "relax-driver" / "driver-000.json").read_text())["optimizer"]["B"]
    # The cubic average of the atomic curvatures (7 + 70 + 70) / 3 on every axis; the strain block was invariant.
    assert [hessian[i][i] for i in range(6)] == pytest.approx([49.0] * 6)
    assert [hessian[i][i] for i in range(6, 12)] == pytest.approx([5.0] * 6)


# -- R8: trust bounds and cells ---------------------------------------------------------------------------------------


def test_the_trust_bounds_hold_in_the_current_cell(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "trust", poscar=_CUBIC_POSCAR)
    _step(module, a, "prepare")
    driver = json.loads((a.workdir / "relax-driver" / "driver-000.json").read_text())
    reference = driver["reference"]
    L = reference["L"]
    # The current point is 10 % isotropically expanded.
    current = [0.0] * 6 + [0.1 * L] * 3 + [0.0] * 3
    driver["optimizer"]["current"] = [current, 0.0, [0.0] * 12]
    opt = module.optimizer(a, driver)
    move = [0.03, 0, 0, -0.03, 0, 0] + [0.0] * 6
    # 0.03 A in the reference cell is 0.033 A in the current one: the step is cut to the 0.03 A bound.
    assert opt.step_scale(move, 1.0) == pytest.approx(1 / 1.1)
    strain = [0.0] * 6 + [0.01 * L] + [0.0] * 5
    # A 0.01 strain increment of the reference is 0.01 / 1.1 of the current cell: it may grow to the 0.01 bound.
    assert opt.step_scale(strain, 1.0) == pytest.approx(1.1)

    assert module.cell_problem(reference, [0.0] * 6 + [0.05 * L] * 3 + [0.0] * 3) is None
    assert "positive definite" in module.cell_problem(reference, [0.0] * 6 + [-1.0 * L] + [0.0] * 5)
    assert "volume" in module.cell_problem(reference, [0.0] * 6 + [0.5 * L] * 3 + [0.0] * 3)


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


# -- Altair recheck: B1, how a run that did not complete ended --------------------------------------------------------


@pytest.mark.parametrize("when", ("before", "after"))
@pytest.mark.parametrize("failure", ("exit", "nelm"))
def test_an_ordinary_failure_at_the_deadline_goes_to_the_remedy_ladder(tmp_path: Path, when: str, failure: str) -> None:
    # The first run takes 1 s and then exits with status 7 (no diagnosis), or reaches NELM; the attempt's deadline
    # falls just after it (3 s) or while it runs (0.3 s). Neither is an interruption.
    module = _module()
    model = {"sleep_on": {"0": 1.0}, **({"exit_on": {"0": 7}} if failure == "exit" else {"nelm_on": [0]})}
    a = _attempt(tmp_path / f"{failure}-{when}", model=model)
    _step(module, a, "prepare")
    a.deadline = time.time() + (3.0 if when == "before" else 0.3)
    outcome = _step(module, a, "evaluate")
    if failure == "exit":
        assert outcome[:2] == ("fail", "vasp.failed")
    else:
        assert outcome[0] == "retry" and a.state["remedies"] == 1
    driver = json.loads((a.workdir / "relax-driver" / "driver-000.json").read_text())
    assert driver["interruptions"] == [] and driver["slowest"] == {}
    (record,) = _incomplete(a.workdir)
    assert (record["evaluation"], record["termination"], record["handling"]) == (1, "exit", "remedy ladder")
    assert record["returncode"] == (7 if failure == "exit" else 0)
    if failure == "nelm":
        assert "electronic_nonconvergence" in [d["code"] for d in record["diagnostics"]]


def _signalled(module, monkeypatch: pytest.MonkeyPatch, termination: str = "signal_15") -> None:
    """Make every VASP run end as if the runner had received a signal while it ran."""

    import dataclasses

    original = module.run_vasp

    def run_vasp(*arguments, **options):
        report = original(*arguments, **options)
        process = dataclasses.replace(report.process, termination=termination, returncode=-15)
        return dataclasses.replace(report, process=process, classification="process_failure")

    monkeypatch.setattr(module, "run_vasp", run_vasp)


@pytest.mark.parametrize("left", (2.0, 3600.0, None))
def test_a_signal_is_an_interruption_only_at_the_deadline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, left: float | None
) -> None:
    module = _module()
    _signalled(module, monkeypatch)
    a = _attempt(tmp_path / f"signal-{left}")
    _step(module, a, "prepare")
    a.deadline = None if left is None else time.time() + left
    outcome = _step(module, a, "evaluate")
    driver = json.loads((a.workdir / "relax-driver" / "driver-000.json").read_text())
    (record,) = _incomplete(a.workdir)
    assert "remedies" not in a.state and record["termination"] == "signal_15"
    if left == 2.0:
        # At the deadline: deferred with no remedy, the stopped run's time as the estimate, and counted.
        assert outcome[:2] == ("advance", "evaluate") and record["handling"] == "deferred"
        assert [i["evaluation"] for i in driver["interruptions"]] == [1]
        assert driver["slowest"]["coarse:loose"] == pytest.approx(record["seconds"])
    else:
        # Far from the deadline, or without one: a cancellation or a stop from outside, for the operator.
        assert outcome[:2] == ("fail", "vasp.interrupted") and record["handling"] == "vasp.interrupted"
        assert driver["interruptions"] == []


def test_more_interruptions_than_allowed_fail_the_job(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    module = _module()
    _signalled(module, monkeypatch, "signal_15_killed")
    a = _attempt(tmp_path / "interruptions", parameters={"maximum_interruptions": 1})
    _step(module, a, "prepare")
    a.deadline = time.time() + 2.0
    assert _step(module, a, "evaluate")[:2] == ("advance", "evaluate")
    outcome = _step(module, a, "evaluate")
    assert outcome[:2] == ("fail", "vasp.interrupted_repeatedly")
    assert len(outcome[3]["interruptions"]) == 2
    # Both stopped runs of evaluation 1 are kept, the second next to the first.
    assert [(r["evaluation"], r["run"]) for r in _incomplete(a.workdir)] == [(1, 1), (1, 2)]
    assert len(_runs(a.workdir)) == 2


def test_the_runners_own_timeout_fails_without_a_remedy(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "timeout", model={"sleep_on": {"0": 30.0}}, parameters={"timeout": 0.5})
    _step(module, a, "prepare")
    a.deadline = time.time() + 3600
    outcome = _step(module, a, "evaluate")
    assert outcome[:2] == ("fail", "vasp.evaluation_timeout")
    (record,) = _incomplete(a.workdir)
    assert record["termination"].startswith("timeout") and record["handling"] == "vasp.evaluation_timeout"
    assert "remedies" not in a.state


# -- Altair recheck: N2, N3, N4 ---------------------------------------------------------------------------------------


def test_a_replayed_remedy_outside_the_electronic_ones_stops_before_staging(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "history")
    history = {
        "format": "httk-vasp-remedy-history",
        "format_version": 2,
        "attempts": {"electronic_nonconvergence": 1, "symmetry": 1},
        "events": [
            {"problem": "electronic_nonconvergence", "changes": [{"operation": "incar.AMIX", "value": 0.1}]},
            {"problem": "symmetry", "changes": [{"operation": "incar.ISYM", "value": 0}]},
        ],
    }
    (a.payload / ".httk-job" / "vasp-remedies.json").write_text(json.dumps(history), encoding="utf-8")
    _step(module, a, "prepare")
    outcome = _step(module, a, "evaluate")
    assert outcome[:2] == ("fail", "vasp.unsupported_remedy")
    assert "incar.ISYM" in outcome[2] and "incar.AMIX" not in outcome[2]
    assert _runs(a.workdir) == [] and not (a.workdir / "INCAR").exists()


@pytest.mark.parametrize("cap", (1, 2))
def test_the_coarse_cap_holds_at_the_ediff_switch(tmp_path: Path, cap: int) -> None:
    # Switch limits so wide that the first point already qualifies: the tight repeat is due after evaluation 1.
    module = _module()
    parameters = {"maximum_evaluations": cap, "switch_force": 1e3, "switch_stress": 1e6}
    a = _attempt(tmp_path / f"cap-{cap}", parameters=parameters)
    outcome = _drive(module, a)
    assert outcome[:2] == ("fail", "vasp.relax_not_converged")
    assert outcome[3]["next"] == ("tighten" if cap == 1 else "step")
    assert len(_runs(a.workdir)) == cap
    if cap == 2:
        assert float(_runs(a.workdir)[1]["tags"]["EDIFF"]) == 1e-8


def test_the_correction_discrepancy_is_in_physical_units(tmp_path: Path) -> None:
    module = _module()
    reference = module.read_poscar(_write(tmp_path / "POSCAR", _POSCAR))
    reference["V0"] = abs(module.det3(reference["cell"]))
    reference["L"] = reference["V0"] ** (1 / 3)
    L = reference["L"]
    # Altair's case: 10 % isotropic contraction, a generalized force of 0.00029 is 0.00029 / 0.9 eV/A.
    contracted = [0.0] * 6 + [-0.1 * L] * 3 + [0.0] * 3
    g = [0.00029, 0, 0, -0.00029, 0, 0] + [0.0] * 6
    assert module.physical_size(reference, contracted, g)[0] == pytest.approx(0.00029 / 0.9, rel=1e-12)
    # Forces and a general stress, through gradient() and back, at a strained and sheared cell.
    y = [0.01, -0.02, 0.005, -0.01, 0.02, -0.005, -0.08 * L, 0.03 * L, 0.05 * L, 0.02 * L, -0.04 * L, 0.03 * L]
    forces = [[0.1, -0.04, 0.02], [-0.1, 0.04, -0.02]]
    stress = [[1.0, 0.3, -0.7], [0.3, -2.5, 0.4], [-0.7, 0.4, 0.8]]
    g = module.gradient(reference, y, forces, stress)
    assert module.physical_size(reference, y, g) == pytest.approx(module.maxima(forces, stress), rel=1e-10)


# -- O1 (J-0210): the acceptance test when the energy and the gradients disagree -------------------------------------


def _biased(y: list[float], bias: float) -> tuple[float, list[float]]:
    """E = y^2/2 + bias y with the gradient of y^2/2 alone: an energy its gradient does not belong to."""

    return 0.5 * y[0] ** 2 + bias * y[0], [y[0]]


def test_the_gradient_acceptance_judges_a_step_by_the_rise_its_gradients_predict(tmp_path: Path) -> None:
    module = _module()

    def told(acceptance: str, ceiling: float | None, to: float, bias: float) -> tuple[bool, dict]:
        opt = module.TrustBFGS([[1.0]], lambda p, t: 1.0, 2e-4, ceiling)
        opt.tell([0.1], *_biased([0.1], bias), acceptance)
        accepted = opt.tell([to], *_biased([to], bias), acceptance)
        return accepted, opt.history[-1]

    # To the gradient's minimum on a biased energy: the gradients predict a fall of 0.005, the energy rises by 0.005.
    accepted, entry = told("energy", None, 0.0, -0.1)
    assert accepted is False and entry["trust_factor"] == 0.5
    assert entry["rise"] == pytest.approx(0.005) and entry["predicted_rise"] == pytest.approx(-0.005)
    accepted, entry = told("gradient", None, 0.0, -0.1)
    assert accepted is True and entry["trust_factor"] == 1.0 and entry["acceptance"] == "gradient"
    # A recorded rise above the ceiling rejects the step all the same.
    assert told("gradient", 2e-3, 0.0, -0.1)[0] is False
    assert told("gradient", 1e-2, 0.0, -0.1)[0] is True
    # An overshoot on a consistent quadratic rises on both counts, by the same amount (the trapezoid is exact there).
    accepted, entry = told("gradient", None, -0.3, 0.0)
    assert accepted is False and entry["predicted_rise"] == pytest.approx(entry["rise"], rel=1e-12) and entry["rise"] > 0
    with pytest.raises(ValueError):
        told("forces", None, 0.0, 0.0)


@pytest.mark.parametrize("acceptance", ("energy", "gradient_tight"))
def test_an_energy_that_disagrees_with_the_stress_is_judged_by_the_gradients(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, acceptance: str
) -> None:
    # The model's energy falls by 20 kbar times the volume beyond the surface its stress belongs to. Judged by the
    # energy, the steps that shrink the cell towards the stress zero rise and are rejected, past a cap of 16 (20
    # coarse evaluations without the cap); judged by the gradients from the EDIFF switch on, none is.
    workspace, job_id, _ = _campaign(
        tmp_path / acceptance,
        monkeypatch,
        model={"pulay": 20.0},
        parameters={"acceptance": acceptance, "maximum_evaluations": 16},
    )
    kind, payload = _job(workspace, job_id)
    coarse = [r for r in _records(payload / "run") if r["kind"] == "coarse"]
    steps = [r for r in coarse if "energy_gradient_mismatch" in r]
    assert steps
    kbar = 1 / 1602.1766208
    for r in steps:
        assert r["acceptance"] == ("energy" if acceptance == "energy" or r["level"] == "loose" else "gradient")
        # The mismatch is the model's term, up to the loose evaluations' noise.
        assert r["energy_gradient_mismatch"] == pytest.approx(-20.0 * kbar * r["volume_change"], abs=5e-5)
    if acceptance == "energy":
        assert kind == "failed"
        failure = _failure(payload)
        assert failure["code"] == "vasp.relax_not_converged"
        assert any(r["accepted"] is False for r in coarse if r["level"] != "loose")
        fits = failure["details"]["energy_gradient_consistency"]
    else:
        assert kind == "succeeded"
        report = _report(payload / "run")
        assert report["within_limits"] is True
        assert all(r["accepted"] for r in coarse if r["level"] != "loose")
        assert len(coarse) <= 13
        fits = report["energy_gradient_consistency"]
    assert fits["tight"]["sigma_kbar"] == pytest.approx(20.0, rel=1e-3)
    assert fits["tight"]["rms_after_fit_eV"] < 1e-6 < fits["tight"]["rms_mismatch_eV"]


def test_an_unknown_acceptance_fails_before_anything_runs(tmp_path: Path) -> None:
    module = _module()
    a = _attempt(tmp_path / "acceptance", parameters={"acceptance": "forces"})
    outcome = _step(module, a, "prepare")
    assert outcome[:2] == ("fail", "vasp.input_invalid") and "'forces'" in outcome[2]
    assert not (a.workdir / "relax-driver").exists()


def test_the_tight_repeat_of_the_switch_point_does_not_set_the_tight_step_estimate(tmp_path: Path) -> None:
    # J-0210: the repeat took 149 s (4 SCF steps from a WAVECAR at its own geometry), the next tight step 441 s.
    module = _module()
    a = _attempt(tmp_path / "estimate", parameters={"expected_seconds": 100.0})
    assert _drive(module, a) == ("succeed",)
    switch = _report(a.workdir)["switched_at"]
    records = _report(a.workdir)["evaluations"]
    assert records[switch]["purpose"] == "tighten" and records[switch + 1]["purpose"] == "step"

    def snapshot(index: int) -> dict:
        return json.loads(module.driver_path(a, index).read_text(encoding="utf-8"))

    after_repeat = snapshot(switch + 1)
    assert "coarse:tight" not in after_repeat["slowest"]
    assert after_repeat["slowest"]["coarse:tighten"] == pytest.approx(records[switch]["seconds"])
    loose = after_repeat["slowest"]["coarse:loose"]
    repeat = after_repeat["slowest"]["coarse:tighten"]
    assert module.expected_seconds(a, after_repeat, "coarse:tight") == pytest.approx(max(1.35 * loose, repeat))
    assert snapshot(switch + 2)["slowest"]["coarse:tight"] == pytest.approx(records[switch + 1]["seconds"])


@pytest.mark.parametrize(
    "repeat, expected",
    # Altair's numbers: slowest loose 416 s; a cheap repeat (J-0210: 149 s) leaves 1.35 x 416, a slow one (900 s,
    # hypothetical) is the estimate; no loose time and no expected_seconds leave the repeat alone.
    ((None, 561.6), (149.0, 561.6), (900.0, 900.0)),
)
def test_a_slow_tight_repeat_is_a_lower_bound_on_the_first_tight_step(
    tmp_path: Path, repeat: float | None, expected: float
) -> None:
    module = _module()
    a = _attempt(tmp_path / "bound", parameters={})
    slowest = {"coarse:loose": 416.0} | ({} if repeat is None else {"coarse:tighten": repeat})
    driver = {"slowest": slowest, "overhead": {"staging": 10.0, "archive": 20.0}}
    plan = {"kind": "coarse", "purpose": "step", "ediff": 1e-8}
    assert module.expected_seconds(a, driver, "coarse:tight") == pytest.approx(expected)
    assert module.time_needed(a, driver, plan) == pytest.approx(1.4 * expected + 30.0)
    assert module.time_needed(a, driver, plan, staged=True) == pytest.approx(1.4 * expected + 20.0)
    # The first final falls back on the same estimate; an ordinary tight step, once run, replaces it.
    assert module.expected_seconds(a, driver, "final") == pytest.approx(8.0 * expected)
    driver["slowest"]["coarse:tight"] = 441.0
    assert module.expected_seconds(a, driver, "coarse:tight") == 441.0
    if repeat is not None:
        # Without a loose time or expected_seconds the repeat is the only estimate.
        assert module.expected_seconds(a, {"slowest": {"coarse:tighten": repeat}}, "coarse:tight") == repeat


@pytest.mark.parametrize("repeat", (1.0, 30.0))
def test_the_first_tight_step_waits_for_the_time_a_slow_repeat_took(tmp_path: Path, repeat: float) -> None:
    # Slowest loose 10 s, overheads 1.5 s (staging) and 2 s (archive), 40 s left in the attempt. After a cheap repeat
    # the step needs 1.4 x 13.5 + 3.5 = 22.4 s and runs; after a 30 s one it needs 1.4 x 30 + 3.5 = 45.5 s and is
    # deferred with that mintime, rounded up to the second.
    module = _module()
    a = _attempt(tmp_path / f"wait-{repeat:g}", parameters={}, deadline=time.time() + 40)
    driver = {"slowest": {"coarse:loose": 10.0, "coarse:tighten": repeat}, "overhead": {"staging": 1.5, "archive": 2.0}}
    plan = {"kind": "coarse", "purpose": "step", "ediff": 1e-8}
    needed = module.time_needed(a, driver, plan)
    assert needed == pytest.approx(22.4 if repeat == 1.0 else 45.5)
    runs = module.guard(a, needed, plan, "before staging")
    if repeat == 1.0:
        assert runs is True and a.outcome is None
    else:
        assert runs is False and a.outcome == ("advance", "evaluate", {"mintime": "0:46"})


def test_an_interrupted_tight_repeat_keeps_its_time_apart_and_bounds_the_tight_estimate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _module()
    a = _attempt(tmp_path / "cut", parameters={}, deadline=time.time() + 1.0)
    driver = {
        "slowest": {"coarse:loose": 10.0},
        "overhead": {"staging": 0.0, "archive": 0.0},
        "next": {"kind": "coarse", "purpose": "tighten", "ediff": 1e-8, "y": []},
        "interruptions": [],
    }
    monkeypatch.setattr(module, "save_driver", lambda *arguments: None)
    module.interrupted(a, driver, 7, 50.5, "signal_15")
    assert "coarse:tight" not in driver["slowest"] and driver["slowest"]["coarse:tighten"] == 50.5
    assert module.expected_seconds(a, driver, "coarse:tight") == 50.5
    # 1.4 x 50.5 = 70.7 s.
    assert a.outcome == ("advance", "evaluate", {"mintime": "1:11"})


# -- O1 review (Altair, 2026-10-18): the failure's fit, the loose stall, the diagnostic without a volume change -----


@pytest.mark.parametrize("cap", (2, 7))
def test_the_cap_failure_fits_every_evaluation_it_counts(tmp_path: Path, cap: int) -> None:
    # R1: the failure was decided before the last record was written, and its fit left that evaluation out.
    module = _module()
    a = _attempt(tmp_path / f"fit-{cap}", model={"pulay": 20.0}, parameters={"maximum_evaluations": cap})
    outcome = _drive(module, a)
    assert outcome[:2] == ("fail", "vasp.relax_not_converged")
    records = _records(a.workdir)
    assert len(records) == cap == outcome[3]["evaluations"]
    fits = outcome[3]["energy_gradient_consistency"]
    assert fits == module.consistency(records)
    steps = [r for r in records if "energy_gradient_mismatch" in r]
    fitted = [pair[1] for fit in fits.values() for pair in fit["pairs"]]
    assert sorted(fitted) == [r["evaluation"] for r in steps]
    if cap == 2:
        # Altair's case: the step that reached the cap was the only one, and the fit was empty.
        assert fitted == [2]


@pytest.mark.parametrize("pulay, stall", ((5.0, 0), (5.0, 3), (20.0, 3)))
def test_repeated_loose_rejections_tighten_the_ediff_early(tmp_path: Path, pulay: float, stall: int) -> None:
    # R2. A cubic cell 2.5 % too large, no noise, no final-mesh change: only the volume moves, and the energy falls
    # *pulay* kbar x dV faster than the stress says. Judged by the energy at the loose EDIFF, the steps that shrink the
    # cell are rejected before the switch limits are reached, and without the stall switch the job never leaves the
    # loose level. With it, three rejections in a row tighten the EDIFF at the accepted point and the gradients judge
    # from there. At 20 kbar the rise ceiling (2 meV per step, per cell: not size-extensive and no bound on the drift
    # over many steps) still rejects every step large enough to move the volume, and the cap of 24 is reached: a
    # documented limit, not a setting to raise (the 20 kbar is the test model's, not a measured Sn value; J-0210's
    # fit is 2.15 kbar, which converges with or without the switch).
    module = _module()
    model = {**_CUBIC_MODEL, "pulay": pulay, "noise_force": 0.0, "noise_stress": 0.0, "final_metric_scale": 1.0}
    parameters = {"maximum_evaluations": 24, "stall_rejections": stall}
    a = _attempt(tmp_path / f"stall-{pulay:g}-{stall}", poscar=_CUBIC_POSCAR, model=model, parameters=parameters)
    outcome = _drive(module, a)
    coarse = [r for r in _records(a.workdir) if r["kind"] == "coarse"]
    loose = [r for r in coarse if r["level"] == "loose"]
    tight = [r for r in coarse if r["level"] == "tight" and r["purpose"] == "step"]
    if stall == 0:
        assert outcome[:2] == ("fail", "vasp.relax_not_converged") and outcome[3]["level"] == "loose"
        assert len(coarse) == 24 and not tight
        return
    # The last snapshot (a failed job writes none after its last evaluation).
    driver = json.loads(sorted((a.workdir / "relax-driver").glob("driver-*.json"))[-1].read_text(encoding="utf-8"))
    assert driver["switch_reason"] == "stall"
    # The switch came after exactly three rejections in a row, at a point outside the switch limits.
    assert [r["accepted"] for r in loose[-3:]] == [False] * 3
    switch = next(r for r in coarse if r["evaluation"] == driver["switched_at"])
    limits = module.DEFAULTS
    assert switch["max_force"] > limits["switch_force"] or switch["max_stress"] > limits["switch_stress"]
    assert all(r["acceptance"] == "gradient" for r in tight)
    if pulay == 5.0:
        assert outcome == ("succeed",) and _report(a.workdir)["within_limits"] is True
        assert _report(a.workdir)["switch_reason"] == "stall"
        assert len(coarse) <= 12
    else:
        assert outcome[:2] == ("fail", "vasp.relax_not_converged") and outcome[3]["level"] == "tight"
        rejected = [r for r in tight if not r["accepted"]]
        # Every tight rejection is the ceiling's: the gradients predicted no rise.
        assert rejected and all(r["predicted_rise"] <= 2e-4 < 2e-3 < r["energy_rise"] for r in rejected)


def _step_record(evaluation: int, dv: float, mismatch: float, volume: float = 64.0) -> dict:
    side = volume ** (1 / 3)
    cell = [[side, 0.0, 0.0], [0.0, side, 0.0], [0.0, 0.0, side]]
    return {
        "evaluation": evaluation,
        "compared_with": evaluation - 1,
        "level": "tight",
        "cell": cell,
        "volume_change": dv,
        "energy_gradient_mismatch": mismatch,
    }


def test_the_consistency_report_keeps_a_mismatch_no_volume_change_explains(tmp_path: Path) -> None:
    # R4. Altair's case: a tight step at fixed volume whose energy rose 1.95 meV beyond its gradients' prediction.
    module = _module()
    sigma = 2.15 * module.KBAR
    fit = module.consistency([_step_record(5, 0.0, 0.00195)])["tight"]
    assert (fit["steps"], fit["volume_steps"], fit["identifiable"]) == (1, 0, False)
    assert fit["rms_mismatch_eV"] == pytest.approx(0.00195) and fit["max_abs_mismatch_eV"] == pytest.approx(0.00195)
    assert fit["sigma_kbar"] is None and fit["rms_after_fit_eV"] is None and fit["pairs"] == [[4, 5]]
    # Below VOLUME_RESOLUTION (1e-4 of 64 A^3) a step does not count; one volume step is fitted exactly by
    # construction, so it is not enough either.
    fit = module.consistency([_step_record(5, 0.005, 0.0), _step_record(6, 1.0, -sigma)])["tight"]
    assert (fit["volume_steps"], fit["identifiable"], fit["sigma_kbar"]) == (1, False, None)
    assert fit["sum_dV2_A6"] == pytest.approx(1.000025)
    # Two isotropic volume steps give the stress back, with nothing left over.
    fit = module.consistency([_step_record(5, 1.0, -sigma), _step_record(6, -0.5, 0.5 * sigma)])["tight"]
    assert fit["identifiable"] is True and fit["sigma_kbar"] == pytest.approx(2.15)
    assert fit["rms_after_fit_eV"] == pytest.approx(0.0, abs=1e-15)
    # A mismatch that does not follow the volume leaves a residual as large as itself: no isotropic stress explains it.
    fit = module.consistency([_step_record(5, 1.0, 0.002), _step_record(6, -1.0, 0.002)])["tight"]
    assert fit["identifiable"] is True and fit["sigma_kbar"] == pytest.approx(0.0, abs=1e-12)
    assert fit["rms_after_fit_eV"] == pytest.approx(0.002) == fit["rms_mismatch_eV"]
