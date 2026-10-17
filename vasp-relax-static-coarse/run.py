#!/usr/bin/env python3
"""Relax on a coarse k-point mesh with an optimizer in Python, then run one final-mesh static.

A copy of ``vasp-relax-static`` edited as a whole: the staged IBRION relaxation is
replaced by a driver loop in which VASP is only the energy, force and stress
engine. Every evaluation is one fresh single point (``NSW = 0``) at a geometry the
optimizer chose, restarting from the previous evaluation's WAVECAR, and the
optimizer (``relax_optimizer.py`` next to this file: a damped BFGS with a trust
radius, standard library only) proposes the next geometry. The relaxation runs
inside one job for as long as it takes; only the time guard below ever pauses it.

The stages, each a run of the one ``evaluate`` step:

* **coarse, loose** — evaluations on the coarse mesh (``coarse_mesh``, else
  ``coarse_kpoint_density``) from ``INCAR.coarse`` with ``coarse_incar_tags``, at
  ``EDIFF = ediff_loose``, until the accepted point has a maximum force at most
  ``switch_force`` and a maximum stress component at most ``switch_stress``;
* **coarse, tight** — the same at ``EDIFF = ediff_tight``. The switch point is
  first evaluated again at ``ediff_tight``, and that evaluation *replaces* the
  optimizer's accepted point (energy and gradient) without an acceptance test or
  a curvature update: energies, forces and stress at two EDIFF levels are never
  compared. Its ``ediff_shift`` records what the tighter EDIFF changed at one
  geometry. With ``offset_probe`` a final-mesh evaluation follows at that point,
  so the coarse-to-final offset is measured at two geometries of the path. The
  stage ends when an accepted point evaluated at ``ediff_tight`` is within
  ``force_limit`` and ``stress_limit``;
* **final** — one single point of that geometry on the final mesh
  (``final_mesh``, else ``kpoint_density``) from ``INCAR.final`` with
  ``final_incar_tags``, at ``EDIFF = final_ediff``. Its forces and stress minus
  those of the tight coarse evaluation at the same geometry (both checked) are
  the offset. If the final static misses the limits, the coarse relaxation
  continues on a corrected surface (**coarse, corrected**), and a second final
  static ends the job (``maximum_corrections``, default 1).

The corrected surface is one model: the coarse energy plus the linear term
``g_offset . (y - y_anchor)``, whose gradient is the constant ``g_offset``, the
offset's gradient in the optimizer's coordinates at the geometry where it was
measured. The optimizer sees exactly that energy and gradient. The stop test
uses the physical prediction, coarse forces and stress plus the measured
offset; the two differ away from the anchor, and every corrected evaluation
records that difference (``correction_model_discrepancy``) as the largest
force (eV/A) and stress component (kbar) at its own geometry, with the largest
compared with ``correction_bound_fraction`` of the limits in the report. It is
a diagnostic: only the final static's own forces and stress decide.

The job succeeds once its last final static has run, whether or not that static
is within the limits; ``within_limits`` in the job state and in
``relax-driver/report.json`` says which, and ``collect.py`` harvests the
structure and energy only when it is true. More than ``maximum_evaluations``
coarse evaluations fail the job with ``vasp.relax_not_converged``, also when
the budget is used up just as the EDIFF is tightened (the tight repeat of the
switch point is never skipped), so the cap is strict.

Time guard: an evaluation starts only when the attempt's time left
(``Attempt.deadline``: its ``maxtime`` or the manager's drain start) covers
``time_factor`` times the expected VASP time of that evaluation plus the
largest staging and archiving overhead measured so far. It is checked before
the inputs are staged and again right before VASP starts. The expected time is
the slowest evaluation so far of the same kind (coarse at the loose EDIFF,
coarse at the tight one, final); before the first of a kind it is
``expected_seconds`` for coarse-loose, ``tight_cost_factor`` times the
coarse-loose one for coarse-tight, and ``final_cost_factor`` times the slowest
coarse one for final. Without ``expected_seconds`` the first evaluation runs
unguarded, and the log says so. When the attempt does not have the time, the
evaluation is published again with that time as its ``mintime``, so a manager
without enough time left does not claim it, and a later manager resumes the
relaxation where it stopped. When the time needed exceeds the attempt's own
``maxtime``, no attempt could ever run it: the job fails with
``vasp.insufficient_time`` (``needed_seconds`` and ``maxtime_seconds`` in the
details), and is continued once its ``maxtime`` is raised.

A VASP run that does not complete is classified by how its process ended
(``termination`` of the supervised process), never by the clock alone:

* ended by a signal to the runner (``signal_15``, also ``…_killed``) at the
  attempt's deadline (``Attempt.deadline`` minus a few seconds or later): the
  manager stopped it at its maxtime or drain. The evaluation is deferred, with
  no remedy, and the stopped run's time becomes its kind's estimate. At most
  ``maximum_interruptions`` such deferrals per job; one more fails it with
  ``vasp.interrupted_repeatedly``;
* ended by a signal well before the deadline (a cancellation, or a stop from
  outside httk): the job fails with ``vasp.interrupted`` for the operator. A
  cancelled job is committed as cancelled by the manager whatever it publishes;
* ended by the runner's own ``timeout`` (a hang guard): the job fails with
  ``vasp.evaluation_timeout``, with no remedy;
* any other ending (a nonzero exit, a diagnosed stop, nonconvergence), however
  close to the deadline: the remedy ladder below.

Every run that did not complete keeps its files and a ``record.json`` (how it
ended, how it was handled, its diagnostics) under
``relax-driver/incomplete/NNN-KK/`` (evaluation NNN, its KK-th incomplete run),
so the next run of that evaluation never overwrites the evidence. A runner
killed outright (SIGKILL), or stopped by a signal outside the VASP run (while
it parses, archives or writes a snapshot), publishes nothing; the manager
records the attempt's failure and the snapshots make a retry safe.

Every evaluation is archived under ``relax-driver/evaluations/NNN/`` with its
inputs, outputs and ``record.json`` (stage, mesh, EDIFF, evaluation index,
energy, forces, stress, SCF steps, seconds, whether a WAVECAR restart was
requested and whether VASP's output shows it was read, remedies in force,
offset applied, accepted or not); the optimizer's state after evaluation NNN is
``relax-driver/driver-NNN.json``, which also names the WAVECAR each mesh kind
restarts from (``relax-driver/<kind>/WAVECAR-NNN``, written by evaluation NNN,
with the inputs it was written for). A retried attempt continues from exactly
the evaluations the job state counts, or from a later snapshot an interrupted
attempt had already written; a WAVECAR no snapshot names is deleted. Once the
last final static is done the snapshot says so, and a retry only publishes:
it never runs that static again. ``relax/`` holds the coarse evaluation at the
final geometry, the final static's files are at the top of the workdir, and
``relax-driver/report.json`` gathers every record and offset; that is where
``collect.py`` reads.

The reviewed remedy ladder is kept for the remedies that change how VASP
reaches the ground state but not the ground state itself (charge mixing,
algorithm, parallelization): an evaluation that fails in a way it recognizes
gets exactly one such remedy and is retried at the same geometry, and the
remedies stay applied in every later evaluation. A remedy that would change the
energy surface or the structure (k-points, smearing, symmetry, LREAL, bands,
EDIFF, the lattice) stops the job with ``vasp.unsupported_remedy`` rather than
compare evaluations of different surfaces.

Symmetry: the symmetry operations of the starting structure are found once
(``symmetry_tolerance``, A) and every gradient is projected onto the
coordinates they leave invariant, so the path keeps the starting symmetry
whatever VASP's own symmetrization does; a ``hessian_seed`` must be invariant
under them. The report gives the number of operations at the start and at the
final geometry.

The job inputs and parameters are documented in this repository's README.
"""

import itertools
import json
import math
import shlex
import shutil
import sys
import time
from pathlib import Path
from xml.etree import ElementTree

from httk.codes.vasp import (
    VaspPreparationOptions,
    apply_vasp_remedy,
    clean_vasp_outputs,
    job_remedy_history_path,
    last_oszicar_energy,
    normalize_poscar_handedness,
    plan_vasp_remedy,
    prepare_vasp_inputs,
    read_incar,
    reapply_vasp_remedies,
    run_vasp,
    update_incar,
    validate_vasp_workdir,
    write_automatic_kpoints,
)
from httk.workflow import Attempt, Runner

sys.path.insert(0, str(Path(__file__).resolve().parent))
from relax_optimizer import TrustBFGS

COLLECT = "INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# What relax-driver/evaluations/NNN/ keeps of every evaluation.
EVALUATION_ARCHIVE = (
    "POSCAR",
    "INCAR",
    "KPOINTS",
    "OUTCAR",
    "OSZICAR",
    "CONTCAR",
    "vasprun.xml",
    "vasp-run-report.json",
    "vasp.out",
    "vasp.err",
)
TEMPLATES = Path(__file__).resolve().parent / "templates"
DRIVER = "relax-driver"
# Where the coarse evaluation at the final geometry is archived for collect.py.
RELAX_ARCHIVE = "relax"
# eV/A^3 per kbar.
KBAR = 1 / 1602.1766208
SQRT2 = math.sqrt(2.0)
# Remedies that change how VASP reaches the electronic ground state, not the ground state: the only ones the driver
# applies, since every other one changes the surface the optimizer has learned (or the structure it controls).
ELECTRONIC_REMEDIES = frozenset(
    f"incar.{tag}" for tag in ("ICHARG", "AMIX", "BMIX", "AMIN", "ALGO", "TIME", "NPAR", "NEDOS")
)
# The INCAR tags a WAVECAR depends on besides the k-point mesh; a cached one is staged only for the same values.
WAVECAR_TAGS = ("ENCUT", "PREC", "NBANDS", "ISPIN", "LNONCOLLINEAR", "LSORBIT", "ISYM", "LHFCALC", "NGX", "NGY", "NGZ")
# VASP's stdout when it reads a WAVECAR, and when it does not (VASP 6.6.0, jobs 415f7cbc and c9adca60 of
# endpoint-pbe-hse); anything else is reported as unknown.
WAVECAR_READ = "the WAVECAR file was read successfully"
WAVECAR_NOT_READ = "WAVECAR not read"
# A run stopped by a signal this close to the attempt's deadline, or later, was stopped by the manager at that deadline:
# the deadline is a whole second and the report arrives after VASP has exited.
DEADLINE_SLACK = 5.0
# The defaults of every driver parameter; the README gives their origin (the Sn/Ge retro of endpoint-pbe-hse).
DEFAULTS = {
    "coarse_mesh": "",
    "final_mesh": "",
    "coarse_kpoint_density": 10.0,
    "kpoint_density": 20.0,
    "centering": "Gamma",
    "coarse_incar_tags": {"PRECFOCK": "Fast"},
    "final_incar_tags": {"PRECFOCK": "Normal"},
    "ediff_loose": 1e-6,
    "ediff_tight": 1e-8,
    "final_ediff": 1e-8,
    "switch_force": 0.02,
    "switch_stress": 1.0,
    "force_limit": 0.003,
    "stress_limit": 0.2,
    "time_factor": 1.4,
    "expected_seconds": 0.0,
    "tight_cost_factor": 2.0,
    "final_cost_factor": 8.0,
    "maximum_evaluations": 60,
    "maximum_corrections": 1,
    "maximum_interruptions": 3,
    "trust_atom": 0.03,
    "trust_strain": 0.01,
    "rise_tolerance": 2e-4,
    "hessian_atom": 70.0,
    "hessian_bulk_modulus": 50.0,
    "hessian_seed": "",
    "hessian_seed_symmetrize": False,
    "symmetry_tolerance": 1e-3,
    "correction_bound_fraction": 0.1,
    "offset_probe": False,
    "coarse_wavecar": "files/WAVECAR",
    "final_wavecar": "files/WAVECAR.final",
}

run = Runner("vasp.relax-static-coarse")


def parameter(a: Attempt, name: str):
    """Return a driver parameter of the job, or its default."""

    value = a.parameter(name, DEFAULTS[name])
    return DEFAULTS[name] if value is None else value


def flag(value: object) -> bool:
    """Read a boolean job parameter, which may arrive as a string from the command line."""

    return value is True or str(value).strip().lower() in ("true", "1", "yes")


# -- 3x3 and vector algebra, standard library only -------------------------------------------------------------------


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def transpose(m: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*m)]


def det3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def inv3(m: list[list[float]]) -> list[list[float]]:
    d = det3(m)
    cof = [
        [
            m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3]
            - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]
            for j in range(3)
        ]
        for i in range(3)
    ]
    return [[cof[i][j] / d for j in range(3)] for i in range(3)]


def voigt(t: list[list[float]]) -> list[float]:
    """A symmetric tensor as an orthonormal 6-vector (off-diagonals weighted by sqrt 2)."""

    return [t[0][0], t[1][1], t[2][2], SQRT2 * t[1][2], SQRT2 * t[0][2], SQRT2 * t[0][1]]


def tensor(w: list[float]) -> list[list[float]]:
    return [
        [w[0], w[5] / SQRT2, w[4] / SQRT2],
        [w[5] / SQRT2, w[1], w[3] / SQRT2],
        [w[4] / SQRT2, w[3] / SQRT2, w[2]],
    ]


def norm(v: list[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


# -- structures and coordinates ---------------------------------------------------------------------------------------
#
# The optimizer's coordinates are y = [u, L w]: u the Cartesian displacements of the atoms in the reference cell
# (3N, in A) and w the symmetric strain eps of the cell as a sqrt-2-weighted Voigt vector, scaled by L = V0^(1/3) so
# that both parts are lengths. The cell is C = C0 (I + eps) (rows are lattice vectors) and the fractional coordinates
# are s = s0 + u C0^-1, so an atom's physical displacement is u (I + eps). The gradient is dE/du = -F (I + eps) and
# dE/deps = -V sym(sigma (I + eps)^-1), with VASP's stress (positive = compressive) in eV/A^3; the net force is
# removed so the structure does not drift, and the gradient is projected onto the coordinates the starting
# structure's symmetry operations leave invariant.


def read_poscar(path: Path) -> dict:
    """Read a VASP-5 POSCAR into the comment, cell (A), species, counts and fractional coordinates."""

    lines = path.read_text(encoding="utf-8").splitlines()
    scale = float(lines[1].split()[0])
    if scale <= 0:
        raise ValueError("a POSCAR scale that gives the volume (negative) is not supported")
    cell = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    mode = lines[7].strip()
    if mode[:1] in "sS":
        raise ValueError("selective dynamics is not supported by the driver")
    n = sum(counts)
    coordinates = [[float(x) for x in lines[8 + i].split()[:3]] for i in range(n)]
    if mode[:1] in "cCkK":
        inverse = inv3(cell)
        coordinates = matmul([[scale * x for x in row] for row in coordinates], inverse)
    return {"comment": lines[0], "cell": cell, "species": species, "counts": counts, "frac": coordinates}


def write_poscar(path: Path, reference: dict, cell: list[list[float]], frac: list[list[float]]) -> None:
    lines = [reference["comment"], "1.0"]
    lines += [" ".join(f"{x:.16f}" for x in row) for row in cell]
    lines += [" ".join(reference["species"]), " ".join(str(c) for c in reference["counts"]), "Direct"]
    lines += [" ".join(f"{x:.16f}" for x in row) for row in frac]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def strain_of(reference: dict, y: list[float]) -> list[list[float]]:
    n3 = 3 * len(reference["frac"])
    return tensor([v / reference["L"] for v in y[n3:]])


def deformation(reference: dict, y: list[float]) -> list[list[float]]:
    """I + eps of coordinates y."""

    eps = strain_of(reference, y)
    return [[(1.0 if i == j else 0.0) + eps[i][j] for j in range(3)] for i in range(3)]


def geometry(reference: dict, y: list[float]) -> tuple[list[list[float]], list[list[float]]]:
    """The cell and fractional coordinates of coordinates y."""

    cell = matmul(reference["cell"], deformation(reference, y))
    inverse = inv3(reference["cell"])
    frac = []
    for a, s0 in enumerate(reference["frac"]):
        du = matmul([y[3 * a : 3 * a + 3]], inverse)[0]
        frac.append([s0[k] + du[k] for k in range(3)])
    return cell, frac


def cell_problem(reference: dict, y: list[float]) -> str | None:
    """Why the cell of coordinates y must not be written (singular, inverted or far deformed), or None."""

    d = deformation(reference, y)
    if not all(math.isfinite(x) for row in d for x in row):
        return "the deformation is not finite"
    # I + eps is symmetric: positive definite exactly when its leading principal minors are positive.
    minors = (d[0][0], d[0][0] * d[1][1] - d[0][1] * d[1][0], det3(d))
    if min(minors) <= 0:
        return f"the deformation is not positive definite (leading minors {minors})"
    if not 0.5 <= minors[2] <= 2.0:
        return f"the volume would change by a factor {minors[2]:.4g}"
    return None


def gradient(reference: dict, y: list[float], forces: list[list[float]], stress: list[list[float]]) -> list[float]:
    """dE/dy at coordinates y from VASP's forces (eV/A) and stress (kbar), before symmetrization."""

    d = deformation(reference, y)
    gu = [[-x for x in row] for row in matmul(forces, d)]
    mean = [sum(row[k] for row in gu) / len(gu) for k in range(3)]
    gu = [[row[k] - mean[k] for k in range(3)] for row in gu]
    volume = abs(det3(reference["cell"])) * det3(d)
    g = matmul(stress, inv3(d))
    g = [[-volume * KBAR * 0.5 * (g[i][j] + g[j][i]) for j in range(3)] for i in range(3)]
    return [x for row in gu for x in row] + [x / reference["L"] for x in voigt(g)]


def solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """The solution of a small square linear system, by Gaussian elimination with partial pivoting."""

    n = len(rhs)
    rows = [list(row) + [b] for row, b in zip(matrix, rhs)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        for r in range(col + 1, n):
            factor = rows[r][col] / rows[col][col]
            rows[r] = [x - factor * p for x, p in zip(rows[r], rows[col])]
    x = [0.0] * n
    for r in reversed(range(n)):
        x[r] = (rows[r][n] - sum(rows[r][k] * x[k] for k in range(r + 1, n))) / rows[r][r]
    return x


def physical_size(reference: dict, y: list[float], g: list[float]) -> tuple[float, float]:
    """The largest atomic force (eV/A) and stress component (kbar) whose gradient at coordinates y is g.

    ``gradient`` inverted: the forces are -g_u (I + eps)^-1, and the stress S is the symmetric solution of
    sym(S (I + eps)^-1) = -G / (V KBAR), G the strain part of g as a tensor, a 6 x 6 linear system.
    """

    n3 = 3 * len(reference["frac"])
    d = deformation(reference, y)
    inverse = inv3(d)
    forces = [[-x for x in row] for row in matmul([g[i : i + 3] for i in range(0, n3, 3)], inverse)]
    volume = abs(det3(reference["cell"])) * det3(d)
    target = tensor([-x * reference["L"] / (volume * KBAR) for x in g[n3:]])
    pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    columns = []
    for i, j in pairs:
        unit = [[1.0 if {p, q} == {i, j} else 0.0 for q in range(3)] for p in range(3)]
        image = [
            [0.5 * (u + v) for u, v in zip(r1, r2)] for r1, r2 in zip(matmul(unit, inverse), matmul(inverse, unit))
        ]
        columns.append([image[p][q] for p, q in pairs])
    x = solve(transpose(columns), [target[p][q] for p, q in pairs])
    stress = [[0.0] * 3 for _ in range(3)]
    for (i, j), value in zip(pairs, x):
        stress[i][j] = stress[j][i] = value
    return maxima(forces, stress)


def maxima(forces: list[list[float]], stress: list[list[float]]) -> tuple[float, float]:
    """The largest atomic force (eV/A) and the largest stress component (kbar)."""

    return max(norm(row) for row in forces), max(abs(x) for row in stress for x in row)


# -- symmetry ---------------------------------------------------------------------------------------------------------
#
# An operation maps fractional coordinates s to s W + t with an integer W (entries -1, 0, 1, which covers the usual
# reduced cells; an operation outside that set is simply not used) that keeps the metric C C^T, and atom a onto atom
# perm[a] of the same species. In Cartesian row vectors it is x -> x Q + t C with Q = C^-1 W C, so it maps a
# displacement field u to u'[perm[a]] = u[a] Q and a strain eps to Q^T eps Q. Averaging a vector of the optimizer's
# coordinates over the operations projects it onto the part the symmetry leaves invariant.


def symmetry_operations(cell: list[list[float]], frac: list[list[float]], kinds: list[str], tolerance: float) -> list:
    """The symmetry operations of a structure as (Q, perm) pairs, within *tolerance* (A)."""

    metric = matmul(cell, transpose(cell))
    lengths = [math.sqrt(metric[i][i]) for i in range(3)]
    slack = 2 * tolerance * max(lengths) + tolerance**2
    candidates = list(itertools.product((-1, 0, 1), repeat=3))

    def product(u, v):
        return sum(u[i] * metric[i][j] * v[j] for i in range(3) for j in range(3))

    rows = [[w for w in candidates if abs(product(w, w) - metric[i][i]) <= slack] for i in range(3)]
    inverse = inv3(cell)
    n = len(frac)
    operations = []
    for w in itertools.product(*rows):
        if any(abs(product(w[i], w[j]) - metric[i][j]) > slack for i in range(3) for j in range(i + 1, 3)):
            continue
        if abs(det3([list(r) for r in w])) != 1:
            continue
        rotated = [[sum(s[k] * w[k][j] for k in range(3)) for j in range(3)] for s in frac]
        for target in range(n):
            if kinds[target] != kinds[0]:
                continue
            t = [frac[target][j] - rotated[0][j] for j in range(3)]
            perm = []
            for a in range(n):
                image = [rotated[a][j] + t[j] for j in range(3)]
                match = None
                for b in range(n):
                    if kinds[b] != kinds[a]:
                        continue
                    delta = [((image[j] - frac[b][j]) + 0.5) % 1.0 - 0.5 for j in range(3)]
                    if norm(matmul([delta], cell)[0]) <= tolerance:
                        match = b
                        break
                if match is None or match in perm:
                    break
                perm.append(match)
            if len(perm) == n:
                q = matmul(matmul(inverse, [list(map(float, r)) for r in w]), cell)
                operations.append({"Q": q, "perm": perm})
    return operations


def transform(operation: dict, v: list[float], n: int) -> list[float]:
    """A vector of the optimizer's coordinates (or gradient) mapped by one symmetry operation."""

    q, perm = operation["Q"], operation["perm"]
    out = [0.0] * len(v)
    for a in range(n):
        b = perm[a]
        for j in range(3):
            out[3 * b + j] = sum(v[3 * a + k] * q[k][j] for k in range(3))
    eps = tensor(v[3 * n :])
    out[3 * n :] = voigt(matmul(matmul(transpose(q), eps), q))
    return out


def symmetrize(reference: dict, v: list[float]) -> list[float]:
    """The part of v the starting structure's symmetry leaves invariant."""

    operations = reference["symmetry"]
    n = len(reference["frac"])
    total = [0.0] * len(v)
    for operation in operations:
        total = [t + x for t, x in zip(total, transform(operation, v, n))]
    return [t / len(operations) for t in total]


def atom_kinds(reference: dict) -> list[str]:
    return [s for s, c in zip(reference["species"], reference["counts"]) for _ in range(c)]


def initial_hessian(a: Attempt, reference: dict) -> list[list[float]] | None:
    """The optimizer's starting Hessian model: isotropic by default, or the job's hessian_seed."""

    n3 = 3 * len(reference["frac"])
    size = n3 + 6
    seed = parameter(a, "hessian_seed")
    if not seed:
        atom, bulk = float(parameter(a, "hessian_atom")), float(parameter(a, "hessian_bulk_modulus"))
        if not (math.isfinite(atom) and math.isfinite(bulk) and atom > 0 and bulk > 0):
            a.fail("vasp.input_invalid", f"hessian_atom {atom} and hessian_bulk_modulus {bulk} must be positive")
            return None
        strain = reference["V0"] * bulk / 160.21766208 / reference["L"] ** 2
        diagonal = [atom] * n3 + [strain] * 6
        return [[diagonal[i] if i == j else 0.0 for j in range(size)] for i in range(size)]
    path = a.payload / seed
    try:
        matrix = json.loads(path.read_text(encoding="utf-8"))["hessian"]
        matrix = [[float(x) for x in row] for row in matrix]
    except (OSError, ValueError, KeyError, TypeError) as exception:
        a.fail("vasp.hessian_seed_invalid", f"the hessian_seed {seed} could not be read: {exception}")
        return None
    if len(matrix) != size or any(len(row) != size for row in matrix):
        a.fail("vasp.hessian_seed_invalid", f"the hessian_seed {seed} is not {size} x {size} (3N + 6)")
        return None
    if not all(math.isfinite(x) for row in matrix for x in row):
        a.fail("vasp.hessian_seed_invalid", f"the hessian_seed {seed} has entries that are not finite")
        return None
    matrix = [[0.5 * (matrix[i][j] + matrix[j][i]) for j in range(size)] for i in range(size)]
    # Invariant under the structure's symmetry, or the optimizer's steps would break it. With
    # hessian_seed_symmetrize, a seed learned at lower symmetry is projected instead: the average of P H P^T over
    # the operations, positive definite whenever H is.
    scale = max(abs(x) for row in matrix for x in row)
    columns = [[1.0 if i == j else 0.0 for i in range(size)] for j in range(size)]
    rotations = [
        transpose([transform(operation, column, n3 // 3) for column in columns]) for operation in reference["symmetry"]
    ]
    images = [matmul(matmul(p, matrix), transpose(p)) for p in rotations]
    worst = max(abs(x - y) for image in images for row, row0 in zip(image, matrix) for x, y in zip(row, row0))
    if flag(parameter(a, "hessian_seed_symmetrize")):
        matrix = [[sum(image[i][j] for image in images) / len(images) for j in range(size)] for i in range(size)]
        a.log.append("note", f"hessian_seed {seed} projected onto the symmetry (largest change {worst:.3g})")
    elif worst > 1e-3 * scale:
        a.fail(
            "vasp.hessian_seed_invalid",
            f"the hessian_seed {seed} is not invariant under the structure's symmetry "
            f"(largest change {worst:.3g} of entries up to {scale:.3g})",
        )
        return None
    # Cholesky: the BFGS model must start positive definite.
    lower = [[0.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1):
            s = matrix[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if s <= 0:
                    a.fail("vasp.hessian_seed_invalid", f"the hessian_seed {seed} is not positive definite")
                    return None
                lower[i][i] = math.sqrt(s)
            else:
                lower[i][j] = s / lower[j][j]
    return matrix


# -- the driver's persistent state ------------------------------------------------------------------------------------


def driver_path(a: Attempt, count: int) -> Path:
    return a.workdir / DRIVER / f"driver-{count:03d}.json"


def save_driver(a: Attempt, driver: dict, count: int) -> None:
    path = driver_path(a, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(driver, indent=1), encoding="utf-8")
    temporary.replace(path)


def optimizer(a: Attempt, driver: dict) -> TrustBFGS:
    reference = driver["reference"]
    trust_atom, trust_strain = float(parameter(a, "trust_atom")), float(parameter(a, "trust_strain"))
    n = len(reference["frac"])

    def step_scale(p: list[float], trust_factor: float) -> float:
        # The bounds hold in the current cell: an atom's physical displacement is du (I + eps), and the strain
        # relative to the current cell is (I + eps)^-1 d eps.
        here = deformation(reference, opt.current[0])
        displacement = max(norm(matmul([p[3 * i : 3 * i + 3]], here)[0]) for i in range(n))
        strain = max(abs(x) for row in matmul(inv3(here), strain_of(reference, p)) for x in row)
        limits = [math.inf]
        if displacement > 0:
            limits.append(trust_atom * trust_factor / displacement)
        if strain > 0:
            limits.append(trust_strain * trust_factor / strain)
        return min(limits)

    state = driver["optimizer"]
    opt = TrustBFGS(state["B"], step_scale, float(parameter(a, "rise_tolerance")))
    opt.trust_factor = state["trust_factor"]
    opt.current = None if state["current"] is None else tuple(state["current"])
    opt.history = state["history"]
    return opt


def store_optimizer(driver: dict, opt: TrustBFGS) -> None:
    driver["optimizer"] = {
        "B": opt.B,
        "trust_factor": opt.trust_factor,
        "current": None if opt.current is None else list(opt.current),
        "history": opt.history,
    }


def slurm_duration(seconds: float) -> str:
    """A whole number of seconds as a Slurm M:S duration (any number of minutes)."""

    whole = max(1, math.ceil(seconds))
    return f"{whole // 60}:{whole % 60:02d}"


def cost_key(a: Attempt, plan: dict) -> str:
    """Which evaluations an evaluation's cost is estimated from: coarse:loose, coarse:tight or final."""

    if plan["kind"] == "final":
        return "final"
    return "coarse:loose" if plan["ediff"] > float(parameter(a, "ediff_tight")) * (1 + 1e-9) else "coarse:tight"


def expected_seconds(a: Attempt, driver: dict, key: str) -> float | None:
    """The expected VASP seconds of the next evaluation of this cost key, or None without any estimate."""

    slowest = driver["slowest"]
    if slowest.get(key) is not None:
        return slowest[key]
    first = float(parameter(a, "expected_seconds")) or None
    loose = slowest.get("coarse:loose") or first
    if key == "coarse:loose":
        return loose
    tight = None if loose is None else float(parameter(a, "tight_cost_factor")) * loose
    if key == "coarse:tight":
        return tight
    coarse = slowest.get("coarse:tight") or tight
    return None if coarse is None else float(parameter(a, "final_cost_factor")) * coarse


def time_needed(a: Attempt, driver: dict, plan: dict, *, staged: bool = False) -> float | None:
    """The seconds an attempt must have left to start this evaluation: before staging, or (*staged*) after it."""

    expected = expected_seconds(a, driver, cost_key(a, plan))
    if expected is None:
        return None
    overhead = driver["overhead"]
    return (
        float(parameter(a, "time_factor")) * expected + overhead["archive"] + (0.0 if staged else overhead["staging"])
    )


def maxtime_of(a: Attempt) -> int | None:
    maxtime = a.context.resources.get("maxtime")
    return None if maxtime is None else int(maxtime)


def with_mintime(a: Attempt, seconds: float | None) -> dict | None:
    """The attempt's resource requirement with mintime set to *seconds*, which the caller checked fits maxtime."""

    if seconds is None:
        return None
    resources: dict[str, object] = dict(a.context.resources)
    maxtime = maxtime_of(a)
    if maxtime is not None:
        resources["maxtime"] = slurm_duration(maxtime)
    resources["mintime"] = slurm_duration(seconds)
    return resources


def insufficient_time(a: Attempt, needed: float, plan: dict) -> bool:
    """Fail the job when even a fresh attempt could not fit the evaluation into its maxtime; True if it did."""

    maxtime = maxtime_of(a)
    if maxtime is None or needed <= maxtime:
        return False
    a.fail(
        "vasp.insufficient_time",
        f"the next {plan['kind']} evaluation needs {needed:.0f} s, more than the attempt maxtime of {maxtime} s; "
        "raise the job's maxtime and continue it",
        details={
            "needed_seconds": needed,
            "maxtime_seconds": maxtime,
            "kind": plan["kind"],
            "purpose": plan["purpose"],
        },
    )
    return True


def guard(a: Attempt, needed: float | None, plan: dict, when: str) -> bool:
    """True when the attempt has the time for this evaluation; else defer (or fail) it and return False."""

    if needed is None:
        if a.deadline is not None and when == "before staging":
            a.log.append(
                "warning", f"time guard: no estimate for the first {plan['kind']} evaluation (expected_seconds)"
            )
        return True
    if insufficient_time(a, needed, plan):
        return False
    if a.deadline is not None and a.deadline - time.time() < needed:
        a.log.append(
            "note",
            f"time guard ({when}): {a.deadline - time.time():.0f} s left, the next {plan['kind']} evaluation needs "
            f"{needed:.0f} s; deferred to a manager with enough time",
        )
        a.advance("evaluate", resources=with_mintime(a, needed))
        return False
    return True


# -- steps ------------------------------------------------------------------------------------------------------------


@run.step
def prepare(a: Attempt) -> None:
    """Stage the structure, POTCAR and seed WAVECARs, set up the optimizer, and go on to the first evaluation."""

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
    normalize_poscar_handedness(a.workdir / "POSCAR")
    try:
        reference = read_poscar(a.workdir / "POSCAR")
    except (ValueError, IndexError) as exception:
        a.fail("vasp.input_invalid", f"the starting structure cannot be relaxed by the driver: {exception}")
        return
    reference["V0"] = abs(det3(reference["cell"]))
    reference["L"] = reference["V0"] ** (1 / 3)
    reference["symmetry"] = symmetry_operations(
        reference["cell"], reference["frac"], atom_kinds(reference), float(parameter(a, "symmetry_tolerance"))
    )
    hessian = initial_hessian(a, reference)
    if hessian is None:
        return
    wavecars = {}
    for kind in ("coarse", "final"):
        seed = a.payload / str(parameter(a, f"{kind}_wavecar"))
        if seed.is_file():
            (a.workdir / DRIVER / kind).mkdir(parents=True, exist_ok=True)
            shutil.copyfile(seed, a.workdir / DRIVER / kind / "WAVECAR-000")
            # A seed's inputs are not known: it is staged as it is, and VASP's output says whether it was read.
            wavecars[kind] = {"file": "WAVECAR-000", "evaluation": "seed", "signature": None}
    driver = {
        "reference": reference,
        "optimizer": {"B": hessian, "trust_factor": 1.0, "current": None, "history": []},
        "next": {"y": [0.0] * (3 * len(reference["frac"]) + 6), "kind": "coarse", "purpose": "start"},
        "level": "loose",
        "current": None,
        "offset": None,
        "offsets": [],
        "finals": [],
        "corrections": 0,
        "coarse_evaluations": 0,
        "meshes": {},
        "slowest": {},
        "overhead": {"staging": 0.0, "archive": 0.0},
        "wavecars": wavecars,
        "switched_at": None,
        "interruptions": [],
    }
    driver["next"]["ediff"] = float(parameter(a, "ediff_loose"))
    save_driver(a, driver, 0)
    a.advance("evaluate", state={"evaluations": 0, "stage": "coarse", "level": "loose"})


class UnsupportedInputs(Exception):
    """The inputs of an evaluation would differ from what the driver's surface was learned on."""


def kpoints_lines(directory: Path) -> list[str]:
    """The generation mode, centering and mesh lines of an automatic KPOINTS."""

    return [line.strip() for line in (directory / "KPOINTS").read_text(encoding="utf-8").splitlines()[1:4]]


def check_remedy_history(path: Path) -> None:
    """Refuse a job remedy history that would replay anything but an electronic remedy.

    :raises UnsupportedInputs: When the history cannot be read or holds another operation.
    """

    if not path.is_file():
        return
    try:
        events = json.loads(path.read_text(encoding="utf-8")).get("events", [])
        operations = [str(item["operation"]) for event in events for item in event.get("changes", [])]
    except (OSError, ValueError, AttributeError, KeyError, TypeError) as exception:
        raise UnsupportedInputs(f"the job's remedy history {path.name} cannot be read: {exception}") from exception
    surface = sorted({operation for operation in operations if operation not in ELECTRONIC_REMEDIES})
    if surface:
        raise UnsupportedInputs(
            f"the job's remedy history would replay {', '.join(surface)}, which changes the energy surface or the "
            "structure the driver relaxes"
        )


def stage_inputs(a: Attempt, driver: dict, plan: dict) -> dict:
    """Write the INCAR, KPOINTS, POSCAR and WAVECAR of one evaluation, and return what was staged.

    :raises UnsupportedInputs: Before anything is written, when the job's remedy history holds a remedy other than an
        electronic one; after, when a remedy in force changed the k-point mesh all the same.
    """

    check_remedy_history(job_remedy_history_path(a.payload))
    kind = plan["kind"]
    reference = driver["reference"]
    cell, frac = geometry(reference, plan["y"])
    write_poscar(a.workdir / "POSCAR", reference, cell, frac)
    incar = a.workdir / "INCAR"
    override = a.payload / "files" / f"INCAR.{kind}"
    shutil.copyfile(override if override.is_file() else TEMPLATES / f"INCAR.{kind}", incar)
    density = parameter(a, "coarse_kpoint_density" if kind == "coarse" else "kpoint_density")
    centering = str(parameter(a, "centering"))
    library = None
    if not (a.workdir / "POTCAR").is_file():
        library = a.setting("vasp.pseudo_library", a.parameter("pseudopotential_library", None)) or None
    parallel_value = a.parameter("parallel_value", None)
    options = VaspPreparationOptions(
        kpoint_density=float(density),
        centering=centering,
        accuracy_per_atom=a.parameter("accuracy_per_atom", 0.001),
        pseudopotential_library=library,
        parallel_tag=a.parameter("parallel_tag", None) or None,
        parallel_value=None if parallel_value is None else int(parallel_value),
        normalize_handedness=False,
        incar_tags={
            **dict(a.parameter("incar_tags", {}) or {}),
            **dict(parameter(a, f"{kind}_incar_tags") or {}),
            "EDIFF": plan["ediff"],
        },
    )
    prepare_vasp_inputs(options, directory=a.workdir)
    # One base mesh per kind for the whole relaxation: the job's, or the one derived for the first evaluation of the
    # kind, so a cell that changes during the relaxation never changes its mesh.
    mesh = str(parameter(a, f"{kind}_mesh") or "").split() or driver["meshes"].get(kind)
    if mesh:
        write_automatic_kpoints([int(x) for x in mesh], a.workdir / "KPOINTS", centering=centering)
    base = kpoints_lines(a.workdir)
    applied = reapply_vasp_remedies(directory=a.workdir, history_path=job_remedy_history_path(a.payload))
    effective = kpoints_lines(a.workdir)
    if effective != base:
        raise UnsupportedInputs(f"a remedy in force changed the {kind} k-points from {base} to {effective}")
    tags = read_incar(incar)
    signature = {"kind": kind, "kpoints": effective, **{tag: tags.get(tag) for tag in WAVECAR_TAGS}}
    # The WAVECAR of the last good evaluation of this kind, if it was written for these inputs (copied, so a failed
    # run never spoils it).
    cache = driver["wavecars"].get(kind)
    stored = None if cache is None else a.workdir / DRIVER / kind / cache["file"]
    usable = stored is not None and stored.is_file() and cache["signature"] in (None, signature)
    wavecar = a.workdir / "WAVECAR"
    if usable:
        shutil.copyfile(stored, wavecar)
    elif wavecar.exists():
        wavecar.unlink()
    # What makes this one evaluation of the driver wins over everything else, except a remedy's charge reset.
    reset = [value for operation, value in applied if operation == "incar.ICHARG"]
    update_incar(
        {
            "NSW": 0,
            "IBRION": -1,
            "ISIF": 2,
            "LWAVE": ".TRUE.",
            "ISTART": 1 if usable else 0,
            "ICHARG": reset[-1] if reset else (0 if usable else 2),
        },
        incar,
    )
    return {
        "restart_requested": usable,
        "restart_source": cache["evaluation"] if usable else None,
        "base_mesh": base[2].split(),
        "signature": signature,
        "remedies": [[operation, value] for operation, value in applied],
    }


def read_evaluation(directory: Path) -> dict:
    """Energy, forces, stress, cell, positions and SCF steps of the last calculation in vasprun.xml."""

    root = ElementTree.parse(directory / "vasprun.xml").getroot()
    calculation = root.findall("calculation")[-1]

    def table(node, name: str) -> list[list[float]]:
        rows = node.findall(f"varray[@name='{name}']/v")
        if not rows:
            raise ValueError(f"vasprun.xml has no {name}")
        return [[float(x) for x in row.text.split()] for row in rows]

    energy = calculation.find("energy/i[@name='e_fr_energy']")
    if energy is None:
        raise ValueError("vasprun.xml has no free energy")
    return {
        "energy": float(energy.text),
        "forces": table(calculation, "forces"),
        "stress": table(calculation, "stress"),
        "cell": table(calculation.find("structure/crystal"), "basis"),
        "positions": table(calculation.find("structure"), "positions"),
        "scf_steps": len(calculation.findall("scstep")),
    }


def restart_observed(directory: Path, requested: bool) -> bool | str:
    """Whether VASP's stdout shows that it read the staged WAVECAR: True, False, or "unknown"."""

    if not requested:
        return False
    stdout = directory / "vasp.out"
    text = stdout.read_text(encoding="utf-8", errors="replace") if stdout.is_file() else ""
    if WAVECAR_READ in text:
        return True
    if WAVECAR_NOT_READ in text:
        return False
    return "unknown"


def same_geometry(cell: list[list[float]], frac: list[list[float]], result: dict) -> bool:
    cells = max(abs(a - b) for ra, rb in zip(cell, result["cell"]) for a, b in zip(ra, rb))
    sites = max(abs((a - b + 0.5) % 1.0 - 0.5) for ra, rb in zip(frac, result["positions"]) for a, b in zip(ra, rb))
    return len(frac) == len(result["positions"]) and cells < 1e-5 and sites < 1e-5


def load_driver(a: Attempt) -> tuple[int, dict] | None:
    """The newest committed driver snapshot: the job state's, or a later one an interrupted attempt wrote."""

    count = int(a.state.get("evaluations", 0))
    try:
        driver = json.loads(driver_path(a, count).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exception:
        a.fail("vasp.driver_state_missing", f"the driver state after evaluation {count} cannot be read: {exception}")
        return None
    # A snapshot is written whole (renamed into place) after its evaluation's record, so one beyond the job state
    # is a finished evaluation whose attempt stopped before the state was merged.
    while driver_path(a, count + 1).is_file():
        count += 1
        driver = json.loads(driver_path(a, count).read_text(encoding="utf-8"))
        a.log.append("note", f"evaluation {count} was finished by an interrupted attempt; continuing after it")
        a.state.merge(step_state(driver, count))
    return count, driver


def step_state(driver: dict, count: int) -> dict:
    stage = driver["next"]["kind"]
    return {"evaluations": count, "stage": "done" if stage == "publish" else stage, "level": driver["level"]}


def drop_stale_wavecars(a: Attempt, driver: dict) -> None:
    """Delete every cached WAVECAR the snapshot does not name (older ones, or one of an unfinished evaluation)."""

    for kind in ("coarse", "final"):
        keep = (driver["wavecars"].get(kind) or {}).get("file")
        directory = a.workdir / DRIVER / kind
        if directory.is_dir():
            for path in directory.glob("WAVECAR*"):
                if path.name != keep:
                    path.unlink()


@run.step
def evaluate(a: Attempt) -> None:
    """Run one single point at the optimizer's geometry, tell the optimizer, and decide what comes next."""

    loaded = load_driver(a)
    if loaded is None:
        return
    count, driver = loaded
    drop_stale_wavecars(a, driver)
    plan = driver["next"]
    if plan["kind"] == "publish":
        # The last final static is done: only the publication is left, from its archive.
        finish(a, driver, step_state(driver, count))
        return
    kind = plan["kind"]
    if not guard(a, time_needed(a, driver, plan), plan, "before staging"):
        return
    argv = shlex.split(a.setting("vasp.command", a.parameter("vasp_command", None)) or "")
    if not argv:
        a.fail(
            "vasp.command_missing",
            "no VASP command is configured: set it with `httk workspace settings set --key vasp.command --value '...' WORKSPACE`, "
            "or set HTTK_VASP_COMMAND on the machine that runs this job, or give the job a vasp_command parameter",
        )
        return
    started = time.monotonic()
    # A rerun must not read the previous run's outputs; the WAVECAR is staged afresh below.
    clean_vasp_outputs(a.workdir, keep=("WAVECAR",))
    try:
        staged = stage_inputs(a, driver, plan)
    except UnsupportedInputs as exception:
        a.fail("vasp.unsupported_remedy", f"evaluation {count + 1} cannot run on the driver's surface: {exception}")
        return
    staging = time.monotonic() - started
    if not guard(a, time_needed(a, driver, plan, staged=True), plan, "after staging"):
        return
    started = time.monotonic()
    try:
        report = run_vasp(argv, directory=a.workdir, timeout=a.parameter("timeout", 86400.0))
    except OSError as exception:
        a.fail("vasp.failed", f"could not run VASP at all: {exception}")
        return
    seconds = time.monotonic() - started
    finished = time.monotonic()
    a.log.append("note", f"VASP {report.classification} in evaluation {count + 1} ({kind}, {plan['purpose']})")
    if report.classification != "completed":
        stop = stopped_by(a, report)
        handling = {"deadline": "deferred", "signal": "vasp.interrupted", "timeout": "vasp.evaluation_timeout"}
        keep_incomplete(a, driver, count + 1, report, seconds, handling.get(stop, "remedy ladder"))
        termination = report.process.termination
        if stop == "deadline":
            interrupted(a, driver, count, seconds, termination)
        elif stop == "signal":
            left = None if a.deadline is None else a.deadline - time.time()
            a.fail(
                "vasp.interrupted",
                f"VASP was stopped by a signal ({termination}) in evaluation {count + 1}, not at the attempt deadline "
                f"({'none' if left is None else f'{left:.0f} s left'}): a cancellation or a stop from outside httk; "
                "retry the job to run the evaluation again",
                details={"termination": termination, "seconds": seconds, "seconds_to_deadline": left},
            )
        elif stop == "timeout":
            a.fail(
                "vasp.evaluation_timeout",
                f"VASP ran into the runner's timeout ({termination}, {seconds:.0f} s) in evaluation {count + 1}",
                details={"termination": termination, "seconds": seconds, "timeout": a.parameter("timeout", 86400.0)},
            )
        else:
            remedy(a, report)
        return
    try:
        result = read_evaluation(a.workdir)
    except (OSError, ValueError, IndexError, AttributeError, ElementTree.ParseError) as exception:
        a.fail("vasp.no_forces", f"evaluation {count + 1} left no readable forces and stress: {exception}")
        return
    cell, frac = geometry(driver["reference"], plan["y"])
    if not same_geometry(cell, frac, result):
        a.fail("vasp.geometry_mismatch", f"VASP ran evaluation {count + 1} at another geometry than the driver wrote")
        return
    index = count + 1
    # The good WAVECAR becomes the restart of the next evaluation of this kind, under a name of its own that only a
    # snapshot written after this evaluation points to.
    (a.workdir / DRIVER / kind).mkdir(parents=True, exist_ok=True)
    if (a.workdir / "WAVECAR").is_file():
        name = f"WAVECAR-{index:03d}"
        (a.workdir / "WAVECAR").replace(a.workdir / DRIVER / kind / name)
        driver["wavecars"][kind] = {"file": name, "evaluation": index, "signature": staged["signature"]}
    driver["meshes"].setdefault(kind, staged["base_mesh"])
    key = cost_key(a, plan)
    driver["slowest"][key] = max(seconds, driver["slowest"].get(key) or 0.0)
    ran = read_incar(a.workdir / "INCAR")
    fmax, smax = maxima(result["forces"], result["stress"])
    record = {
        "evaluation": index,
        "kind": kind,
        "purpose": plan["purpose"],
        "level": driver["level"] if kind == "coarse" else "final",
        "mesh": " ".join(staged["base_mesh"]),
        "EDIFF": float(ran.get("EDIFF", plan["ediff"])),
        "PRECFOCK": ran.get("PRECFOCK"),
        "ICHARG": ran.get("ICHARG"),
        "remedies": staged["remedies"],
        "energy": result["energy"],
        "max_force": fmax,
        "max_stress": smax,
        "forces": result["forces"],
        "stress": result["stress"],
        "cell": cell,
        "frac": frac,
        "y": plan["y"],
        "scf_steps": result["scf_steps"],
        "seconds": seconds,
        "staging_seconds": staging,
        "restart_requested": staged["restart_requested"],
        "restart_source": staged["restart_source"],
        "restart_observed": restart_observed(a.workdir, staged["restart_requested"]),
        "offset_applied": driver["offset"] is not None and kind == "coarse",
    }
    archive = a.workdir / DRIVER / "evaluations" / f"{index:03d}"
    archive.mkdir(parents=True, exist_ok=True)
    for name in EVALUATION_ARCHIVE:
        if (a.workdir / name).is_file():
            shutil.copyfile(a.workdir / name, archive / name)
    outcome = decide(a, driver, record)
    (archive / "record.json").write_text(json.dumps(record, indent=1), encoding="utf-8")
    if outcome is None:
        return
    overhead = driver["overhead"]
    overhead["staging"] = max(overhead["staging"], staging)
    overhead["archive"] = max(overhead["archive"], time.monotonic() - finished)
    if outcome == "publish":
        driver["next"] = {"kind": "publish", "purpose": "publish"}
    save_driver(a, driver, index)
    state = step_state(driver, index)
    state["last_energy"] = result["energy"]
    if outcome == "publish":
        finish(a, driver, state)
        return
    needed = time_needed(a, driver, driver["next"])
    maxtime = maxtime_of(a)
    if needed is not None and maxtime is not None and needed > maxtime:
        # This evaluation is committed; the job stops before the next one, which no attempt could fit.
        a.state.merge(state)
        insufficient_time(a, needed, driver["next"])
        return
    a.advance("evaluate", state=state, resources=with_mintime(a, needed))


def stopped_by(a: Attempt, report) -> str | None:
    """Who stopped a VASP run that did not complete: "deadline", "signal" or "timeout"; None when it ended by itself.

    Only the process's termination tells: a signal to the runner (the manager's SIGTERM at maxtime, drain or
    cancellation) is "deadline" when it came at the attempt's deadline and "signal" otherwise; the runner's own timeout
    is "timeout". A run that exited, with any status, or that a diagnostic stopped, ended by itself however close to
    the deadline, and goes to the remedy ladder.
    """

    termination = str(report.process.termination)
    if termination.startswith("timeout"):
        return "timeout"
    if termination.startswith("signal_"):
        if a.deadline is not None and time.time() >= a.deadline - DEADLINE_SLACK:
            return "deadline"
        return "signal"
    return None


def keep_incomplete(a: Attempt, driver: dict, index: int, report, seconds: float, handling: str) -> None:
    """Keep the files of a run of evaluation *index* that did not complete, before any later run replaces them."""

    directory = a.workdir / DRIVER / "incomplete"
    number = 1 + sum(1 for path in directory.glob(f"{index:03d}-*")) if directory.is_dir() else 1
    archive = directory / f"{index:03d}-{number:02d}"
    archive.mkdir(parents=True, exist_ok=True)
    for name in EVALUATION_ARCHIVE:
        if (a.workdir / name).is_file():
            shutil.copyfile(a.workdir / name, archive / name)
    plan = driver["next"]
    record = {
        "evaluation": index,
        "run": number,
        "kind": plan["kind"],
        "purpose": plan["purpose"],
        "classification": report.classification,
        "termination": report.process.termination,
        "returncode": report.process.returncode,
        "seconds": seconds,
        "seconds_to_deadline": None if a.deadline is None else a.deadline - time.time(),
        "handling": handling,
        "diagnostics": [
            {"code": item.code, "severity": str(item.severity), "summary": item.summary} for item in report.diagnostics
        ],
    }
    (archive / "record.json").write_text(json.dumps(record, indent=1), encoding="utf-8")


def interrupted(a: Attempt, driver: dict, count: int, seconds: float, termination: str) -> None:
    """Defer an evaluation the manager stopped at the attempt's deadline: no remedy, the same evaluation again later.

    The stopped run took at least *seconds*, which becomes its kind's estimate in the committed snapshot (only the
    estimate and the list of interruptions change), so the next attempt must have the time for it. More than
    maximum_interruptions of them fail the job: every deferral is a new activation, which no other cap bounds.
    """

    plan = driver["next"]
    key = cost_key(a, plan)
    driver["slowest"][key] = max(seconds, driver["slowest"].get(key) or 0.0)
    interruptions = driver.setdefault("interruptions", [])
    interruptions.append({"evaluation": count + 1, "seconds": seconds, "termination": termination})
    save_driver(a, driver, count)
    maximum = int(parameter(a, "maximum_interruptions"))
    if len(interruptions) > maximum:
        a.fail(
            "vasp.interrupted_repeatedly",
            f"VASP was stopped at the attempt deadline {len(interruptions)} times (maximum_interruptions {maximum}), "
            f"last in evaluation {count + 1} after {seconds:.0f} s",
            details={"interruptions": interruptions, "maximum_interruptions": maximum},
        )
        return
    a.log.append(
        "note",
        f"VASP was stopped at the attempt deadline ({termination}) after {seconds:.0f} s; evaluation {count + 1} waits "
        f"(interruption {len(interruptions)} of at most {maximum})",
    )
    needed = time_needed(a, driver, plan)
    if not insufficient_time(a, needed, plan):
        a.advance("evaluate", resources=with_mintime(a, needed))


def remedy(a: Attempt, report) -> None:
    """Apply one electronic remedy of the reviewed ladder and retry the same evaluation, or fail with the diagnosis."""

    applied = int(a.state.get("remedies", 0))
    maximum = int(a.parameter("maximum_remedies", 8))
    history = job_remedy_history_path(a.payload)
    # Planned before the budget is consulted, so a job that stops here still says which remedy it would have applied.
    try:
        decision = plan_vasp_remedy(
            report.diagnostics,
            directory=a.workdir,
            history_path=history,
            policy=a.parameter("remedy_policy", "reviewed-v1"),
        )
    except ValueError as exception:
        a.fail("vasp.failed", f"planning a VASP remedy failed: {exception}")
        return
    if decision.give_up or applied >= maximum:
        if applied >= maximum:
            message = f"VASP {report.classification} after {applied} remedies"
        else:
            message = f"VASP {report.classification} with no remaining remedy"
        a.fail("vasp.failed", message, details=decision.as_mapping())
        return
    surface = [operation for operation, _ in decision.changes if operation not in ELECTRONIC_REMEDIES]
    if surface:
        a.fail(
            "vasp.unsupported_remedy",
            f"VASP {report.classification}: the {decision.policy} remedy for {decision.problem} changes "
            f"{', '.join(surface)}, which would change the energy surface or the structure the driver relaxes",
            details=decision.as_mapping(),
        )
        return
    apply_vasp_remedy(decision, directory=a.workdir, history_path=history)
    a.state.merge({"remedies": applied + 1, "classification": report.classification})
    a.log.append("note", f"applied a remedy for {decision.problem}")
    # The retry prepares the evaluation again, re-applying every remedy, at the driver's geometry.
    a.retry(f"applied the {decision.policy} remedy for {decision.problem}")


def corrected(driver: dict, record: dict) -> tuple[list[list[float]], list[list[float]]]:
    """The coarse forces and stress of a record with the current offset added: the prediction of the final mesh."""

    offset = driver["offset"]
    if offset is None or record["kind"] != "coarse":
        return record["forces"], record["stress"]
    forces = [[f + d for f, d in zip(row, drow)] for row, drow in zip(record["forces"], offset["forces"])]
    stress = [[s + d for s, d in zip(row, drow)] for row, drow in zip(record["stress"], offset["stress"])]
    return forces, stress


def model(driver: dict, record: dict) -> tuple[float, list[float]]:
    """The energy and gradient the optimizer sees for a coarse record: on the corrected surface once there is one."""

    reference = driver["reference"]
    g = gradient(reference, record["y"], record["forces"], record["stress"])
    energy = record["energy"]
    offset = driver["offset"]
    if offset is not None:
        g = [x + d for x, d in zip(g, offset["gradient"])]
        energy += sum(d * (y - y0) for d, y, y0 in zip(offset["gradient"], record["y"], offset["y"]))
        # How far the physical prediction (offset forces and stress carried with the cell) has moved from the model,
        # both in the symmetric subspace the optimizer works in (so it is zero at the anchor), as forces and stress.
        moved = symmetrize(reference, gradient(reference, record["y"], offset["forces"], offset["stress"]))
        force, stress = physical_size(reference, record["y"], [m - d for m, d in zip(moved, offset["gradient"])])
        record["model_energy"] = energy
        record["correction_model_discrepancy"] = {"force": force, "stress": stress}
    return energy, symmetrize(reference, g)


def tell(a: Attempt, driver: dict, record: dict) -> bool:
    """Give a coarse evaluation to the optimizer and keep the accepted point's checks; False after a failure."""

    current = driver["current"]
    if current is not None and not math.isclose(record["EDIFF"], current["EDIFF"], rel_tol=1e-9):
        a.fail(
            "vasp.driver_inconsistent",
            f"evaluation {record['evaluation']} ran at EDIFF {record['EDIFF']:g}, the accepted point at "
            f"{current['EDIFF']:g}; the driver does not compare evaluations of two EDIFF levels",
        )
        return False
    energy, g = model(driver, record)
    opt = optimizer(a, driver)
    accepted = opt.tell(record["y"], energy, g)
    store_optimizer(driver, opt)
    accept(driver, record, accepted, opt.trust_factor)
    return True


def accept(driver: dict, record: dict, accepted: bool, trust_factor: float) -> None:
    forces, stress = corrected(driver, record)
    fmax, smax = maxima(forces, stress)
    record.update(accepted=accepted, trust_factor=trust_factor, corrected_max_force=fmax, corrected_max_stress=smax)
    if accepted:
        driver["current"] = {
            "evaluation": record["evaluation"],
            "max_force": fmax,
            "max_stress": smax,
            "EDIFF": record["EDIFF"],
            "y": record["y"],
        }


def rebase(a: Attempt, driver: dict, record: dict) -> bool:
    """Make the tight-EDIFF evaluation of the accepted point the optimizer's accepted point; False after a failure."""

    current = driver["current"]
    tight = float(parameter(a, "ediff_tight"))
    if record["y"] != current["y"] or record["EDIFF"] > tight * (1 + 1e-9) or driver["offset"] is not None:
        a.fail(
            "vasp.driver_inconsistent",
            f"evaluation {record['evaluation']} is not a tight-EDIFF evaluation of the accepted point "
            f"{current['evaluation']}",
        )
        return False
    loose = evaluation_record(a, current["evaluation"])
    energy, g = model(driver, record)
    opt = optimizer(a, driver)
    # No acceptance test and no curvature update: the two energies and gradients belong to different EDIFF levels.
    opt.current = (list(record["y"]), energy, g)
    opt.history.append(
        {"accepted": True, "energy": energy, "trust_factor": opt.trust_factor, "damped": False, "rebased": True}
    )
    store_optimizer(driver, opt)
    record["rebased"] = True
    record["ediff_shift"] = {
        "from_evaluation": loose["evaluation"],
        "energy": record["energy"] - loose["energy"],
        "max_force": max(norm([f - c for f, c in zip(r, rc)]) for r, rc in zip(record["forces"], loose["forces"])),
        "max_stress": max(abs(s - c) for r, rc in zip(record["stress"], loose["stress"]) for s, c in zip(r, rc)),
    }
    accept(driver, record, True, opt.trust_factor)
    return True


def within(a: Attempt, fmax: float, smax: float) -> bool:
    return fmax <= float(parameter(a, "force_limit")) and smax <= float(parameter(a, "stress_limit"))


def plan_next(a: Attempt, driver: dict, y: list[float], kind: str, purpose: str, ediff: float) -> str:
    driver["next"] = {"y": list(y), "kind": kind, "purpose": purpose, "ediff": ediff}
    return "evaluate"


def exhausted(a: Attempt, driver: dict, purpose: str) -> bool:
    """Fail the job when the coarse evaluations are used up before the next one (*purpose*); True if it did."""

    if driver["coarse_evaluations"] < int(parameter(a, "maximum_evaluations")):
        return False
    current = driver["current"]
    a.fail(
        "vasp.relax_not_converged",
        f"the coarse relaxation did not converge in {driver['coarse_evaluations']} evaluations"
        + (" (the budget ran out before the tight repeat of the switch point)" if purpose == "tighten" else ""),
        details={
            "evaluations": driver["coarse_evaluations"],
            "max_force": current["max_force"],
            "max_stress": current["max_stress"],
            "level": driver["level"],
            "next": purpose,
        },
    )
    return True


def next_coarse(a: Attempt, driver: dict) -> str | None:
    """Go to the final static when the accepted point is converged at the tight EDIFF, else ask for a step."""

    tight = float(parameter(a, "ediff_tight"))
    current = driver["current"]
    converged = within(a, current["max_force"], current["max_stress"])
    if driver["level"] != "loose" and current["EDIFF"] <= tight * (1 + 1e-9) and converged:
        return plan_next(a, driver, current["y"], "final", "final", float(parameter(a, "final_ediff")))
    if exhausted(a, driver, "step"):
        return None
    reference = driver["reference"]
    opt = optimizer(a, driver)
    y0 = opt.current[0]
    # The step is kept in the symmetric subspace (it is already, up to rounding, with a symmetric model).
    step = symmetrize(reference, [y - x for y, x in zip(opt.ask(), y0)])
    y = [x + s for x, s in zip(y0, step)]
    problem = cell_problem(reference, y)
    if problem is not None:
        a.fail("vasp.cell_invalid", f"the optimizer proposed a cell that cannot be written: {problem}")
        return None
    ediff = float(parameter(a, "ediff_loose" if driver["level"] == "loose" else "ediff_tight"))
    return plan_next(a, driver, y, "coarse", "step", ediff)


def evaluation_record(a: Attempt, index: int) -> dict:
    path = a.workdir / DRIVER / "evaluations" / f"{index:03d}" / "record.json"
    return json.loads(path.read_text(encoding="utf-8"))


def measure_offset(a: Attempt, driver: dict, record: dict) -> dict | None:
    """The final-mesh minus the coarse-mesh forces and stress at one geometry, or None after a failure."""

    coarse = evaluation_record(a, driver["current"]["evaluation"])
    tight = float(parameter(a, "ediff_tight"))
    if coarse["kind"] != "coarse" or coarse["y"] != record["y"] or coarse["EDIFF"] > tight * (1 + 1e-9):
        a.fail(
            "vasp.offset_pairing",
            f"final evaluation {record['evaluation']} has no tight-EDIFF coarse partner at its geometry "
            f"(accepted point {coarse['evaluation']}, EDIFF {coarse['EDIFF']:g})",
        )
        return None
    forces = [[f - c for f, c in zip(row, crow)] for row, crow in zip(record["forces"], coarse["forces"])]
    stress = [[f - c for f, c in zip(row, crow)] for row, crow in zip(record["stress"], coarse["stress"])]
    fmax, smax = maxima(forces, stress)
    offset = {
        "coarse_evaluation": coarse["evaluation"],
        "final_evaluation": record["evaluation"],
        "purpose": record["purpose"],
        "coarse_EDIFF": coarse["EDIFF"],
        "final_EDIFF": record["EDIFF"],
        "energy": record["energy"] - coarse["energy"],
        "max_force": fmax,
        "max_stress": smax,
        "forces": forces,
        "stress": stress,
        "y": record["y"],
    }
    driver["offsets"].append({k: v for k, v in offset.items() if k != "y"})
    return offset


def decide(a: Attempt, driver: dict, record: dict) -> str | None:
    """Update the driver with one evaluation and plan the next: "evaluate", "publish", or None after a failure."""

    if record["kind"] == "coarse":
        driver["coarse_evaluations"] += 1
        if record["purpose"] == "tighten":
            if not rebase(a, driver, record):
                return None
            if flag(parameter(a, "offset_probe")):
                return plan_next(a, driver, record["y"], "final", "probe", float(parameter(a, "final_ediff")))
            return next_coarse(a, driver)
        if not tell(a, driver, record):
            return None
        current = driver["current"]
        resolvable = current["max_force"] <= float(parameter(a, "switch_force")) and current["max_stress"] <= float(
            parameter(a, "switch_stress")
        )
        if driver["level"] == "loose" and resolvable:
            # The tight repeat is mandatory: without budget for it the job stops instead of skipping it.
            if exhausted(a, driver, "tighten"):
                return None
            driver["level"] = "tight"
            driver["switched_at"] = current["evaluation"]
            a.log.append("note", f"EDIFF tightened after evaluation {record['evaluation']}")
            return plan_next(a, driver, current["y"], "coarse", "tighten", float(parameter(a, "ediff_tight")))
        return next_coarse(a, driver)

    offset = measure_offset(a, driver, record)
    if offset is None:
        return None
    record["offset"] = {k: offset[k] for k in ("coarse_evaluation", "max_force", "max_stress", "energy")}
    if record["purpose"] == "probe":
        return next_coarse(a, driver)
    driver["finals"].append(record["evaluation"])
    record["within_limits"] = within(a, record["max_force"], record["max_stress"])
    if record["within_limits"] or driver["corrections"] >= int(parameter(a, "maximum_corrections")):
        return "publish"
    # Correct once: the coarse surface plus the offset's linear term, anchored where it was measured.
    driver["corrections"] += 1
    driver["level"] = "corrected"
    reference = driver["reference"]
    coarse = evaluation_record(a, offset["coarse_evaluation"])
    offset["gradient"] = symmetrize(reference, gradient(reference, coarse["y"], offset["forces"], offset["stress"]))
    driver["offset"] = offset
    energy, g = model(driver, coarse)
    opt = optimizer(a, driver)
    opt.current = (list(coarse["y"]), energy, g)
    store_optimizer(driver, opt)
    fmax, smax = maxima(*corrected(driver, coarse))
    driver["current"].update(max_force=fmax, max_stress=smax)
    a.log.append("note", f"final static {record['evaluation']} missed the limits; correcting the coarse surface once")
    return next_coarse(a, driver)


def finish(a: Attempt, driver: dict, state: dict) -> None:
    """Put the final static's files back at the top from its archive, archive the coarse evaluation at its geometry
    under relax/, write the report, and go on to publish. Repeating it repeats the same result."""

    last = driver["finals"][-1]
    final = evaluation_record(a, last)
    evaluations = a.workdir / DRIVER / "evaluations"
    for name in EVALUATION_ARCHIVE:
        if (evaluations / f"{last:03d}" / name).is_file():
            shutil.copyfile(evaluations / f"{last:03d}" / name, a.workdir / name)
    coarse_index = final["offset"]["coarse_evaluation"]
    archive = a.workdir / RELAX_ARCHIVE
    if archive.is_dir():
        shutil.rmtree(archive)
    shutil.copytree(evaluations / f"{coarse_index:03d}", archive)
    records = [json.loads(path.read_text(encoding="utf-8")) for path in sorted(evaluations.glob("*/record.json"))]
    reference = driver["reference"]
    end = symmetry_operations(
        final["cell"], final["frac"], atom_kinds(reference), float(parameter(a, "symmetry_tolerance"))
    )
    discrepancies = [r["correction_model_discrepancy"] for r in records if "correction_model_discrepancy" in r]
    bound = {
        "force": max((d["force"] for d in discrepancies), default=0.0),
        "stress": max((d["stress"] for d in discrepancies), default=0.0),
    }
    fraction = float(parameter(a, "correction_bound_fraction"))
    bound["within"] = bound["force"] <= fraction * float(parameter(a, "force_limit")) and bound[
        "stress"
    ] <= fraction * float(parameter(a, "stress_limit"))
    report = {
        "workflow": "vasp.relax-static-coarse",
        "optimizer": "relax_optimizer.TrustBFGS",
        "parameters": {name: parameter(a, name) for name in DEFAULTS},
        "within_limits": final["within_limits"],
        "final_evaluation": last,
        "final_max_force": final["max_force"],
        "final_max_stress": final["max_stress"],
        "coarse_evaluation_at_final_geometry": coarse_index,
        "coarse_evaluations": driver["coarse_evaluations"],
        "corrections": driver["corrections"],
        "switched_at": driver["switched_at"],
        "meshes": driver["meshes"],
        "offsets": driver["offsets"],
        "correction_model_bound": bound,
        "symmetry_operations": {"start": len(reference["symmetry"]), "final": len(end)},
        "restarts": {
            "requested": sum(1 for r in records if r["restart_requested"]),
            "observed": sum(1 for r in records if r["restart_observed"] is True),
            "unknown": sum(1 for r in records if r["restart_observed"] == "unknown"),
        },
        "interruptions": driver.get("interruptions", []),
        "incomplete_runs": [
            json.loads(path.read_text(encoding="utf-8"))
            for path in sorted((a.workdir / DRIVER / "incomplete").glob("*/record.json"))
        ],
        "evaluations": records,
    }
    (a.workdir / DRIVER / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    if len(end) < len(reference["symmetry"]):
        a.log.append("warning", f"the final geometry has {len(end)} of the {len(reference['symmetry'])} operations")
    if discrepancies and not bound["within"]:
        a.log.append("warning", f"the corrected-surface model moved from the prediction by {bound}")
    a.log.append(
        "note",
        f"final static {last} {'within' if final['within_limits'] else 'OUTSIDE'} the limits: "
        f"F {final['max_force']:.4g} eV/A, sigma {final['max_stress']:.4g} kbar, after "
        f"{driver['coarse_evaluations']} coarse evaluations and {driver['corrections']} corrections",
    )
    state.update(
        within_limits=final["within_limits"],
        coarse_evaluations=driver["coarse_evaluations"],
        corrections=driver["corrections"],
        relax_energy=evaluation_record(a, coarse_index)["energy"],
        static_max_force=final["max_force"],
        static_max_stress=final["max_stress"],
    )
    a.advance("publish", state=state)


@run.step
def publish(a: Attempt) -> None:
    """Publish the final static, the coarse evaluation at its geometry and the driver report, and complete the job."""

    base = a.parameter("data_prefix", "") or ""
    # Results stay in the persistent workdir; publish_data opts a job into a curated data/ copy too.
    to_data = a.parameter("publish_data", False) is True
    for stage_name, directory in (("static", a.workdir), (RELAX_ARCHIVE, a.workdir / RELAX_ARCHIVE)):
        prefix = f"{base}/{stage_name}" if base else stage_name
        published = []
        for name in (a.parameter("collect", COLLECT) or "").split() + (
            ["record.json"] if stage_name == "relax" else []
        ):
            if (directory / name).is_file():
                if to_data:
                    a.put(directory / name, f"{prefix}/{name}")
                published.append(name)
        if to_data:
            a.log.append("note", f"published to data/{prefix}: {', '.join(published) or 'nothing'}")
        else:
            # Without publish_data the persistent workdir is the result.
            a.log.append("note", f"kept in the workdir: {', '.join(published) or 'nothing'}")
    if to_data:
        a.put(a.workdir / DRIVER / "report.json", f"{base}/{DRIVER}/report.json" if base else f"{DRIVER}/report.json")
    oszicar = a.workdir / "OSZICAR"
    a.state["static_energy"] = last_oszicar_energy(oszicar) if oszicar.is_file() else None
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
