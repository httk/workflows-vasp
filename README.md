# workflows-vasp

This repository is the home of the *httk₂* VASP workflows, one self-contained
workflow package directory per workflow, referenced directly by a git commit.
They use the `httk.codes.vasp` helper library of an installed
*httk-workflow-vasp* (`pip install httk-workflow-vasp`).

| Directory | Short name | What it does |
| --- | --- | --- |
| `vasp-relax` | `vasp.relax` | Relax the geometry of a structure in the staged, robust *httk* v1 way, from editable INCAR templates, with the reviewed remedy ladder (Python runner). |
| `vasp-relax-bash` | `vasp.relax-bash` | The same staged relaxation, authored in Bash. |
| `vasp-static` | `vasp.static` | One single-point total-energy calculation of a fixed structure, from an editable INCAR template. |
| `vasp-relax-static` | `vasp.relax-static` | The staged relaxation of `vasp-relax`, then a single point of the relaxed structure. |
| `vasp-relax-static-coarse` | `vasp.relax-static-coarse` | A relaxation on a coarse k-point mesh driven by a BFGS optimizer in Python, VASP giving one single point per step, then one final-mesh static, correcting the coarse surface once by the measured offset. |
| `vasp-relax-static-hse` | `vasp.relax-static-hse` | An HSE relaxation from a PBE-relaxed start: calls `vasp.relax-static-coarse` with PBE (unless the structure is already PBE-relaxed), then with HSE from its result, seeded with its meshes, WAVECARs and optionally its Hessian. |
| `vasp-general` | `vasp.general` | Any VASP calculation, run on input files the job supplies; only the POTCAR may be assembled. |
| `vasp-band-structure` | `vasp.band-structure` | The electronic band structure of a structure along its HPKOT path: a self-consistent run, then a non-self-consistent line-mode run, from editable INCAR templates. |

Each package directory — its runner entry (`run.py`, or `run.sh` in
`vasp-relax-bash`) with its `templates/` — is a starting point to copy and
edit as a whole: the runner's steps spell out their work — reading the job
parameters, staging inputs, building the preparation options, running VASP
under supervision, planning and applying one remedy of the reviewed ladder, and
publishing — directly on the primitives of `httk.codes.vasp`. The
Python runners (`vasp-relax`, `vasp-static`, `vasp-relax-static`, `vasp-general`,
`vasp-band-structure`) are
independent of each other and deliberately repeat that code, templates
included; `vasp-relax-bash` is the same staged relaxation, step for step, in
Bash against the sibling Bash VASP API, with its own copies of the templates.

Each `collect.py` is just as explicit about where the results are: it names
the files its workflow's outputs come from, both where the runner leaves them
in the workdir and where `publish` puts them in published data, and reads
them with the `read_structure` and `read_total_energy` helpers of
`httk.codes.vasp.collect`, locating the files with `record.result_file`. A copied runner that keeps more results — say,
a second relaxation that archives the first one's CONTCAR — collects them by
adding one line per output to its copied `collect.py` (and declaring the output
in `httk_workflow.toml`). Running a job still requires an installed
*httk-workflow-vasp*, since the runners import `httk.codes.vasp` for inputs,
remedies, diagnostics, and collection.

## Templates

The runners other than `vasp-general` take their INCAR settings from template files, the way
*httk* v1 task templates did, so the settings are edited as INCAR files rather
than as code or job parameters. Each stage of a workflow has one template in
its package's `templates/` directory:

| Stage | Template | What it is |
| --- | --- | --- |
| `preprerelax` | `INCAR.preprerelax` | A cheap, coarse relaxation that may also change the cell shape (`ISIF = 7`). |
| `prerelax` | `INCAR.prerelax` | A loose relaxation of positions and cell (its own `EDIFF = 1E-3`, `EDIFFG = -0.3`) that writes the charge density. |
| `relax1`, `relax2`, … | `INCAR.relax` | Accurate relaxation rounds; from `relax2` on, each reads the previous round's charge density. |
| `static` | `INCAR.static` | A single point with the electronic settings of `INCAR.relax` and the ionic loop off, reading the final relax round's charge density when there is one. |
| `scf` | `INCAR.scf` | `vasp-band-structure`'s self-consistent stage: `INCAR.static` from atomic charges (`ICHARG = 2`), writing the charge density (`LCHARG = .TRUE.`). |
| `bands` | `INCAR.bands` | `vasp-band-structure`'s non-self-consistent band stage (`ICHARG = 11`, `ISMEAR = 0`), with the density-defining settings of `INCAR.scf`. |

`vasp-relax` and `vasp-relax-bash` run `preprerelax`, `prerelax` and then relax rounds: from `relax2`
on, a round whose final energy is within ten times `EDIFFG` of the previous
round's ends the relaxation, and a round that ran out of `NSW` ionic steps, or
moved further, is followed by another. The pre-stages and `relax1` also accept a
run that only ran out of ionic steps. Each finished stage is archived in the
workdir under `stages/<stage>/` and its CONTCAR becomes the next stage's POSCAR;
the final round stays at the top of the workdir. `vasp-relax-static` runs the
same stages, archives the final round under `relax/`, and then runs `static` at
the top of the workdir. `vasp-static` runs `static` alone. The relaxation
templates and their stage logic are those of the *httk* v1
`vasp-robust-relax-formenrg` template, version 1.0 by Christopher Tholander.

Every stage prepares its inputs the same way. The template is the starting
INCAR: it is copied into the workdir, and preparation then rewrites the tags it
derives and the job's `incar_tags`, so those lines (an `EDIFF` and `EDIFFG`
included) may move to the end of the INCAR VASP runs with. An `EDIFF` and
`EDIFFG` the template sets are kept, while absent ones are
derived from `accuracy_per_atom` (`EDIFF = accuracy·atoms/margin` with a margin
of 500 when the template's `NSW` is above 1 and 33 otherwise, and
`EDIFFG = accuracy·atoms`). `MAGMOM` and `NBANDS` are derived when absent, and
the KPOINTS grid comes from the stage's k-point density. Every remedy the job
has already applied is re-applied, so a remedy found necessary in one stage
stays in force in the later ones, as in v1. Finally, a template's
`ICHARG = 1` stays only when the workdir holds a non-empty CHGCAR from a run
with the same `PREC`, `ENCUT` and `ENAUG` (compared case-insensitively, an
absent tag equal only to an absent one), and becomes `ICHARG = 2` otherwise.
So `relax1` starts from atomic charges, since `prerelax` ran at `PREC = Normal`
and the default cutoff, while later relax rounds and the static stage reuse
the density.

There are three ways to change the settings:

* `incar_tags` (a job parameter) — tags applied to every stage, winning over
  the templates and every derived value.
* `files/INCAR.<stage>` (a payload file of one job, e.g. `files/INCAR.relax`) —
  replaces that stage's template wholesale for that job.
* Edit the templates of a copied package — the v1 way: copy the package
  directory, edit its `templates/INCAR.*`, and install the copy with
  `--workflow-dir`. A runner installed as a single file
  (`httk job new --from-runner .../run.py`) has no `templates/` and fails with
  `vasp.template_missing`; reference the package directory instead (a git URI
  ending in `#vasp-relax`, or `--workflow-dir`). `vasp-relax-bash` and `vasp-relax-static` carry their own
  copies of the `vasp-relax` (and `vasp-static`) templates; this repository's tests keep them
  identical, and a copied package is free to diverge.

### Differences from httk v1

* **ICHARG.** v1's templates ask relax to read the previous charge density, but
  v1's own test (`[ ! -e CHGCAR -o -s CHGCAR ]`) set `ICHARG = 2` whether or not
  a CHGCAR was there, so v1 never read one. Reading prerelax's density into
  relax1 would also put it on a different FFT grid, which was never validated.
  The rule above reads a density only from a run on the same grid.
* **Remedies carry over the way the ladder records them.** Every `bump_kpoints`
  step adds one to the grid again in every later stage (v1 kept a single flag),
  and `scale_ediff` scales a template's own `EDIFF` too.
* **No adjustments after a successful stage.** v1 also adjusted the inputs after
  a stage that succeeded (for example DENTET → `ISMEAR = 0`, k-point shifts →
  Gamma centering); here those diagnostics are warnings and change nothing.
* **Electronic nonconvergence is remedied.** A stage whose electronic loop did
  not converge goes to the remedy ladder and restarts from that stage's starting
  POSCAR, where v1 accepted it like an ionic non-convergence.
* **No `fix_too_close` stage.** The ladder's `ions_too_close` remedy, which
  scales the lattice by 5 %, replaces v1's separate stages, which scaled it by
  10 % in the pre-stages and by 5 % in relax1.

## vasp-relax-static-coarse

`vasp-relax-static-coarse` is `vasp-relax-static` with its staged IBRION
relaxation replaced by a driver loop in the runner. VASP only supplies energy,
forces and stress. Each evaluation is a fresh single point (`NSW = 0`) at a
geometry chosen by the optimizer, restarting from the WAVECAR of the previous
evaluation on the same mesh. The optimizer is `relax_optimizer.py`, a damped
BFGS with a trust radius, standard library only. It runs on 3N atomic
displacements plus 6 cell strains. The symmetry operations of the starting
structure are found once (`symmetry_tolerance`), and every gradient and step
is projected onto the coordinates they leave invariant, so the path keeps the
starting symmetry whatever VASP's own symmetrization does. The report gives the
number of operations at the start and at the final geometry.

Stages:

1. **Coarse, loose.** Evaluations on the coarse mesh at `ediff_loose`, until
   the accepted point is within `switch_force` / `switch_stress`.
2. **Coarse, tight.** The switch point is evaluated again at `ediff_tight`,
   and that evaluation *replaces* the optimizer's accepted point without an
   acceptance test or a curvature update, so energies and gradients at two
   EDIFF levels are never compared (its `ediff_shift` records what the tighter
   EDIFF changed). Then evaluations at `ediff_tight` until an accepted point is
   within `force_limit` / `stress_limit`.
3. **Final static.** One final-mesh static of that geometry. Its forces and
   stress minus those of the tight coarse evaluation at the same geometry (the
   pairing is checked; otherwise the job fails with `vasp.offset_pairing`) give
   the offset. If the static misses the limits, the coarse relaxation continues
   on a corrected surface, and a second final static ends the job
   (`maximum_corrections`).

The corrected surface is one model: the coarse energy plus the linear term
`g_offset · (y − y_anchor)`, with the constant gradient `g_offset` (the
offset's gradient where it was measured). The stop test uses the physical
prediction, coarse forces and stress plus the measured offset. Each corrected
evaluation records how far the two differ (`correction_model_discrepancy`,
the largest force in eV/Å and stress component in kbar at that evaluation's
geometry, in the symmetric subspace the optimizer works in), and the report
compares the largest with `correction_bound_fraction` of the limits
(`correction_model_bound`). It is a diagnostic, logged as a warning when out
of bound; only the final static's own forces and stress decide.

The job succeeds once its last final static has run. `within_limits` in the
job state and the report says whether that static met the limits, and
`collect.py` harvests the structure and energy only when it did: an
out-of-limits result stays in the workdir and its report. More than
`maximum_evaluations` coarse evaluations fail the job with
`vasp.relax_not_converged`, also when the budget runs out just as the EDIFF is
tightened (the tight repeat of the switch point is never skipped, so the cap is
strict); a proposed cell that is singular, inverted or
changes the volume by more than a factor 2 fails it with `vasp.cell_invalid`.

**Time guard.** An evaluation starts only if the attempt's time left
(`Attempt.deadline`) covers `time_factor` times its expected VASP time plus the
largest staging and archiving overhead measured so far. It is checked before
staging and again right before VASP starts. The expected time is the slowest
evaluation so far of the same kind (coarse at the loose EDIFF, coarse at the
tight one, final); before the first of a kind it is `expected_seconds` (first
coarse), `tight_cost_factor` times the loose one (first tight), or
`final_cost_factor` times the slowest coarse one (first final). Without
`expected_seconds` the first evaluation runs unguarded and the log says so;
give it for any run with a time limit. An evaluation without the time is
published again with that `mintime`, so a manager with less time left does not
claim it and a later one resumes. When the time needed exceeds the attempt's
own `maxtime`, no attempt could run it: the job fails with
`vasp.insufficient_time` (`needed_seconds`, `maxtime_seconds`) after
committing what it has, and is continued once `maxtime` is raised.

**Runs that do not complete** are classified by how the process ended (the
supervised process's `termination`), never by the clock alone. Ended by a
signal to the runner (`signal_15`, also `…_killed`) at the attempt's deadline
(within 5 s of it, or later): the manager stopped it at its maxtime or drain;
no remedy is used, the evaluation waits for a manager with 1.4 times the time
the stopped run took, and the relaxation resumes from the last committed
evaluation. At most `maximum_interruptions` such deferrals per job; one more
fails it with `vasp.interrupted_repeatedly`. Ended by a signal well before the
deadline (a cancellation, or a stop from outside httk): the job fails with
`vasp.interrupted` (a cancelled job is committed as cancelled whatever it
publishes). Ended by the runner's own `timeout`: the job fails with
`vasp.evaluation_timeout`, without a remedy. Any other ending (a nonzero exit,
a diagnosed stop, nonconvergence), however close to the deadline, goes to the
remedy ladder. Every such run keeps its files and a `record.json` (how it
ended, how it was handled, its diagnostics) under
`relax-driver/incomplete/NNN-KK/`. A runner killed outright, or stopped
outside the VASP run, publishes nothing; the manager records the failure and
the snapshots make a retry safe.

**Remedies.** Only remedies that change how VASP reaches the electronic ground
state are applied (`ICHARG`, `AMIX`, `BMIX`, `AMIN`, `ALGO`, `TIME`, `NPAR`,
`NEDOS`); they stay applied in later evaluations, and a remedy's `ICHARG`
reset is kept even when a WAVECAR is read. A remedy that would change the
surface or the structure (k-points, centering, smearing, symmetry, `LREAL`,
bands, EDIFF, the lattice) stops the job with `vasp.unsupported_remedy` rather
than compare evaluations of different surfaces; so does a job remedy history
that holds any other operation (checked before an evaluation's inputs are
written, so an older or edited history cannot bring one in). `rattle_amplitude` is not used.

**Restarts and recovery.** A WAVECAR is cached as
`relax-driver/<kind>/WAVECAR-NNN` (written by evaluation NNN) with the inputs
it was written for (k-points and `ENCUT`, `PREC`, `NBANDS`, `ISPIN`, ...), and
staged only for the same inputs; a payload seed is staged as it is. A failed
run works on a copy, so the cache stays intact. Each record says whether a
restart was requested (`restart_requested`, `restart_source`) and whether
VASP's stdout shows the WAVECAR was read (`restart_observed`: true for "the
WAVECAR file was read successfully", false for "WAVECAR not read", otherwise
`"unknown"`). The optimizer state after evaluation NNN is
`relax-driver/driver-NNN.json`, written after the evaluation's record; it names
the cached WAVECAR of each kind, and a WAVECAR no snapshot names is deleted. A
retried attempt continues from the job state's snapshot or from a later one an
interrupted attempt already wrote. After the last final static the snapshot
says only publication is left, and a retry republishes from the archive
without running VASP.

**Provenance.** Each evaluation is archived with its inputs and outputs and a
`record.json` under `relax-driver/evaluations/NNN/`. The record holds the
stage, mesh, EDIFF, PRECFOCK, ICHARG and remedies in force, energy, forces,
stress, SCF steps, VASP and staging seconds, the restart evidence, and whether
the offset was applied. `relax-driver/report.json` gathers every record and
every measured offset, the correction-model bound, the symmetry counts, a
restart summary, the interruptions and the records of incomplete runs, and is published with `publish_data`. `relax/` holds the
coarse evaluation at the final geometry and the final static is at the top, the
layout of `vasp-relax-static` that `collect.py` reads.

The INCAR templates are `templates/INCAR.coarse` and `templates/INCAR.final`
(payload overrides `files/INCAR.coarse` and `files/INCAR.final`); they are PBE,
and an HSE run gives its own (as `vasp-relax-static-hse` does). After the
templates and every tag, the runner sets the tags that make each run one
evaluation: `NSW = 0`, `IBRION = -1`, `ISIF = 2`, `LWAVE = .TRUE.`, the
stage's `EDIFF`, and `ISTART`/`ICHARG` from whether a fitting WAVECAR is there.
The general job parameters above apply, except `maximum_relax_rounds`,
`prerelax_kpoint_density` and `rattle_amplitude`; `centering` defaults to
`Gamma` here, as tetrahedron integration (`ISMEAR = -5`) needs. The driver's
own parameters are listed below. Their defaults come from the HSE campaign
that motivated this workflow (KSnF3/KGeF3; for the measurements behind them,
see `notebook/retro/owen-replay-2026-10-17.md` and
`notebook/retro/owen-fixes-2026-10-17.md` in the researchbotter project):

| Parameter | Default | Meaning and origin |
| --- | --- | --- |
| `coarse_mesh`, `final_mesh` | empty | Explicit `"n1 n2 n3"` meshes. When empty, the mesh is derived once from `coarse_kpoint_density` (10) or `kpoint_density` (20) and then kept for the whole relaxation. The campaign used 5³ coarse and 6³ Γ final for KGeF3 TET. KSnF3 TET used 6³ for both, with only PRECFOCK coarser. |
| `coarse_incar_tags`, `final_incar_tags` | `{"PRECFOCK": "Fast"}`, `{"PRECFOCK": "Normal"}` | Tags for one kind of evaluation, applied over `incar_tags`. |
| `ediff_loose`, `ediff_tight`, `final_ediff` | 1e-6, 1e-8, 1e-8 | At 1e-6, repeated evaluations at the same geometry differed by up to 0.002 eV/Å and 0.28 kbar, which is the size of the limits. At 1e-8 they differed by ≤ 0.0007 eV/Å and 0.017 kbar. |
| `switch_force`, `switch_stress` | 0.02 eV/Å, 1 kbar | Where the tight EDIFF starts. |
| `force_limit`, `stress_limit` | 0.003 eV/Å, 0.2 kbar | The campaign's acceptance limits (largest atomic force, largest stress component). |
| `time_factor` | 1.4 | Time-guard factor on the expected VASP time. |
| `expected_seconds` | 0 (none) | Expected VASP seconds of the first coarse evaluation; without it the first evaluation is unguarded. |
| `tight_cost_factor` | 2 | Expected cost of the first tight-EDIFF coarse evaluation relative to the slowest loose one, until one has run. Not yet measured; deliberately generous. |
| `final_cost_factor` | 8 | Expected cost of the first final-mesh evaluation relative to the slowest coarse one. It is used only before any final-mesh evaluation has run. It is a rough estimate: mesh 5³ → 6³ cost 3.4× per SCF step, and Normal cost 2.3–2.8× Fast. |
| `maximum_evaluations`, `maximum_corrections` | 60, 1 | Cap on coarse evaluations (strict, including the tight repeat), and the number of offset corrections. |
| `maximum_interruptions` | 3 | Evaluations deferred because the manager stopped VASP at the attempt's deadline; one more fails the job. Every deferral is a new activation, which no other driver cap bounds. |
| `trust_atom`, `trust_strain`, `rise_tolerance` | 0.03 Å, 0.01, 2e-4 eV | Largest step per atom (physical displacement in the current cell) and per strain component (relative to the current cell) at trust factor 1. An evaluation whose energy rises more than `rise_tolerance` above the accepted point at the same EDIFF is rejected and the trust factor halves. |
| `hessian_atom`, `hessian_bulk_modulus` | 70 eV/Å², 50 GPa | The isotropic initial Hessian model; both must be positive. In the offline replay, scaling it by ×0.5 to ×2 changed the KGeF3 TET count from 15 to 19 evaluations. |
| `hessian_seed` | empty | Optional payload JSON `{"hessian": [[...]]}`, (3N + 6)², in the driver's coordinates: Cartesian displacements in Å, then the √2-weighted Voigt strain times V^(1/3). It must be finite, positive definite and invariant under the structure's symmetry operations (checked), or the job fails with `vasp.hessian_seed_invalid`. Off by default: the PBE Hessians fitted from the campaign's trajectories were indefinite. |
| `hessian_seed_symmetrize` | false | Project the seed onto the structure's symmetry (average of P H Pᵀ over the operations, logged with the largest change) instead of refusing a seed that is not invariant. `vasp-relax-static-hse` sets it for a PBE Hessian, learned along a path that may have started at lower symmetry. |
| `symmetry_tolerance` | 0.001 Å | Tolerance of the symmetry search (operations with integer matrices of entries −1, 0, 1 in the given cell). |
| `correction_bound_fraction` | 0.1 | The corrected-surface model and the physical prediction are reported as within bound when they differ by at most this fraction of the limits. |
| `offset_probe` | false | Adds a final-mesh evaluation at the tight re-evaluation of the switch point, so the offset is measured at two geometries. |
| `coarse_wavecar`, `final_wavecar` | `files/WAVECAR`, `files/WAVECAR.final` | Optional seed WAVECARs for the first evaluation of each kind. |

## vasp-relax-static-hse

`vasp-relax-static-hse` runs no VASP itself. It declares
`[workflow.calls] coarse = "vasp.relax-static-coarse"` and runs that driver as
child jobs (`Attempt.call`, then `Attempt.gather`), once per functional, so the
optimizer, the offset correction, the time guard and the remedies are the
driver's own. Installing it installs the driver too, and a manager does not
claim it until the driver is installed in the workspace.

- `mode = "no_pbe"` (default): child `pbe` runs the driver on the input
  structure with the driver's own (PBE) templates; then child `hse` runs it
  from the PBE result's `relax/CONTCAR`.
- `mode = "start_from_pbe"`: the input structure is PBE-relaxed and only the
  `hse` child runs. `pbe_provenance` (where the structure came from) is
  required, or the job fails with `vasp.pbe_provenance_missing`. Optional
  seeds in the payload: `files/WAVECAR.pbe`, `files/WAVECAR.pbe.final` (used
  only when the HSE stage's mesh of that kind is set), `files/hessian.pbe.json`.

The HSE child gets the PBE child's coarse and final meshes unless the job sets
them, so the WAVECARs fit; with `seed_wavecar` (default true) the PBE child's
last coarse- and final-mesh WAVECARs as its first restarts; with
`hessian_from_pbe` (default false) the PBE child's last BFGS Hessian as its
`hessian_seed`. The HSE INCARs are `templates/INCAR.hse.coarse` and
`templates/INCAR.hse.final` (the driver's templates plus ALGO Damped, LHFCALC,
AEXX 0.25, HFSCREEN 0.2, TIME 0.4, ISYM 3), overridden by
`files/INCAR.hse.coarse` / `files/INCAR.hse.final`; the PBE child takes
`files/INCAR.pbe.coarse` / `files/INCAR.pbe.final` when given. Every other job
parameter goes to both children, then `pbe_parameters` / `hse_parameters` (job
parameter dictionaries) on top.

`stages.json` records per stage the child job key and id, the called
workflow's installed id, the structure, INCAR and seed sha256, meshes,
evaluation counts, final energy, force and stress, and `within_limits`, and the
lattice and volume shift of the HSE stage from its start. `after_hse` copies
the HSE child's `relax/` and final static into the workdir in the driver's
layout, plus its report as `hse-report.json`, and `collect.py` reads them
there; the children collect as child runs of their own. A child that fails
fails the job with `vasp.pbe_stage_failed` or `vasp.hse_stage_failed`, its
failure recorded in `stages.json`. A PBE child that ends outside its limits is
still used as the HSE start.

Limits: a called child's `maxtime` is capped at the calling attempt's, so do
not give this job a `maxtime` below one HSE evaluation; the first HSE
evaluation has no cost estimate yet, so the driver's time guard does not
cover it; all children travel with the job (an exchange eject returns the
whole tree, WAVECARs included).

## vasp-general

`vasp-general` runs any VASP calculation: the job supplies every input file
but the POTCAR, and VASP runs on them exactly as given. Nothing is derived, so
INCAR and KPOINTS reach VASP byte for byte. A job payload looks like:

```
files/INCAR      required
files/POSCAR     required (or the structure input, which lands here too)
files/KPOINTS    optional: VASP can use KSPACING instead
files/POTCAR     optional: assembled from pseudopotential_library otherwise
files/WAVECAR    optional, like any other file VASP reads (CHGCAR, ...)
```

Every regular file directly under `files/` is copied into the workdir;
subdirectories are not, nor are `vasp-run-report.json` and
`POTCAR.provenance.json`, which are the runner's own bookkeeping (a note in the
job log names any it skipped). The `structure` input is optional and is staged as
`files/POSCAR`, so it and a supplied `files/POSCAR` are one file and a job gives
one of them. A missing INCAR, POSCAR or POTCAR (with no library to assemble
one) fails the job with `vasp.input_missing`.

Remedies are off by default (`maximum_remedies = 0`): a run that does not
complete fails the job at once with `vasp.failed`, and the failure details
carry the remedy the reviewed ladder would have applied. Set
`maximum_remedies` above zero to let the runner apply that many remedies and
rerun. The collected outputs are `final_structure` (the CONTCAR, when the
calculation wrote one) and `total_energy` (from the OUTCAR); an output whose
file is missing is reported as unfulfilled rather than failing the collection.
A directory holding the input files is staged as one job's `files/` with:

```console
httk job new --install --workflow-dir vasp-general --files my-calculation/
```

## vasp-band-structure

`vasp-band-structure` computes the electronic band structure of a structure
that is already chosen, normally the `relaxed_structure` of a `vasp.relax` job
chained before it; it does not relax. `prepare` recognizes the structure's
space group within `symmetry_tolerance` with `httk.atomistic.band_path`, and
replaces it by its HPKOT primitive cell (the convention of SeeK-path): the
band k-points refer to that cell's reciprocal basis, so both VASP runs compute
it, and it is the collected `band_structure_cell`. Recognition symmetrizes:
sites may move within the tolerance, so the cell is the standardized structure,
not merely another basis of the input. The written POSCAR keeps the input
POSCAR's species order, so a supplied POTCAR and `LDAU*` arrays keep their
meaning. `prepare` also writes the line-mode `KPOINTS.bands` and
`band_path.json`, the exact labels, coordinates (as fractions) and path, next
to the rounded KPOINTS. Its `transform` is the integer basis change from the
`httk.atomistic.primitive_cell` of the input to the band cell, not from the
input POSCAR's own basis.

The `scf` stage (`INCAR.scf`) runs on the automatic mesh of `kpoint_density`
and writes the CHGCAR; it is archived under `scf/` and its CHGCAR stays at the
top of the workdir. The `bands` stage (`INCAR.bands`) then runs there with
`ICHARG = 11` and `LCHARG = .FALSE.` (it must not overwrite the scf density)
on `KPOINTS.bands`. Both stages share the job's remedy budget and
history, so a remedy of the scf stage is re-applied in the band stage, except
that the band stage always runs with `ICHARG = 11`. Remedies that move the
structure (`scale_lattice`, `contcar_to_poscar`) are never applied, nor, in the
band stage, an `ICHARG` change; the ladder goes on to its next rung instead.
`NBANDS` may differ between the stages (more bands show more of the
conduction range); its derived default is the same in both.

The job fails by name when (the first five in `prepare`, before any VASP run,
where the INCAR checks read `incar_tags` and each stage's payload INCAR, or its
template when there is none):

* `vasp.cell_dependent_tag_unsupported` — `MAGMOM`, `NELECT`, `NUPDOWN` or an
  absolute FFT grid (`NGX NGY NGZ NGXF NGYF NGZF`) is set, or the POSCAR
  comment has a `[MAGMOM=...]`: these are per site, per cell or sized for the
  cell, and the HPKOT cell's sites and size differ from the input's.
* `vasp.band_path_failed` — the structure has no band path: site moments, an
  irrational metric, or spglib is missing (`pip install httk-atomistic[default]`).
* `vasp.unsupported_calculation` — `LHFCALC`, `METAGGA`, `LNONCOLLINEAR` or
  `LSORBIT` is set (checked again on each stage's prepared INCAR). The scope is
  collinear LDA/GGA(+U): a hybrid functional needs zero-weight explicit
  k-points in a self-consistent run, and a meta-GGA the kinetic-energy density,
  which `ICHARG = 11` does not provide.
* `vasp.invalid_band_inputs` — `incar_tags` sets `ICHARG` at all (they apply
  to both stages, so any value is wrong for one), the band INCAR sets it to
  anything but 11, or `band_divisions` is not an integer of at least 2; also
  when the band INCAR does not read `ICHARG = 11` at a launch.
* `vasp.input_mismatch` — the HPKOT cell cannot be written in the input's
  species order.
* `vasp.charge_density_missing` — the band stage finds no non-empty CHGCAR.
* `vasp.charge_density_incompatible` — the band INCAR differs from the scf
  run's in a density-defining tag (`PREC ENCUT ENAUG NGX NGY NGZ NGXF NGYF NGZF
  ISPIN MAGMOM LMAXMIX GGA METAGGA LASPH LDAU LDAUTYPE LDAUL LDAUU LDAUJ NELECT
  LNONCOLLINEAR LSORBIT GGA_COMPAT LUSE_VDW AGGAC AGGAX ALDAC PARAM1 PARAM2
  LMAXPAW XC`, compared case-insensitively, an absent tag equal only
  to an absent one); the details list each difference.

The outputs are `band_structure_cell` (the band-stage POSCAR), `total_energy`
(of the scf run, from `scf/OUTCAR`), and two `files` entries: `band_eigenvalues`
(the band run's EIGENVAL) and `band_path` (`band_path.json`). With
`publish_data` the stages are published under `scf/` and `bands/` below
`data_prefix` (default empty), always including the files collect reads. This
runner imports *httk-atomistic*, and recognition needs spglib, which the
manifest requires.

## Job parameters

Every packaged VASP runner reads the same `inputs` object, and every member is
optional. Paths are relative to the job payload; lists are space-separated
strings, so the Bash runner and the Python runners read one contract.

* `poscar` (default `files/POSCAR`) — the starting structure; any VASP-5
  POSCAR or CONTCAR file will do.
* `potcar` (default `files/POTCAR`) — a pre-assembled POTCAR. When the file
  is absent and `pseudopotential_library` is set, the POTCAR is assembled per
  species and a provenance record is written next to it.
* `pseudopotential_library` (default none) — root of a VASP pseudopotential
  library, one directory per variant.
* `kpoint_density` (default `20.0`), `centering` (default `Monkhorst-Pack`),
  `accuracy_per_atom` (default `0.001`) — passed to
  `httk.codes.vasp.inputs.VaspPreparationOptions`. `kpoint_density` is the
  density of the relax rounds and the static stage, and of the scf stage of
  `vasp-band-structure`.
* `prerelax_kpoint_density` (default `10.0`) — the k-point density of the
  `preprerelax` and `prerelax` stages.
* `maximum_relax_rounds` (default `5`) — a relax round numbered above this that
  still has not converged fails the job with `vasp.relax_not_converged`; with
  the default the last round is `relax6`, as v1's `RELAXSTEP > 5`.
* `parallel_tag` and `parallel_value` (default none) — one of `NPAR`, `NCORE`,
  or `KPAR`, and its value.
* `incar_tags` (default empty) — explicit INCAR tags for every stage. They win
  over the templates and every derived value.
* `timeout` (default `86400`) — seconds one VASP execution may take before its
  process group is terminated.
* `maximum_remedies` (default `8`, and `0` in `vasp-general`) — how many
  remedies this job may apply in total, over all its stages, before it fails.
  The ladder is bounded per problem as well.
* `remedy_policy` (default `reviewed-v1`) — the registered remedy policy the
  runner plans with, so a group with its own reviewed practice registers a
  policy with `httk.codes.vasp.register_remedy_policy` and names it here
  instead of editing a runner.
* `rattle_amplitude` (default `0.0`) — when positive, the POSCAR is rattled by
  this amplitude after every applied remedy, with a seed derived from the
  attempt, so two retries never repeat one structure. `vasp-band-structure`
  does not use it: its structure is the HPKOT cell its k-points refer to.
* `band_divisions` (default `40`, `vasp-band-structure` only) — k-points per
  segment of the line-mode band path.
* `symmetry_tolerance` (default `1e-3`, in Å, `vasp-band-structure` only) —
  the tolerance the structure's symmetry is recognized within, also the
  precision its POSCAR is read with; recorded in `band_path.json`.
* `publish_data` (default `false`) — when `true`, `publish` also copies the
  `collect` files into the job's `data/` directory; by default outputs stay in
  the workdir only.
* `collect` (default `INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml
  vasp-run-report.json POTCAR.provenance.json`) — space-separated file names
  copied to `data/` only when `publish_data` is `true`. `vasp-relax-static`
  also uses this list to archive the relaxation before the static stage.
  `vasp-band-structure` defaults to `INCAR KPOINTS POSCAR OUTCAR OSZICAR
  EIGENVAL vasprun.xml band_path.json vasp-run-report.json
  POTCAR.provenance.json` and also archives these into `scf/`.
  Missing files are skipped. The `stages/` archives are not published.
* `data_prefix` (default `vasp`) — directory below the job's data the
  collected files are published under when `publish_data` is `true`;
  ignored for workdir results. `vasp-relax-static` defaults to an empty prefix
  and publishes its stages under `relax/` and `static/`, and
  `vasp-band-structure` under `scf/` and `bands/`.
* `vasp_command` (default empty) — the VASP command as one argv string, split
  the way a shell splits it. It names the program (`vasp_std`); the run helper
  prepends the attempt's launch prefix (`HTTK_WORKFLOW_LAUNCH`), so do not put
  `srun` or `mpirun` in it. A leftover launcher (for example `vasp.command = "srun -n 32 vasp_std"`) is refused by `run_vasp` with a `ValueError` when a prefix applies; the Python workflows catch only `OSError`, so the attempt ends as a runner error whose message explains the fix, and the Bash runner reports "could not run VASP at all (status 2)" with the explanation on stderr. The environment variable `HTTK_VASP_COMMAND`
  overrides it, which is how a deployment — or a test — chooses the executable
  without touching any job.

The workdir (`run/`) persists across attempts, so the inputs a remedy rewrites
are the inputs the next attempt reads. By default the persistent workdir is
the result, with no `data/` copy, and collectors read it directly. Pass
`--parameter publish_data=true` to `httk job new` (or
`parameters={"publish_data": True}` to `new_job`) to also publish the curated
files into `data/`.

## Declaration identity

Each package's `declaration.json` is identified by its `$id`, the GitHub raw
URL of that file at a `wf-decl-vX.Y.Z` tag, and `httk_workflow.toml` names the
same URL as `declaration_uri`. The version is the declaration's own
`x-httk-definition.version`. Declaration versions across the repository form
one increasing sequence shared with the tags. A `wf-decl-*` tag is created only
when a declaration's content changes and is never moved. To change a
declaration, edit the document and give it a version higher than every existing
`wf-decl-*` tag (normally the next minor), set `$id` and `declaration_uri` to
that tag's URL, commit, tag `wf-decl-v<version>`, and push the tag with
`git push origin wf-decl-v<version>` (a plain `git push` does not push tags).
Declarations that did not change keep their version and `$id`, so versions may
skip numbers. The URL names the declaration document, not the workflow code,
which jobs reference by git commit (next section). The earlier
`https://schemas.httk.org/defs/v0.1/workflows/…` identities remain valid for
runs already stored under them.

## Reference a workflow by commit

A workflow is installed in a workspace before jobs of it are created;
`--install` installs it first (`httk workflow install --workspace WS SOURCE`
does it on its own):

```console
httk job new --install --workflow 'git+https://github.com/httk/workflows-vasp@<ref>#vasp-relax' --input structure=POSCAR
```

`<ref>` may be a commit, branch, or tag, or omitted for the default branch; it
is canonicalized to the full commit hash the repository is cloned at, and the
named subdirectory's `httk_workflow.toml` package is installed. Once a workflow
has been installed this way, its short name (e.g. `vasp.relax`) also resolves.

## Install as a plugin

```console
httk plugin install 'git+https://github.com/httk/workflows-vasp'
```

installs all the workflow packages listed in `httk_plugin.toml` at once on
this machine; each is still installed into a workspace (`httk job new
--install --workflow vasp.relax`) before its jobs are created.

## Install with pip

From a clone, `pip install -e .` or `pip install .` installs the workflows'
dependencies and makes every workflow of this repository known to *httk₂* in
that Python environment, at HEAD and at every `wf-v*` tag, each under its
GitHub commit URI. Only committed content is installed; use a package path
for a working copy. An editable install follows `git pull` (restart running
httk processes); a non-editable install at an unchanged version picks up new
tags only with `pip install --force-reinstall --no-deps .`. Shallow clones are
refused (`git fetch --unshallow --tags`).
