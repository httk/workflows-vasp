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
| `vasp-relax-ions` | `vasp.relax-ions` | One accurate ionic relaxation of a structure in its fixed cell (`ISIF = 2`) or at fixed volume (`ISIF = 4`), from an editable INCAR template, with the reviewed remedy ladder. |
| `vasp-dos` | `vasp.dos` | The electronic density of states and band gaps of a structure: a self-consistent run (or a given charge density), then a non-self-consistent run on a dense uniform mesh. |
| `vasp-dielectric` | `vasp.dielectric` | The static and high-frequency dielectric tensors, Born effective charges and piezoelectric tensor of an insulator by density-functional perturbation theory (one `IBRION = 8` run). |
| `vasp-elastic` | `vasp.elastic` | The elastic tensor and moduli of a structure: an optional tight relaxation, then one `vasp.relax-ions` child per strained cell, fitted with *httk-analyse*. |
| `vasp-eos` | `vasp.eos` | The Birch–Murnaghan equation of state of a structure: an optional tight relaxation, then one fixed-volume `vasp.relax-ions` child per volume (optionally a `vasp.static` each), fitted with *httk-analyse*. |
| `vasp-phonons` | `vasp.phonons` | Harmonic phonons of a structure with phonopy: one `vasp.static` child per displaced supercell, then force constants, band structure, DOS and thermodynamics, with an optional non-analytical correction. |
| `vasp-full-characterization` | `vasp.full-characterization` | Characterize one material with the other packages as child jobs: one tight relaxation, then band structure, density of states, elastic tensor, equation of state, dielectric response (insulators) and phonons, sharing every calculation. |

Each package directory — its runner entry (`run.py`, or `run.sh` in
`vasp-relax-bash`) with its `templates/` — is a starting point to copy and
edit as a whole: the runner's steps spell out their work — reading the job
parameters, staging inputs, building the preparation options, running VASP
under supervision, planning and applying one remedy of the reviewed ladder, and
publishing — directly on the primitives of `httk.codes.vasp`. The
Python runners (`vasp-relax`, `vasp-static`, `vasp-relax-static`, `vasp-general`,
`vasp-band-structure`, `vasp-relax-ions`) are
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
| `nscf` | `INCAR.nscf` | `vasp-dos`' non-self-consistent DOS stage: its `INCAR.scf` (a copy of `vasp-band-structure`'s) with `ICHARG = 11`, `LCHARG = .FALSE.`, `ISMEAR = -5`, `LORBIT = 11` and `NEDOS = 3001` (an httk choice), as atomate2's `NonSCFSetGenerator` in uniform mode. |
| `relaxions` | `INCAR.relaxions` | `vasp-relax-ions`' fixed-cell relaxation (`ISIF = 2`, `IBRION = 2`, `NSW = 99`): the electronic settings of `INCAR.relax` made tight as atomate2's `ElasticRelaxMaker` (`EDIFF = 1E-7`, `EDIFFG = -0.001`, `LREAL = .FALSE.`), with `ISMEAR = 0`, `SIGMA = 0.05`, from atomic charges, writing no CHGCAR or WAVECAR; `NELM = 100` is an httk choice. |
| `dielectric` | `INCAR.dielectric` | `vasp-dielectric`'s perturbation-theory run (`IBRION = 8`, `LEPSILON`, `LPEAD`, `NSW = 1`, `ISMEAR = -5`, `LREAL = .FALSE.`, `EDIFF = 1E-7`): `INCAR.static`'s electronic settings with atomate2's `lepsilon` choices; insulators only. |
| `phonon` | `INCAR.phonon` | `vasp-phonons`' displaced-supercell single points: `INCAR.static` with atomate2 `PhononDisplacementMaker`'s `EDIFF = 1E-7`, `LREAL = .FALSE.` and `LCHARG = .FALSE.`, and the httk choices `ISMEAR = 0`, `SIGMA = 0.05`, `ADDGRID = .TRUE.`, from atomic charges. |

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
`EDIFFG = accuracy·atoms`). `MAGMOM` and `NBANDS` are derived when absent,
unless `incar_tags` sets them to `null`, which leaves them to VASP's own
defaults. A `null` `NBANDS` may still end up explicit: at launch the VASP
launch wrapper writes VASP's default rounded to the band-group count when the
launch has more than one band group and the rounded default stays within
`2 × NELECT`, unless the run is one it leaves to VASP (DFPT, hybrid, GW-type,
`IMAGES`) or the setting `vasp.nbands_rounding` is `off` (see
*httk-workflow-vasp*'s helper details). The KPOINTS grid comes from the stage's k-point density. Every remedy the job
has already applied is re-applied, so a remedy found necessary in one stage
stays in force in the later ones, as in v1. Finally, a template's
`ICHARG = 1` stays only when the workdir holds a non-empty CHGCAR from a run
with the same `PREC`, `ENCUT` and `ENAUG` (compared case-insensitively, an
absent tag equal only to an absent one), and becomes `ICHARG = 2` otherwise.
So `relax1` starts from atomic charges, since `prerelax` ran at `PREC = Normal`
and the default cutoff, while later relax rounds and the static stage reuse
the density.

### Changing INCAR tags

There are several ways to change individual INCAR settings:

* `incar_tags` (a job parameter) — tags applied to every stage, winning over
  the templates and every derived value. A `null` value removes the tag, the
  template's line included, and keeps the runner from deriving one, so VASP
  uses its own default: `httk job new ... --parameter 'incar_tags={"NBANDS": null, "MAGMOM": null}'`
  runs with VASP's band count and moments.
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
* In a copied Bash runner, `httk_vasp_set_tag TAG VALUE [INCAR]` and
  `httk_vasp_unset_tag TAG [INCAR]` (VASP Bash API version 4) edit the INCAR
  after `httk_vasp_prepare`; unsetting an absent tag is no error. In a copied
  Python runner, `httk.codes.vasp.update_incar({"TAG": value}, "INCAR")` does
  the same, and a `None` value removes the tag.

*httk* v1 edited tags the same two ways: in the task template's INCAR, and with
its Bash INCAR helpers inside the runner script. Here the template route is
unchanged, the Bash helpers live on as `httk_vasp_set_tag`/`httk_vasp_unset_tag`,
and `incar_tags` adds a per-job route that needs no copied package.

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
   the accepted point is within `switch_force` / `switch_stress`, or, with
   `acceptance = gradient_tight`, after `stall_rejections` rejected steps in a
   row (the stall switch: `switch_reason` in the report is `limits` or `stall`).
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
coarse), `tight_cost_factor` times the loose one but at least the tight repeat
of the switch point (first tight: the repeat restarts at its own geometry and
is usually cheap, so it is timed apart and only a slow one raises the
estimate), or
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
every measured offset, the correction-model bound, the energy-gradient consistency fit, the symmetry counts, a
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
| `tight_cost_factor` | 1.35 | Expected cost of the first tight-EDIFF coarse step relative to the slowest loose evaluation, until one has run; the tight repeat of the switch point, timed apart since it restarts at its own geometry, is a lower bound (a cheap repeat, 149 s against 441 s for the next step in J-0210, cannot lower the estimate; a slow one raises it). Measured mean tight/loose warm cost: 1.15 (KBeF3 TET, J-0209), 1.34 (KSnF3 TET, J-0210). |
| `final_cost_factor` | 8 | Expected cost of the first final-mesh evaluation relative to the slowest coarse one. It is used only before any final-mesh evaluation has run. It is a rough estimate: mesh 5³ → 6³ cost 3.4× per SCF step, and Normal cost 2.3–2.8× Fast. |
| `maximum_evaluations`, `maximum_corrections` | 60, 1 | Cap on coarse evaluations (strict, including the tight repeat), and the number of offset corrections. |
| `maximum_interruptions` | 3 | Evaluations deferred because the manager stopped VASP at the attempt's deadline; one more fails the job. Every deferral is a new activation, which no other driver cap bounds. |
| `trust_atom`, `trust_strain`, `rise_tolerance` | 0.03 Å, 0.01, 2e-4 eV | Largest step per atom (physical displacement in the current cell) and per strain component (relative to the current cell) at trust factor 1. An evaluation whose energy rises more than `rise_tolerance` above the accepted point at the same EDIFF is rejected and the trust factor halves; which energy rise is set by `acceptance`. |
| `acceptance`, `rise_ceiling` | `gradient_tight`, 2e-3 eV | `energy`: the recorded energy rise. `gradient`: the rise the gradients predict along the step, ½ (g₀ + g₁)·s (exact on a quadratic), and a recorded rise above `rise_ceiling` still rejects. `gradient_tight`: `energy` at the loose EDIFF, `gradient` from the switch on. KSnF3 TET (J-0210, HSE, PRECFOCK Fast) had an energy that fell 2.15 kbar × ΔV faster with volume than its stress said (a fit; that it is a basis/Pulay term is a hypothesis, supported by restart jumps but not isolated); judged by the energy, 6 of 11 tight steps that lowered forces and stress were rejected and the job hit its cap at 0.0037 eV/Å. Gradient acceptance does not prove that forces and stress improve or that the gradients are right; the final static's own forces and stress decide. The ceiling is a per-step, per-cell heuristic: not size-extensive, and no bound on a drift over many smaller rises or on a downward jump. It also limits the cure: with a mismatch of σ kbar a step can change the volume by only about 1602 × `rise_ceiling` / σ Å³ (0.16 Å³ at 20 kbar, where the test model's cubic relaxation reaches a 24-evaluation cap; 1.5 Å³ at 2.15 kbar). Each step records `energy_rise`, `predicted_rise`, `energy_gradient_mismatch` and `volume_change`; the report (and a `vasp.relax_not_converged` failure, over the same evaluations it counts) gives per level the steps compared, the rms and largest mismatch, ΣΔV², and the isotropic stress that best explains the mismatch (`energy_gradient_consistency`). That stress is null unless at least two steps change the volume by 1e-4 of the cell (`identifiable`); it is a fit to the surface, not a measured Pulay stress, and the driver does not use it. |
| `stall_rejections` | 3 | With `gradient_tight`: this many rejected loose steps in a row tighten the EDIFF at the accepted point before the switch limits are reached, so the gradients judge from there (0: off). An energy that disagrees with the stress can reject the loose steps that would reach the limits: in the test model's cubic relaxation at 5 kbar the job otherwise stays loose until a 24-evaluation cap, and converges in 10 evaluations with the switch. It is one switch, and its tight repeat counts in the cap (Altair O1 review R2, 2026-10-18). |
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

## vasp-relax-ions

`vasp-relax-ions` relaxes the atomic positions of a structure in its fixed
cell, once and accurately, from `templates/INCAR.relaxions` (or a payload
`files/INCAR.relaxions`); the template's `EDIFF` and `EDIFFG` are kept. The
`isif` parameter (2, the default, or 4) sets `ISIF`: 4 also relaxes the cell
shape at fixed volume (an httk choice where atomate2's constant-volume
relaxations use 2). `ISIF`, `IBRION` and `NSW` in `incar_tags` fail
`vasp.unsupported_incar_tag` before any VASP run. The cell is never scaled by a
remedy (`scale_lattice` is excluded). A run that uses up its `NSW` ionic steps
is not accepted, unlike in `vasp-relax`'s pre-stages: the ladder continues it
from its CONTCAR (`contcar_to_poscar`), once, and a second such run fails
`vasp.failed`. The outputs are `relaxed_structure` (the CONTCAR) and
`total_energy`. It is the deformation leaf of the coming `vasp.elastic` and
`vasp.eos`, which call it once per strained cell with one `kpoint_grid` and one
`excluded_remedies` for all of them.

## vasp-dos

`vasp-dos` computes the electronic density of states and the band gaps of a
structure that is already chosen (normally relaxed), in its own cell: no
standardization, so `MAGMOM` and the other per-site tags keep their meaning.
The `scf` stage (`INCAR.scf`, a copy of `vasp-band-structure`'s) runs on the
automatic mesh of `kpoint_density` (default 20) and writes the CHGCAR; it is
archived under `scf/`. The `nscf` stage (`INCAR.nscf`) then reads it with
`ICHARG = 11` and `LCHARG = .FALSE.` on the denser automatic mesh of
`dos_kpoint_density` (default 30), with the tetrahedron method
(`ISMEAR = -5`), `LORBIT = 11` and `NEDOS = 3001`. Both meshes are
Γ-centred by default here (`centering = Gamma`, for the tetrahedron method).
`analyse` then reads the nscf EIGENVAL occupations and the DOSCAR Fermi energy
and writes `band_edges.json`. Remedies, the shared remedy history, the
`ICHARG = 11` enforcement and the density-tag guard are those of
`vasp-band-structure`; `scale_lattice` and `contcar_to_poscar` are never
applied, nor, in the nscf stage, an `ICHARG` change.

A given density skips the scf stage: a payload `files/CHGCAR` with
`files/scf-INCAR`, the INCAR it was computed with (and optionally
`files/scf-OUTCAR`, kept as `scf/OUTCAR` for provenance), as
`vasp.full-characterization` passes its band-structure child's. The nscf INCAR
is checked against `scf-INCAR` in the density-defining tags before the run.
No `total_energy` is collected then: that energy belongs to the job that
computed the density.

The band edges come from `httk.analyse.matsci.electronic.band_edges` with
occupations per state on 0..1, `occupation_tolerance` (default `1e-3`) and the
Fermi energy as the reference; a partially occupied state anywhere makes the
material metallic, and both gaps are then recorded as 0. The tetrahedron method's Blöchl correction gives states at the Fermi level of a metal occupations slightly outside 0..1 (say 1.0012 or −0.0007): any beyond `occupation_tolerance` counts as partial (0.5), and `band_edges.json` records how many as `occupations_out_of_range`. With `ISPIN = 2` (the
templates' default) the material's edges and indirect gap take both spin
channels together (minimum CBM over both minus maximum VBM over both; band
indices count the spin-up bands first), while the direct gap is the smallest
per-channel direct gap (spin-conserving transitions); `band_edges.json` also
holds each channel's own edges under `per_spin` (`null` for a channel without
both occupied and empty states).

The job fails by name with `vasp.unsupported_calculation` (`LHFCALC`,
`METAGGA`, `LNONCOLLINEAR`, `LSORBIT`, as in `vasp-band-structure`),
`vasp.invalid_dos_inputs` (`ICHARG` in `incar_tags`, or an nscf INCAR whose
`ICHARG` is not 11), `vasp.input_missing` (also `files/CHGCAR` without
`files/scf-INCAR`), `vasp.charge_density_missing`,
`vasp.charge_density_incompatible`, `vasp.eigenvalues_without_occupations`
(an EIGENVAL without occupation columns), `vasp.band_edges_failed` (unreadable
EIGENVAL or DOSCAR, or no band edges) and `vasp.analysis_unavailable`
(*httk-analyse* is not installed).

The outputs are `total_energy` (from `scf/OUTCAR`, only when the scf ran in
this job), the typed `electronic_density_of_states` (the total, spin-summed
DOS) and `fermi_energy` (from the DOSCAR), `band_gap` and `direct_band_gap`
(from `band_edges.json`, never re-analysed at collect), and three `files`
entries: `dos` (DOSCAR), `eigenvalues` (EIGENVAL) and `band_edges`
(`band_edges.json`). Collecting needs *httk-analyse*. With `publish_data` the
stages are published under `scf/` and `nscf/` below `data_prefix` (default
empty). Differences from atomate2's `NonSCFSetGenerator` uniform mode: the
nscf mesh is `dos_kpoint_density` (an httk choice; `kpoint_density = 20` is
already about atomate2's reciprocal density 100), `NEDOS = 3001` is an httk
choice, and no projected DOS is collected (`LORBIT = 11` writes it into the
DOSCAR and PROCAR, but only the total DOS is read).

## vasp-dielectric

One density-functional perturbation theory run (`IBRION = 8`, `LEPSILON`) of a
fixed structure: the static and high-frequency dielectric tensors, the Born
effective charges and the piezoelectric tensor. **Insulators only.** The
template `INCAR.dielectric` has the electronic settings of `INCAR.static` plus
atomate2's `StaticSetGenerator(lepsilon=True)` choices (`IBRION = 8`,
`LEPSILON`, `LPEAD`, `NSW = 1`, `ISMEAR = -5`, `LREAL = .FALSE.`) and the *httk*
choices `EDIFF = 1E-7`, atomic starting charges (`ICHARG = 2`), and no WAVECAR or CHGCAR.

`incar_tags` may not set `IBRION`, `LEPSILON`, `NCORE` or `NPAR`
(`vasp.unsupported_incar_tag`), and `parallel_tag` may not be `NPAR` or `NCORE`
(`vasp.unsupported_parallel_tag`): VASP forbids both with `IBRION = 8`; use `KPAR`.
A Cartesian POSCAR is refused (`vasp.input_unsupported`). Parameters: `kpoint_grid = "n1 n2 n3"` fixes the
mesh, `excluded_remedies` names ladder operations never to apply, `centering`
defaults to `Gamma` (the tetrahedron method wants it); the general parameters apply.
The structure is never changed: `scale_lattice`, `contcar_to_poscar` and `incar.NPAR` are never applied. `rattle_amplitude` is not used, because the response refers to the structure as given.

A completed run's OUTCAR is read into `dielectric.json` before publishing; a
run without a dielectric tensor, or with a Born charge count that is not the atom
count, fails `vasp.dielectric_missing`. The file holds `electronic` and `ionic`
(3x3; `static` is their sum), `born_charges` (`[site][field direction][displacement direction]`, rows as VASP prints
them, in the POSCAR's site order), `piezoelectric_electronic` and `piezoelectric_ionic`
(3x6 in VASP's column order, named by `piezoelectric_column_order`: XX YY ZZ XY YZ ZX),
`lattice` (3x3, Å, scale applied) and `positions` (direct coordinates) of the POSCAR VASP ran,
`species`, `counts`, `poscar_sha256` and `energy`. Absent blocks are `null`.

Outputs: `dielectric_cell`, `static_relative_permittivity` and
`high_frequency_relative_permittivity` (isotropic means), `born_effective_charges`,
`piezoelectric_tensor` (electronic plus ionic, reordered to Voigt xx yy zz yz xz xy),
`total_energy`, and `dielectric_response` (`dielectric.json`). Collecting needs *httk-analyse*.

## vasp-elastic

`vasp-elastic` computes the elastic stiffness tensor of one structure and runs
no VASP itself. It declares `[workflow.calls] relax = "vasp.relax"` and
`relax_ions = "vasp.relax-ions"` and runs them as child jobs:

- `start` validates the parameters (anything invalid fails
  `vasp.input_invalid` before any child is called). With `relax` (default
  true) it calls `vasp.relax` (child `relax`) with this package's
  `templates/INCAR.relax` for its relax rounds: `vasp-relax`'s template with
  atomate2 `TightRelaxSetGenerator`'s `EDIFF = 1E-7`, `EDIFFG = -0.001` and
  `LREAL = .FALSE.` (a payload `files/INCAR.relax` replaces it); the pre-stages
  stay loose. A failed relaxation fails `vasp.relax_failed`. With
  `relax = false` the input structure is strained as given (for a structure
  that is already relaxed tightly).
- `deform` copies the relaxed structure to `relaxed/POSCAR` and writes one
  strained cell per strain state `s` and magnitude `m`,
  `deformations/d00/POSCAR`, `d01`, … (state-major): the Voigt strain
  `e = m s` (`xx yy zz yz xz xy`, engineering shear strains, as httk-analyse's
  fit takes them), the symmetric small strain `eps`, and the lattice rows
  `a_i' = (I + eps) a_i`. Only the lattice rows are rewritten, so the
  structure must be a VASP-5 POSCAR with scale 1 and `Direct` coordinates,
  without `Selective dynamics` (otherwise `vasp.input_unsupported`). It calls
  `vasp.relax-ions` with `isif = 2` once per cell and waits for all of them.
- `fit` reads every child's stress (OUTCAR, tensile-positive GPa) and fits
  the 21 stiffness components and a stress offset by least squares
  (httk-analyse `fit_stress_strain(..., fit_offset=True)`).

Parameters: `relax` (default true); `strain_states` (default the six unit
states; a list of six-number lists, or a JSON string of one; together they
must span all six components); `strain_magnitudes` (default
`-0.01 -0.005 0.005 0.01`, atomate2's; at least two distinct nonzero
numbers); `deformation_kpoint_density` (default 40, an httk choice denser
than the relaxation's 20, as atomate2's `ElasticRelaxMaker` uses a denser
mesh); `centering` (default `Gamma`); `publish_data` (default false: also put
`elastic.json` and `relaxed/POSCAR` under `data/<data_prefix>/`). The shared
parameters `incar_tags`, `accuracy_per_atom`, `maximum_remedies`,
`pseudopotential_library`, `parallel_tag`, `parallel_value`, `timeout`,
`vasp_command` and `remedy_policy` are passed to every child; the POTCAR
travels as a file.

Every deformation runs with one numerical setting: the mesh is computed once
on `relaxed/POSCAR` and passed as `kpoint_grid` with `centering`, and the
remedies that would change the numerics of one cell but not the others are
excluded (`excluded_remedies = "bump_kpoints equal_kpoints centering
incar.ISMEAR incar.SIGMA incar.LREAL scale_lattice"`). Before fitting, `fit`
compares the children's KPOINTS and INCARs and fails
`vasp.inconsistent_children` (the differing tags per label in the details)
when the KPOINTS or any tag other than the electronic-path tags (`ALGO TIME
AMIX BMIX AMIN ICHARG ISTART NELM NELMIN NBANDS NPAR NCORE KPAR ISYM
SYMPREC EDIFF NEDOS`) differ. Other failures: `vasp.dependency_failed`
(any deformation did not succeed; every failed label with its code in the
details; no fit, rerun with `httk job continue`), `vasp.child_workdir_missing`,
`vasp.stress_missing`, `vasp.analysis_unavailable` (no httk-analyse where the
job runs) and `vasp.elastic_fit_failed`.

`elastic.json` records the reference lattice, the mesh, each deformation's
strain, and per child its label, job key and id, outcome, INCAR tags, KPOINTS,
applied remedy operations and stress; then the fit: `stiffness` (6×6, GPa),
`stress_offset`, `residuals`, `rmse`, `condition_number`, `strain_range`,
`stable`, `stability_eigenvalues`, `compliance` (null when singular) and, for a
stable tensor, the six Voigt/Reuss/Hill moduli and `universal_anisotropy`.
`collect.py` reads only it and `relaxed/POSCAR` (it needs httk-analyse) and
rebuilds the records from the stored stiffness, never re-fitting. The outputs
are `relaxed_structure`, `elastic_tensor`, `compliance_tensor`,
`bulk_modulus_voigt`/`_reuss`/`_hill`, `shear_modulus_voigt`/`_reuss`/`_hill`,
`universal_anisotropy_index` (the roles and refs of `atomate2.elastic`) and
`elastic_result` (`elastic.json`). An unstable tensor still succeeds, says so
in the log and in `stable`, and collects only `elastic_tensor` and
`compliance_tensor`; the moduli roles are unfulfilled.

Differences from atomate2's `ElasticMaker`: the deformation set is not reduced
by symmetry (the full set is correct, only more expensive; a follow-up on
httk-atomistic's symmetry); `F = I + eps` is symmetric where atomate2's is
upper-triangular, an O(eps²) difference the offset fit absorbs; the tight
relaxation keeps `ENCUT = 600` where atomate2 sets 700; there is no
`max_failed_deformations` (any failed deformation fails the job). The
deformations run as concurrently as the manager has slots; a child's `maxtime`
is capped at the parent attempt's, so this job declares none.

## vasp-eos

`vasp-eos` computes the third-order Birch–Murnaghan equation of state of one
structure and runs no VASP itself. It declares `[workflow.calls] relax =
"vasp.relax"`, `relax_ions = "vasp.relax-ions"` and `static = "vasp.static"`
and runs them as child jobs:

- `start` validates the parameters (anything invalid fails
  `vasp.input_invalid` before any child is called). With `relax` (default
  true) it calls `vasp.relax` (child `relax`) with this package's tight
  `templates/INCAR.relax`, as `vasp-elastic` does (a payload
  `files/INCAR.relax` replaces it; the pre-stages stay loose); a failed
  relaxation fails `vasp.relax_failed`. With `relax = false` the input
  structure is the reference as given.
- `deform` copies the reference to `relaxed/POSCAR` and, for each of
  `number_of_frames` linear strains `ε` evenly spread over `linear_strain`,
  scales every lattice row by `1 + ε` (volume `V₀(1 + ε)³`) into
  `volumes/v00/POSCAR`, `v01`, … Only the lattice rows are rewritten, so the
  structure must be a VASP-5 POSCAR with scale 1 and `Direct` coordinates,
  without `Selective dynamics` (otherwise `vasp.input_unsupported`). It calls
  `vasp.relax-ions` with `isif = 4` once per volume (the atoms and the cell
  shape relax at fixed volume) and waits for all of them.
- `statics` (only with `static = true`) calls `vasp.static` (children `s00`,
  `s01`, …) on each relaxed CONTCAR.
- `fit` takes each volume from its `vasp.relax-ions` CONTCAR and each energy
  (`energy(sigma->0)`) from the OUTCAR of the static, or with `static = false`
  of the relaxation itself, and fits httk-analyse's `fit_birch_murnaghan`.

Parameters: `relax` (default true); `linear_strain` (default `[-0.05, 0.05]`,
two numbers or a space-separated string, below zero and above it, within
±0.3, so the curve has an interior minimum); `number_of_frames` (default 6, at
least 5); `static` (default false); `kpoint_density` (default 20) and
`centering` (default `Gamma`) of the one mesh; `publish_data` (default false:
also put `eos.json` and `relaxed/POSCAR` under `data/<data_prefix>/`). The
shared parameters `incar_tags`, `accuracy_per_atom`, `maximum_remedies`,
`pseudopotential_library`, `parallel_tag`, `parallel_value`, `timeout`,
`vasp_command` and `remedy_policy` are passed to every child; the POTCAR
travels as a file.

Every volume and static runs with one numerical setting, by the rule of
`vasp-elastic`: one mesh computed on `relaxed/POSCAR`, passed as
`kpoint_grid` with `centering`, and the same `excluded_remedies`. `fit`
compares the KPOINTS and INCARs of the children that provide the energies and
fails `vasp.inconsistent_children` when the KPOINTS or any tag differ other
than the electronic-path tags `ALGO TIME AMIX BMIX AMIN ICHARG ISTART NELM
NELMIN NBANDS NPAR NCORE KPAR ISYM SYMPREC EDIFF NEDOS` (a remedy may tighten
`EDIFF` or raise `NEDOS` in one child without changing its energy). Other
failures:
`vasp.dependency_failed` (any child did not succeed; every failed label with
its code in the details; rerun with `httk job continue`),
`vasp.child_workdir_missing`, `vasp.energy_missing`,
`vasp.analysis_unavailable` and `vasp.eos_fit_failed` (for instance no
interior minimum; the data stay in `eos.json`).

`eos.json` records the plan (`strains`, `reference_volume`, `kpoint_grid`,
`centering`, `excluded_remedies`), per child its label, job key and id,
outcome, INCAR tags, KPOINTS and applied remedy operations, per volume its
strain, volume and energy, `energy_source` (`relax_ions` or `static`), and
every field of the fit: `equilibrium_volume` (Å³), `equilibrium_energy` (eV),
`bulk_modulus` (GPa), `bulk_modulus_derivative`, `volumes`, `energies`,
`residuals`, `rmse` and `condition_number`. `collect.py` reads only it and
`relaxed/POSCAR` (it needs httk-analyse) and rebuilds the fit from the stored
values, never re-fitting. The outputs are `relaxed_structure`,
`equilibrium_volume`, `equilibrium_energy`, `bulk_modulus`,
`bulk_modulus_pressure_derivative` and `eos_result` (`eos.json`); the energy
RMSE is only in the JSON.

Differences from atomate2's `EosMaker`: the fixed-volume relaxations use
`ISIF = 4` (an httk choice; atomate2 uses 2); `static` defaults to false, as
in the generic `EosMaker` (`static_maker = None`; the MP variants set one);
the tight relaxation keeps `ENCUT = 600` where atomate2 sets 700; the fit is
always Birch–Murnaghan. The volumes run as concurrently as the manager has
slots; a child's `maxtime` is capped at the parent attempt's, so this job
declares none.

## vasp-phonons

`vasp-phonons` computes the harmonic phonons of one structure with phonopy and
runs no VASP itself. It declares `[workflow.calls] relax = "vasp.relax"` and
`static = "vasp.static"` and runs them as child jobs:

- `start` validates the parameters (anything invalid fails
  `vasp.input_invalid` before any child is called), refuses `NSW` and `IBRION`
  in `incar_tags` (`vasp.unsupported_incar_tag`) and the cell-dependent
  `MAGMOM`, `NELECT`, `NUPDOWN`, `NG[XYZ]`, `NG[XYZ]F` or a POSCAR comment's
  `[MAGMOM=...]` (`vasp.cell_dependent_tag_unsupported`: magnetic phonons are
  out of scope for now). With `relax` (default true) it calls the tight
  `vasp.relax` exactly as `vasp-elastic` does; with `relax = false` the input
  structure is taken as relaxed.
- `displace` copies the relaxed structure to `relaxed/POSCAR`. With
  `standardize` (default true) the phonon unit cell `unitcell/POSCAR` is its
  HPKOT primitive cell from `httk.atomistic.band_path`, in the input's species
  order (`vasp.input_mismatch` otherwise), as in `vasp-band-structure`; with
  `standardize = false` the given cell is used as it is (so a `dielectric.json`
  of that cell still matches), and `band_path`'s cell `L_path` must be another
  basis of the same lattice in the same frame: `T = L_path L⁻¹` (rows as
  vectors) integer within 1e-6 with `|det T| = 1`, else `vasp.input_invalid`
  (pass `standardize = true`). The unit cell's handedness is normalized once,
  and `unitcell/band_path.json` holds the q-point path in that very cell's
  reciprocal basis, `q = q_path T⁻ᵀ`, with `T` as `basis_change` (the
  identity, or its negation after the handedness normalization, when
  standardized). The supercell is diagonal, `n_i = ceil(min_length / w_i)`
  with `w_i` the cell's width perpendicular to the other two vectors; more
  than `max_atoms` atoms fail `vasp.supercell_too_large` (the multiples and the
  atoms in the details). phonopy (`primitive_matrix = "P"`: the unit cell is
  the cell diagonalized) writes `phonopy_disp.yaml` and the displaced
  supercells `displacements/000/POSCAR`, … (`vasp.input_mismatch` if one does
  not list the unit cell's species line). One mesh is computed on the
  undisplaced `supercell/POSCAR`, and `vasp.static` runs once per displacement
  (children `p000`, `p001`, …) with `templates/INCAR.phonon` as its
  `INCAR.static`, that mesh as `kpoint_grid` with `centering`, and the
  `excluded_remedies` of `vasp-elastic`.
- `analyse` reads every child's forces (the last OUTCAR frame), sets them on
  phonopy's displacement dataset from `phonopy_disp.yaml`, writes `FORCE_SETS`
  and produces the force constants, then runs the mesh, the total DOS, the
  band structure and httk-analyse's `harmonic_from_phonopy` thermodynamics.

The non-analytical correction: a payload `files/dielectric.json` (as
`vasp-dielectric` writes it) is checked in `displace`, before any child: it
must have Born charges and `electronic`, the unit cell's `species` and
`counts`, and its `lattice` and `positions` (modulo 1) within 1e-5 Å and 1e-5
of the unit cell's, else `vasp.dielectric_mismatch`. It is then copied to
`unitcell/dielectric.json`, and phonopy gets `nac_params = {"born",
"dielectric": electronic, "factor": 14.39965}` (Born rows as VASP prints them,
the field direction, as phonopy reads them; the factor is e²/(4πε₀) in eV Å,
which phonopy 4.8 requires). Since Born charges and ε∞ must come from the
same cell and site order, a dielectric run on the cell `vasp-phonons`
standardizes to is the reliable source; `vasp.full-characterization` passes
its own HPKOT cell with `standardize = false`.

Parameters: `relax` (default true); `standardize` (default true);
`min_length` (default 20.0 Å, atomate2's); `max_atoms` (default 300);
`displacement` (default 0.01 Å); `phonon_kpoint_density` (default 20) and
`centering` (default `Gamma`) of the one supercell mesh; `mesh_length`
(default 100.0, phonopy's length-based q mesh); `band_points` (default 51
q-points per path segment); `temperatures` (default 0 to 1000 K in steps of
10, a list or a space-separated string); `cutoff_frequency` (default 0.05 THz:
modes below it in magnitude are omitted as zero modes); `imaginary` (`raise`,
the default, or `omit`); `imaginary_tolerance` (default 0.05 THz: frequencies
below its negative are reported as imaginary); `symmetry_tolerance` (default
1e-3 Å); `publish_data` (default false: also put `phonons.json`,
`relaxed/POSCAR`, `unitcell/POSCAR`, `phonopy_disp.yaml` and `FORCE_SETS` under
`data/<data_prefix>/`). The shared parameters are passed to every child as in
`vasp-elastic`; the POTCAR travels as a file.

Units: frequencies in THz; the DOS in modes per THz per unit cell (it
integrates to three times the unit cell's sites); the thermodynamics in eV and
eV/K per unit cell, the zero-point energy in eV. The band structure's q-points
are fractional coordinates of the unit cell's reciprocal basis, `band_points`
per segment with the HPKOT labels at both ends of each segment, and its
frequencies are sorted at each q-point (`frequencies[mode][q]`, as
`atomate2.phonons`' collect stores them).

Imaginary modes: the count of mesh frequencies below `-imaginary_tolerance`
and the lowest frequency are always in `phonons.json` and the log. With
`imaginary = "raise"` any negative frequency beyond the zero-mode cutoff
leaves the thermodynamics out (`null`, with the reason, in `phonons.json`;
the job still succeeds, and `vibrational_thermodynamics` and
`zero_point_energy` are unfulfilled); `omit` computes them without those
modes.

`phonons.json` records the plan (the validated parameters, the unit cell, the
supercell matrix and atoms, each displacement's atom and vector, the mesh,
`excluded_remedies`), per child its label, job key and id, outcome, INCAR
tags, KPOINTS, applied remedies and forces, and the results: `nac`, `mesh`,
`imaginary_modes`, every field of httk-analyse's `HarmonicThermodynamics`,
`dos` and `band_structure`. `collect.py` reads only it, `unitcell/POSCAR` and
`relaxed/POSCAR` (it needs httk-analyse) and never recomputes. The outputs
are `relaxed_structure` (only with `relax`), `phonon_cell`,
`phonon_band_structure`, `phonon_density_of_states`,
`vibrational_thermodynamics` and `zero_point_energy` (the roles and refs of
`atomate2.phonons`), and the files `phonons_result` (`phonons.json`),
`phonopy_disp` and `force_sets`. Failures besides the above:
`vasp.relax_failed`, `vasp.band_path_failed`, `vasp.dependency_failed`,
`vasp.child_workdir_missing`, `vasp.inconsistent_children` (as in
`vasp-elastic`, but `EDIFF` and `NEDOS`, which the `scale_ediff` and
`incar.NEDOS` remedies change, may also differ), `vasp.forces_missing` (no forces on every supercell atom)
and `vasp.analysis_unavailable` (no phonopy or httk-analyse where the job
runs).

Differences from atomate2's `PhononMaker`: the supercell is diagonal by
perpendicular width (atomate2 searches cubic-like non-diagonal supercells, a
follow-up here); the unit cell is always the HPKOT cell and the path
`httk.atomistic.band_path`'s, not phonopy's seekpath; Born charges come from a
given `dielectric.json` rather than a dielectric run of its own; the
displacement runs keep `ENCUT = 600` where atomate2 sets 700, and add
`ISMEAR = 0`, `SIGMA = 0.05` and `ADDGRID = .TRUE.`; no static energy run. The
displacements run as concurrently as the manager has slots; a child's
`maxtime` is capped at the parent attempt's, so this job declares none.

## vasp-full-characterization

`vasp-full-characterization` characterizes one material with the other
packages of this repository as child jobs, in parallel where their
dependencies allow, and runs no VASP itself. It declares `[workflow.calls]
relax = "vasp.relax"`, `band_structure = "vasp.band-structure"`, `dos =
"vasp.dos"`, `elastic = "vasp.elastic"`, `eos = "vasp.eos"`, `dielectric =
"vasp.dielectric"` and `phonons = "vasp.phonons"`:

```
start ── call relax ──► plan ── standardize, decide insulator or metal, call band_structure,
         (tight)                elastic, eos, and dielectric (insulator) or phonons (metal) ──► progress
progress ── whenever a child ends: once band_structure has ended, call dos on its density;
            once dielectric has ended, call phonons with its Born charges; when none is live, finish
```

- `start` validates the parameters (`vasp.input_invalid`) and refuses the
  cell-dependent `MAGMOM`, `NELECT`, `NUPDOWN`, `NG[XYZ]`, `NG[XYZ]F` in
  `incar_tags` or any `<label>_parameters.incar_tags`, or a POSCAR comment's
  `[MAGMOM=...]` (`vasp.cell_dependent_tag_unsupported`), before any child.
  It calls `vasp.relax` (child `relax`) with this package's tight
  `templates/INCAR.relax`, as `vasp-elastic` does; a failed relaxation fails
  `vasp.relax_failed`.
- `plan` copies the relaxed CONTCAR to `relaxed/POSCAR` and, with
  `standardize` (default true), writes its HPKOT primitive cell to
  `standardized/POSCAR` exactly as `vasp-band-structure` builds it (within
  `symmetry_tolerance`, default 1e-3 Å; `vasp.band_path_failed`,
  `vasp.input_mismatch`); with `standardize = false` the relaxed structure is
  copied (phonons then needs it to be its own HPKOT cell). Every property child
  starts from that cell, with the POTCAR as a file.
- `progress` keeps the ended children and rejoins the live ones until none is
  live, then writes `characterization.json` and finishes.

No calculation runs twice. The tight relaxation is the only one: `elastic`,
`eos` and `phonons` get `relax = false` (and phonons `standardize = false`).
`dos` gets the band-structure child's HPKOT `POSCAR`, `POTCAR` and `CHGCAR`
with its `scf/INCAR` and `scf/OUTCAR` as `files/scf-INCAR` and
`files/scf-OUTCAR`, so it runs only its nscf stage and the band-structure scf
is the only self-consistent run on the cell (if band structure fails, dos runs
its own scf; a note says so). `phonons` gets the dielectric child's
`dielectric.json` for the non-analytical correction, computed on the very cell
phonons runs on (without it, if dielectric fails: a note says so, and the job
cannot succeed, see below). `characterization.json` records `phonons_nac` among its decisions.

The insulator decision is taken from the relaxation, as atomate2's parent
takes it: the band edges (httk-analyse's `band_edges`, both spin channels
stacked) of the occupations in the final relax round's `EIGENVAL`. That round
runs with the tetrahedron method (`ISMEAR = -5`), whose Blöchl-corrected
occupations at the Fermi level of a metal can lie slightly outside [0, 1]: an
occupation outside [−1e-3, 1 + 1e-3] counts as 0.5, a partial one
(`occupations_out_of_range` records how many). An insulator has no partial occupation and an indirect gap above
`band_gap_tolerance` (default 0.1 eV, an httk choice; atomate2 uses its
`BANDGAP_TOL`). Only an insulator gets a dielectric child. No readable
EIGENVAL with occupations fails `vasp.band_gap_unknown`: nothing is guessed.
The dos child's gap is recorded too, with a note when it disagrees.

Parameters: `parts` (default `"band_structure dos elastic eos dielectric
phonons"`, any subset, space-separated or a list; an unknown name is
`vasp.input_invalid`; dos without band_structure runs its own scf);
`standardize`; `symmetry_tolerance`; `band_gap_tolerance`;
`<label>_parameters` (`relax_parameters`, `band_structure_parameters`,
`dos_parameters`, `elastic_parameters`, `eos_parameters`,
`dielectric_parameters`, `phonons_parameters`), dictionaries of job parameters
for that child, applied over the shared parameters (a member replaces the
shared one whole), with this job's own structural settings winning over them;
`poscar` and `potcar`. The shared parameters `incar_tags`, `accuracy_per_atom`,
`maximum_remedies`, `pseudopotential_library`, `parallel_tag`,
`parallel_value`, `timeout`, `vasp_command` and `remedy_policy` reach every
child. The fan-out children compute their own meshes; dos reads the density
only when both children ran with the same density-defining tags.

The job yields `relaxed_structure`, `standardized_structure` and
`characterization` (`characterization.json`: each child's job id, key, outcome
and failure, and the decisions and notes); all three are also in its data.
Every property is collected on the child runs, each under its own package and
declaration, with parent → child run edges: band structure
(`band_structure_cell`, `total_energy`, `band_eigenvalues`, `band_path`), dos
(`band_gap`, `direct_band_gap`, `fermi_energy`,
`electronic_density_of_states`, …), elastic (`elastic_tensor` and the moduli),
eos (`bulk_modulus`, `equilibrium_volume`, …), dielectric
(`static_relative_permittivity`, `high_frequency_relative_permittivity`,
`born_effective_charges`, …), phonons (`phonon_band_structure`,
`zero_point_energy`, …) and the relaxation (`relaxed_structure`,
`total_energy`). The job succeeds when every child it called succeeded, and
otherwise fails `vasp.characterization_incomplete` naming the failed ones
(their failures in the details); succeeded children are collected either way.
An insulator whose phonons ran without `dielectric.json` (dielectric in
`parts` but failed) also fails `vasp.characterization_incomplete`, listing
`phonons (without the non-analytical correction: dielectric did not succeed)`,
even after `httk job continue` repairs the dielectric. The remedy: run
`vasp.phonons` with `relax = false` and `standardize = false` on the job's
`standardized/POSCAR`, with the dielectric child's `dielectric.json` as
`files/dielectric.json`.
Magnetic materials are out of scope for now. The job declares no `maxtime`,
which would cap every child's.

The easiest launch on a daemon-run workspace installs every package of a tag,
with calls between them bound to the same commit:

```console
httk workflow install --workspace main 'git+https://github.com/httk/workflows-vasp@<tag>'
httk job new --workspace main --workflow vasp.full-characterization --input-from structure Si.cif
```

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
  cell, and the HPKOT cell's sites and size differ from the input's. A `null`
  in `incar_tags` only removes the tag (VASP's default) and is allowed.
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
* `kpoint_grid` (default empty, `vasp-static` and `vasp-relax-ions` only) — a
  fixed automatic mesh `"n1 n2 n3"` (or a list of three integers), written with
  the job's `centering` instead of the mesh of `kpoint_density`, so every child
  of a fan-out runs on one mesh. Anything but three positive integers fails
  `vasp.input_invalid`.
* `excluded_remedies` (default empty, `vasp-static` and `vasp-relax-ions` only)
  — ladder operations never to apply, space-separated (for example
  `"bump_kpoints equal_kpoints centering incar.ISMEAR incar.LREAL"`); the ladder
  goes on to its next rung instead. Each must be an operation name
  (`bump_bands bump_kpoints centering contcar_to_poscar equal_kpoints scale_ediff
  scale_lattice` or `incar.<TAG>`), or the job fails `vasp.input_invalid`; names
  are not checked against the ladder, so one it never uses excludes nothing.
* `isif` (default `2`, `vasp-relax-ions` only) — `2` relaxes the positions in
  the fixed cell, `4` also the cell shape at fixed volume; anything else fails
  `vasp.input_invalid`.
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
  attempt, so two retries never repeat one structure. `vasp-band-structure`,
  `vasp-dos` and `vasp-dielectric` do not use it, and the composites
  (`vasp-elastic`, `vasp-eos`, `vasp-phonons`, `vasp-full-characterization`) do not
  forward it to their children: their structure is the one their results refer to.
* `band_divisions` (default `40`, `vasp-band-structure` only) — k-points per
  segment of the line-mode band path.
* `symmetry_tolerance` (default `1e-3`, in Å, `vasp-band-structure` only) —
  the tolerance the structure's symmetry is recognized within, also the
  precision its POSCAR is read with; recorded in `band_path.json`.
* `dos_kpoint_density` (default `30.0`, `vasp-dos` only) — the k-point
  density of the nscf stage's uniform mesh; `kpoint_density` is the scf
  stage's.
* `occupation_tolerance` (default `1e-3`, `vasp-dos` only) — occupations
  within this of 0 or 1 count as empty or full; any other is partial, making
  the material metallic.
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
