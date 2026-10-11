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
| `vasp-general` | `vasp.general` | Any VASP calculation, run on input files the job supplies; only the POTCAR may be assembled. |

Each package directory — its runner entry (`run.py`, or `run.sh` in
`vasp-relax-bash`) with its `templates/` — is a starting point to copy and
edit as a whole: the runner's steps spell out their work — reading the job
parameters, staging inputs, building the preparation options, running VASP
under supervision, planning and applying one remedy of the reviewed ladder, and
publishing — directly on the primitives of `httk.codes.vasp`. The
Python runners (`vasp-relax`, `vasp-static`, `vasp-relax-static`, `vasp-general`) are
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
  density of the relax rounds and the static stage.
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
  attempt, so two retries never repeat one structure.
* `publish_data` (default `false`) — when `true`, `publish` also copies the
  `collect` files into the job's `data/` directory; by default outputs stay in
  the workdir only.
* `collect` (default `INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml
  vasp-run-report.json POTCAR.provenance.json`) — space-separated file names
  copied to `data/` only when `publish_data` is `true`. `vasp-relax-static`
  also uses this list to archive the relaxation before the static stage.
  Missing files are skipped. The `stages/` archives are not published.
* `data_prefix` (default `vasp`) — directory below the job's data the
  collected files are published under when `publish_data` is `true`;
  ignored for workdir results. `vasp-relax-static` defaults to an empty prefix
  and publishes its stages under `relax/` and `static/`.
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
