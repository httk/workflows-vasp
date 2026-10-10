# workflows-vasp

This repository is the home of the *httk₂* VASP workflows, one self-contained
workflow package directory per workflow, referenced directly by a git commit.
They use the `httk.codes.vasp` helper library of an installed
*httk-workflow-vasp* (`pip install httk-workflow-vasp`).

| Directory | Short name | What it does |
| --- | --- | --- |
| `vasp-relax` | `vasp.relax` | Relax the geometry of a structure with the reviewed remedy ladder (Python runner). |
| `vasp-relax-bash` | `vasp.relax-bash` | The same relaxation, authored in Bash. |
| `vasp-static` | `vasp.static` | One single-point total-energy calculation of a fixed structure. |
| `vasp-relax-static` | `vasp.relax-static` | Relax, promote the relaxed structure, then run it statically. |

Each runner entry (`run.py`, or `run.sh` in `vasp-relax-bash`) is a starting
point to copy and edit: its steps spell out their work — reading the job
parameters, staging inputs, building the preparation options, running VASP
under supervision, planning and applying one remedy of the reviewed ladder, and
publishing — directly on the primitives of `httk.codes.vasp`. The
three Python runners (`vasp-relax`, `vasp-static`, `vasp-relax-static`) are
independent of each other and deliberately repeat that code; `vasp-relax-bash`
is the same relaxation, step for step, in Bash against the sibling Bash VASP
API.

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

## Job parameters

Every packaged VASP runner reads the same `inputs` object, and every member is
optional. Paths are relative to the job payload; lists are space-separated
strings, so the Bash runner and the Python runners read one contract.

* `poscar` (default `files/POSCAR`) — the starting structure; any VASP-5
  POSCAR or CONTCAR file will do.
* `incar` (default `files/INCAR`) — the starting INCAR. When the file is
  absent the runner starts from an empty INCAR and derives everything.
* `potcar` (default `files/POTCAR`) — a pre-assembled POTCAR. When the file
  is absent and `pseudopotential_library` is set, the POTCAR is assembled per
  species and a provenance record is written next to it.
* `pseudopotential_library` (default none) — root of a VASP pseudopotential
  library, one directory per variant.
* `kpoint_density` (default `20.0`), `centering` (default `Monkhorst-Pack`),
  `accuracy_per_atom` (default `0.001`) — passed to
  `httk.codes.vasp.inputs.VaspPreparationOptions`.
* `parallel_tag` and `parallel_value` (default none) — one of `NPAR`, `NCORE`,
  or `KPAR`, and its value.
* `incar_tags` (default empty) — explicit INCAR tags. They are applied before
  anything is derived and win over every derived value.
* `static_incar_tags` (default `IBRION = -1` and `NSW = 0`) — only in
  `vasp-static` and `vasp-relax-static`: the tags that turn the calculation
  into a single point.
* `timeout` (default `86400`) — seconds one VASP execution may take before its
  process group is terminated.
* `maximum_remedies` (default `8`) — how many remedies this job may apply in
  total before it fails. The ladder is bounded per problem as well.
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
  Missing files are skipped.
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

installs all four workflow packages listed in `httk_plugin.toml` at once on
this machine; each is still installed into a workspace (`httk job new
--install --workflow vasp.relax`) before its jobs are created.
