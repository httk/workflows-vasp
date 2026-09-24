# workflows-vasp

This repository is the home of the *httk₂* VASP workflows, one self-contained
workflow package directory per workflow, referenced directly by a git commit.
They use the `httk.workflow.vasp` helper library of an installed
*httk-workflow*.

| Directory | Short name | What it does |
| --- | --- | --- |
| `vasp-relax` | `vasp.relax` | Relax the geometry of a structure with the reviewed remedy ladder (Python runner). |
| `vasp-relax-bash` | `vasp.relax-bash` | The same relaxation, authored in Bash. |
| `vasp-static` | `vasp.static` | One single-point total-energy calculation of a fixed structure. |
| `vasp-relax-static` | `vasp.relax-static` | Relax, promote the relaxed structure, then run it statically. |

Three of the four `run` files (`vasp-relax`, `vasp-static`, `vasp-relax-static`)
are thin step declarations built on the Python VASP step API of
`httk.workflow.vasp` (`httk.workflow.vasp.steps`): each step delegates to one
library function — staging inputs, running VASP under the reviewed remedy
ladder, promoting a relaxation to a single point, and publishing a stage — and
only reads its own job parameters. `vasp-relax-bash` is the same relaxation
authored in Bash against the sibling Bash VASP API. Each `collect.py` calls the
matching `httk.workflow.vasp.collect` function. Running a job still requires an
installed *httk-workflow*, since the runners import `httk.workflow.vasp` for
inputs, remedies, diagnostics, steps, and collection.

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
  `httk.workflow.vasp.inputs.VaspPreparationOptions`.
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
  policy with `httk.workflow.vasp.register_remedy_policy` and names it here
  instead of editing a runner.
* `rattle_amplitude` (default `0.0`) — when positive, the POSCAR is rattled by
  this amplitude after every applied remedy, with a seed derived from the
  attempt, so two retries never repeat one structure.
* `collect` (default `INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml
  vasp-run-report.json POTCAR.provenance.json`) — space-separated file names
  copied to transactional data only when opted in. The default `data.mode`
  `none` keeps outputs in the workdir. `vasp-relax-static` also uses this list
  to archive the relaxation before the static stage. Missing files are skipped.
* `data_prefix` (default `vasp`) — directory below the job's data the
  collected files are published under when transactional data is enabled;
  ignored for workdir results. `vasp-relax-static` defaults to an empty prefix
  and publishes its stages under `relax/` and `static/`.
* `vasp_command` (default empty) — the VASP command as one argv string, split
  the way a shell splits it. The environment variable `HTTK_VASP_COMMAND`
  overrides it, which is how a deployment — or a test — chooses the executable
  without touching any job.

A job running one of these workflows needs `workdir.mode` `persistent`: the
inputs a remedy rewrites have to be the inputs the next attempt reads. All four
workflows default to `data.mode` `none`: the persistent workdir is the result,
with no `data/` copy. Collectors read it directly. Pass `--data-mode
transactional` to `httk job new` (or `data_mode="transactional"` to
`new_job`) to also publish the curated files into `data/`.

## Reference a workflow by commit

```console
httk job new --workflow 'git+https://github.com/httk/workflows-vasp@<ref>#vasp-relax' --input structure=POSCAR
```

`<ref>` may be a commit, branch, or tag, or omitted for the default branch; it
is canonicalized to the full commit hash the repository is cloned at, and the
named subdirectory's `httk_workflow.toml` package is used. Once a workflow has
been referenced this way, its short name (e.g. `vasp.relax`) also resolves.

## Install as a plugin

```console
httk plugin install 'git+https://github.com/httk/workflows-vasp'
```

installs all four workflow packages listed in `httk_plugin.toml` at once.
