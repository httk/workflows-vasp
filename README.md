# workflows-vasp

Packaged copies of the four built-in VASP workflows of *httk-workflow*, one
self-contained workflow package directory per workflow, so each can be
referenced directly by a git commit instead of only by the *httk-workflow*
release that bundles it.

| Directory | Short name | What it does |
| --- | --- | --- |
| `vasp-relax` | `vasp.relax` | Relax the geometry of a structure with the reviewed remedy ladder (Python runner). |
| `vasp-relax-bash` | `vasp.relax-bash` | The same relaxation, authored in Bash. |
| `vasp-static` | `vasp.static` | One single-point total-energy calculation of a fixed structure. |
| `vasp-relax-static` | `vasp.relax-static` | Relax, promote the relaxed structure, then run it statically. |

These are package copies, not new implementations: each `run` is a byte-for-byte
copy of the matching built-in runner in *httk-workflow*
(`httk.workflow.vasp.runners`), and each `collect.py` calls the matching
built-in `httk.workflow.vasp.collect` function. Running a job still requires an
installed *httk-workflow*, since the runners import `httk.workflow.vasp` for
inputs, remedies, diagnostics, and collection.

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
httk plugin install 'git+https://github.com/httk/workflows-vasp.git'
```

installs all four workflow packages listed in `httk_plugin.toml` at once.
