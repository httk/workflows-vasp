"""The packaged VASP workflow directories, driven end to end through the real manager.

Ported from httk-workflow's ``tests/test_vasp_runners.py`` (``git show
15b0b36:tests/test_vasp_runners.py``), now that the four workflows live here as
package directories instead of built-in registrations. A fake ``vasp`` writes
plausible OSZICAR, OUTCAR, CONTCAR, and ``vasprun.xml`` files (and a CHGCAR when
the INCAR asks for one), logs every INCAR it ran with to ``fake-vasp-runs.jsonl``,
and can be scripted: a sequence of final energies, runs that use up all ``NSW``
ionic steps, and a first run that fails with a diagnosable error, which is what
exercises the remedy ladder. Nothing here simulates the workflow protocol:
every test submits a real job to a real workspace and lets a real
:class:`httk.workflow.TaskManager` run it, resolving the runner straight from
this repository's package directories with
:func:`httk.workflow.scaffold.new_job`, installing the package first — the same
call a deployment makes with ``--install --workflow-dir``.

Dropped: ``test_the_installed_package_form_resolves_the_packaged_runner``. It
exercised the reserved ``pkg:`` reference for a *pip-installed* runner package,
which does not apply to a package driven straight from its checkout directory
(a directory-sourced provider has no installed distribution to pin); the
workspace installation every other test already makes covers job resolution. ``test_every_packaged_runner_describes_itself_and_is_referenceable``
is folded into ``test_packages.py``, which already checks every package's
``--describe`` steps and workflow name against its manifest; the Bash/Python
parity this test also checked is covered more strongly below, by
``test_the_bash_runner_and_the_python_runner_publish_the_same_result``
comparing the entire job outcome rather than only the declared steps.
"""

import json
import re
import runpy
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest
from httk.workflow import TaskManager, Workspace
from httk.workflow.introspection import read_state, resolve_job
from httk.workflow.scaffold import new_job

REPO_ROOT = Path(__file__).resolve().parent.parent

_POSCAR = """silicon
1.0
2.0 0.0 0.0
0.0 2.0 0.0
0.0 0.0 2.0
Si
2
Direct
0.0000000000 0.0000000000 0.0000000000
0.5000000000 0.5000000000 0.5000000000
"""
_POTCAR = "  TITEL  = PAW_PBE Si 05Jan2001\n   ZVAL   =    4.000    mass and valenz\n"

_FAKE_VASP = '''#!/usr/bin/env python3
"""A fake VASP that writes plausible outputs, scripted by the test that installed it."""

import json
from pathlib import Path

from httk.codes.vasp import read_incar

# The final energy of each run in turn, the last one repeating.
ENERGIES = {energies!r}
# The runs that use up every NSW ionic step without converging.
IONIC_NONCONVERGED = {ionic!r}
# The runs that also report too few bands, a non-warning problem the ladder remedies.
TOO_FEW_BANDS = {too_few_bands!r}

attempts = Path("fake-vasp-attempts")
count = int(attempts.read_text()) if attempts.is_file() else 0
attempts.write_text(str(count + 1))
tags = read_incar("INCAR")
with Path("fake-vasp-runs.jsonl").open("a") as stream:
    stream.write(json.dumps(dict(INCAR=Path("INCAR").read_text(), tags=tags, CHGCAR=Path("CHGCAR").is_file())) + "\\n")
structure = Path("POSCAR").read_text().splitlines()

if {fail_when}:
    print({failure!r})
    Path("OUTCAR").write_text(" fake vasp 6.4.1\\n   NELM   =     60\\n   NSW    =     99\\n")
    Path("OSZICAR").write_text("DAV:   1    -0.100000000000E+02\\n")
    raise SystemExit(1)

if count in TOO_FEW_BANDS:
    print("TOO FEW BANDS")
energy = ENERGIES[min(count, len(ENERGIES) - 1)]
nsw = int(tags.get("NSW", "0"))
toten = {static_energy} if tags.get("NSW") == "0" else energy
Path("OUTCAR").write_text(
    " fake vasp 6.4.1\\n"
    "   NELM   =     60;   NELMIN=  2; NELMDL= -5\\n"
    "   NSW    =     %d    number of steps for IOM\\n"
    "   maximum number of plane-waves:    1234\\n"
    "   FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)\\n"
    "   free  energy   TOTEN  =       %.8f eV\\n"
    "   energy  without entropy=      %.8f  energy(sigma->0) =      %.8f\\n"
    " General timing and accounting informations for this job:\\n" % (nsw, toten, toten, toten)
)
oszicar = ["       N       E                     dE             d eps       ncg     rms"]
for step in range(1, (nsw if count in IONIC_NONCONVERGED else 1) + 1):
    oszicar.append("DAV:   1    -0.100000000000E+02   -0.10000E+02   -0.30000E+01   128   0.500E+01")
    oszicar.append("DAV:   2    %.12E   -0.50000E+00   -0.10000E+00   128   0.100E+00" % energy)
    oszicar.append("   %d F= %.8E E0= %.8E  d E =-.105000E+02" % (step, energy, energy))
Path("OSZICAR").write_text("\\n".join(oszicar) + "\\n")
relaxed = list(structure)
relaxed[0] = "relaxed by the fake vasp"
relaxed[-1] = "0.5100000000 0.5100000000 0.5100000000"
Path("CONTCAR").write_text("\\n".join(relaxed) + "\\n")
Path("vasprun.xml").write_text(
    '<modeling><structure name="finalpos"><crystal>'
    '<i name="volume">      8.00000000 </i></crystal></structure></modeling>\\n'
)
if tags.get("LCHARG", "").strip(".").upper().startswith("T"):
    Path("CHGCAR").write_text("fake charge density\\n")
'''

_COLLECTED = (
    "INCAR",
    "KPOINTS",
    "OUTCAR",
    "CONTCAR",
    "OSZICAR",
    "vasprun.xml",
    "vasp-run-report.json",
)


def _fake_vasp(
    root: Path,
    *,
    fail_once: bool = False,
    fail_when: str | None = None,
    failure: str = "LAPACK: Routine ZPOTRF failed!",
    energies: tuple[float, ...] = (-10.5,),
    ionic: tuple[int, ...] = (),
    too_few_bands: tuple[int, ...] = (),
    static_energy: float = -10.5,
) -> Path:
    """Install the fake VASP executable and return its path.

    The first run fails with *failure* when *fail_once* is set; *fail_when* is a
    Python condition over the run index ``count`` and the INCAR ``tags`` instead.
    Run ``count`` ends at ``energies[count]`` (the last one repeating), and the runs
    in *ionic* use up all their ``NSW`` ionic steps, and those in *too_few_bands*
    also print VASP's TOO FEW BANDS warning. A single point (``NSW = 0``)
    reports *static_energy* in its OUTCAR.
    """

    path = root / ("fail-once-vasp" if fail_once else "fake-vasp")
    path.write_text(
        _FAKE_VASP.format(
            fail_when=fail_when or ("count == 0" if fail_once else "False"),
            failure=failure,
            energies=list(energies),
            ionic=list(ionic),
            too_few_bands=list(too_few_bands),
            static_energy=static_energy,
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _campaign(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    directory: str,
    *,
    parameters: dict[str, object] | None = None,
    fail_once: bool = False,
    fake: dict[str, Any] | None = None,
    files: tuple[str, ...] = ("POSCAR", "POTCAR"),
    extra_files: dict[str, str] | None = None,
    publish_data: bool = True,
    initial_step: str = "prepare",
    command: str | None = None,
    workspace_settings: dict[str, object] | None = None,
    set_command_environment: bool = True,
    bare_runner: bool = False,
) -> tuple[Workspace, str]:
    """Submit and run one job of one packaged workflow directory, and return where it landed.

    ``publish_data`` sets the runners' opt-in for a curated ``data/`` copy of the results.
    ``fake`` scripts the fake VASP (see :func:`_fake_vasp`), and ``extra_files``
    adds payload files by name and content, such as an ``INCAR.<stage>`` override.
    ``bare_runner`` targets the directory's runner entry file directly instead of the
    package directory, which skips the manifest's declared-input checks (a bare
    runner file carries no input metadata) — used by the one test that means to
    reach the runner's own defensive "missing input" check instead.
    """

    root.mkdir(parents=True)
    executable = _fake_vasp(root, fail_once=fail_once, **(fake or {}))
    if set_command_environment:
        monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable) if command is None else command)
    workspace = Workspace.initialize(root / "workspace")
    for key, value in (workspace_settings or {}).items():
        workspace.set_setting(key, value)
    staged = root / "inputs"
    staged.mkdir()
    stage: dict[str, Path] = {}
    contents = {name: text for name, text in (("POSCAR", _POSCAR), ("POTCAR", _POTCAR)) if name in files}
    for name, content in {**contents, **(extra_files or {})}.items():
        path = staged / name
        path.write_text(content, encoding="utf-8")
        stage[name] = path
    target = (
        REPO_ROOT / directory / ("run.sh" if directory == "vasp-relax-bash" else "run.py")
        if bare_runner
        else REPO_ROOT / directory
    )
    job = new_job(
        workspace,
        target,
        files=stage,
        parameters={"publish_data": True, **(parameters or {})} if publish_data else parameters or {},
        step=initial_step,
        name=f"packaged {directory}",
        install=True,
    )
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=300.0)
    return workspace, job.job_id


def _runs(workdir: Path) -> list[dict[str, Any]]:
    """Every fake VASP run of one job, in order: its INCAR text and tags, and whether a CHGCAR was there."""

    return [json.loads(line) for line in (workdir / "fake-vasp-runs.jsonl").read_text(encoding="utf-8").splitlines()]


def _template(directory: str, name: str) -> str:
    return (REPO_ROOT / directory / "templates" / name).read_text(encoding="utf-8")


def _payload_of(workspace: Workspace, job_id: str) -> tuple[str, Path]:
    """Return the state and the current directory of one job."""

    ref = resolve_job(workspace, job_id)
    return ref.state, ref.path


def _failure(workspace: Workspace, job_id: str) -> Any:
    state, damaged = read_state(resolve_job(workspace, job_id))
    assert damaged is None and state is not None and state.failure is not None
    return state.failure


def _job_state(payload: Path) -> dict[str, Any]:
    path = payload / ".httk-job" / "state.json"
    return {} if not path.is_file() else json.loads(path.read_text(encoding="utf-8"))


def _files(root: Path) -> list[str]:
    """Every regular file below *root*, excluding runner-private bookkeeping."""

    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and ".httk-runner" not in path.parts
    )


def test_the_packaged_relax_runner_runs_the_v1_stages_from_their_templates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The regression: without a payload INCAR the relaxation used to start from an
    # empty INCAR, so NSW stayed at VASP's 0 and nothing relaxed.
    workspace, job_id = _campaign(
        tmp_path / "relax",
        monkeypatch,
        "vasp-relax",
        parameters={"incar_tags": {"ISYM": 0}},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"

    # preprerelax, prerelax, relax1 and relax2, each with its template's ionic loop;
    # relax2 landed within 10·EDIFFG of relax1 and stays at the top of the workdir.
    runs = _runs(workdir)
    assert [(run["tags"]["NSW"], run["tags"]["IBRION"], run["tags"]["ISIF"]) for run in runs] == [
        ("20", "2", "7"),
        ("20", "2", "3"),
        ("20", "1", "3"),
        ("20", "1", "3"),
    ]
    assert sorted(path.name for path in (workdir / "stages").iterdir()) == ["preprerelax", "prerelax", "relax1"]
    assert (workdir / "stages" / "prerelax" / "INCAR").read_text(encoding="utf-8") == runs[1]["INCAR"]
    assert (workdir / "INCAR").read_text(encoding="utf-8") == runs[3]["INCAR"]

    # The prerelax template sets its own EDIFF and EDIFFG, which are kept; the relax
    # template sets neither, so they are derived with the relaxation margin of 500.
    assert (runs[1]["tags"]["EDIFF"], runs[1]["tags"]["EDIFFG"]) == ("1E-3", "-0.3")
    for run in runs[2:]:
        assert float(run["tags"]["EDIFF"]) == pytest.approx(0.001 * 2 / 500)
        assert float(run["tags"]["EDIFFG"]) == pytest.approx(0.001 * 2)
    # relax1 does not read prerelax's CHGCAR, written on another grid (PREC = Normal,
    # default ENCUT); relax2 reads relax1's, written on its own grid.
    assert [(run["tags"]["ICHARG"], run["CHGCAR"]) for run in runs[2:]] == [("2", True), ("1", True)]

    # The job's incar_tags reach every stage, rewritten rather than duplicated.
    assert all(run["tags"]["ISYM"] == "0" for run in runs)
    assert runs[3]["INCAR"].count("ISYM") == 1
    for tag in ("MAGMOM", "NBANDS"):
        assert f"{tag} = " in runs[3]["INCAR"]
    # The pre-stages run at prerelax_kpoint_density, the relax rounds at kpoint_density.
    assert (workdir / "stages" / "prerelax" / "KPOINTS").read_text(encoding="utf-8").splitlines()[3] == "6 6 6"
    assert (workdir / "KPOINTS").read_text(encoding="utf-8").splitlines()[2:4] == ["Monkhorst-Pack", "11 11 11"]

    state = _job_state(payload)
    assert state["classification"] == "completed"
    assert float(str(state["energy"])) == pytest.approx(-10.5)
    assert (state["stage"], state["relax_round"]) == ("relax", 2)
    assert state["chgcar_grid"] == "accurate;600;900"
    assert "remedies" not in state
    assert (workdir / "fake-vasp-attempts").read_text(encoding="utf-8") == "4"

    # And the finished calculation was published into data/ (publish_data).
    published = _files(payload / "data")
    assert published == [f"vasp/{name}" for name in sorted(_COLLECTED)]
    assert (payload / "data" / "vasp" / "CONTCAR").read_text(encoding="utf-8").splitlines()[-1].startswith("0.51")


def test_relax_rounds_continue_until_the_energy_settles(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # preprerelax uses up its ionic steps (accepted), relax2 moves by 0.5 eV, relax3
    # uses up its ionic steps, and relax4 lands within 10·EDIFFG = 0.02 eV of relax3.
    workspace, job_id = _campaign(
        tmp_path / "rounds",
        monkeypatch,
        "vasp-relax",
        fake={"energies": (-9.0, -10.0, -10.0, -10.5, -10.505, -10.506), "ionic": (0, 4)},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    assert sorted(path.name for path in (workdir / "stages").iterdir()) == [
        "preprerelax",
        "prerelax",
        "relax1",
        "relax2",
        "relax3",
    ]
    report = json.loads((workdir / "stages" / "preprerelax" / "vasp-run-report.json").read_text(encoding="utf-8"))
    assert report["classification"] == "nonconverged"
    state = _job_state(payload)
    assert (state["relax_round"], float(str(state["energy"]))) == (4, pytest.approx(-10.506))
    assert "remedies" not in state


@pytest.mark.parametrize("directory", ("vasp-relax", "vasp-relax-bash"))
def test_a_relaxation_that_never_settles_fails_by_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, directory: str, test_profile
) -> None:
    if directory == "vasp-relax-bash" and not test_profile.extended:
        pytest.skip("the Bash runner's variant only runs under HTTK_TEST_PROFILE=extended")
    workspace, job_id = _campaign(
        tmp_path / "unsettled",
        monkeypatch,
        directory,
        parameters={"maximum_relax_rounds": 2},
        fake={"energies": (-9.0, -10.0, -11.0, -12.0, -13.0)},
    )

    assert _payload_of(workspace, job_id)[0] == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.relax_not_converged"
    assert failure["message"] == (
        "the relaxation did not converge in rounds relax1 to relax3: relax3 changed the energy by 1 eV (limit 0.02)"
    )
    assert failure["details"]["energy_change"] == pytest.approx(1.0)
    payload = _payload_of(workspace, job_id)[1]
    assert len(_runs(payload / "run")) == 5


@pytest.mark.parametrize("directory", ("vasp-relax", "vasp-relax-bash"))
def test_a_remedy_stays_applied_in_every_later_stage(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, directory: str, test_profile
) -> None:
    if directory == "vasp-relax-bash" and not test_profile.extended:
        pytest.skip("the Bash runner's variant only runs under HTTK_TEST_PROFILE=extended")
    workspace, job_id = _campaign(
        tmp_path / "sticky",
        monkeypatch,
        directory,
        fail_once=True,
        fake={
            "failure": "VERY BAD NEWS! internal error in subroutine IBZKPT: inverse of rotation matrix was not found"
        },
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    runs = _runs(payload / "run")
    # preprerelax failed, was remedied with ISYM = 0 in place and rerun; every stage
    # prepared afresh from its template after that re-applied the remedy.
    assert "ISYM" not in runs[0]["tags"]
    assert [run["tags"].get("ISYM") for run in runs[1:]] == ["0"] * 4
    assert runs[-1]["tags"]["IBRION"] == "1"
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"inverse_rotation": 1}
    assert _job_state(payload)["remedies"] == 1


def test_a_kpoint_remedy_is_replayed_into_every_later_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "kpoints",
        monkeypatch,
        "vasp-relax",
        fail_once=True,
        fake={"failure": "VERY BAD NEWS! Reciprocal lattice and k-lattice belong to different class of lattices."},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"kpoints_class": 1}
    # The ladder's bump_kpoints raised preprerelax's 6 6 6 grid, and every later
    # stage gets its own grid raised by one too: relax's 11 11 11 became 12 12 12.
    assert (workdir / "stages" / "preprerelax" / "KPOINTS").read_text(encoding="utf-8").splitlines()[3] == "7 7 7"
    assert (workdir / "KPOINTS").read_text(encoding="utf-8").splitlines()[3] == "12 12 12"


def test_a_runner_installed_without_its_templates_fails_by_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A bare runner file is installed alone, without the templates/ beside it.
    workspace, job_id = _campaign(tmp_path / "bare", monkeypatch, "vasp-relax", bare_runner=True)

    assert _payload_of(workspace, job_id)[0] == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.template_missing"
    assert failure["details"]["expected"].endswith("templates/INCAR.preprerelax")
    assert "package directory" in failure["message"]


def test_a_retried_stage_step_takes_the_charge_density_grid_from_job_state(tmp_path: Path) -> None:
    # Driven through the step function directly: an attempt that dies between the
    # template copy and its advance cannot be provoked cheaply through the manager.
    # Its retry finds relax's template already in the workdir next to prerelax's
    # CHGCAR, so the grid that CHGCAR was written on has to come from job state.
    from httk.codes.vasp import read_incar

    stage = runpy.run_path(str(REPO_ROOT / "vasp-relax" / "run.py"))["stage"]
    workdir = tmp_path / "run"
    workdir.mkdir()
    (workdir / "POSCAR").write_text(_POSCAR, encoding="utf-8")
    (workdir / "POTCAR").write_text(_POTCAR, encoding="utf-8")
    (workdir / "CHGCAR").write_text("prerelax charge density\n", encoding="utf-8")
    shutil.copyfile(REPO_ROOT / "vasp-relax" / "templates" / "INCAR.relax", workdir / "INCAR")
    # prerelax's grid (PREC = Normal, default ENCUT and ENAUG), then relax's own as a control.
    for grid, icharg in (("normal;;", "2"), ("accurate;600;900", "1")):
        advanced: list[str] = []
        attempt = SimpleNamespace(
            state={"stage": "relax", "relax_round": 1, "chgcar_grid": grid},
            payload=tmp_path,
            workdir=workdir,
            parameter=lambda name, default=None: default,
            setting=lambda name, default=None: default,
            log=SimpleNamespace(append=lambda *arguments: None),
            advance=advanced.append,
        )
        stage(attempt)
        assert advanced == ["run"]
        assert read_incar(workdir / "INCAR")["ICHARG"] == icharg


def test_a_payload_incar_replaces_a_stage_template(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # A prerelax that writes no CHGCAR leaves relax1 nothing to read, so its
    # ICHARG = 1 falls back to ICHARG = 2; relax2 reads the CHGCAR relax1 wrote.
    workspace, job_id = _campaign(
        tmp_path / "override",
        monkeypatch,
        "vasp-relax",
        extra_files={
            "INCAR.prerelax": _template("vasp-relax", "INCAR.prerelax").replace("LCHARG=.TRUE.", "LCHARG=.FALSE."),
            "INCAR.relax": _template("vasp-relax", "INCAR.relax").replace("ENCUT=600", "ENCUT=520"),
        },
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    runs = _runs(payload / "run")
    assert runs[1]["tags"]["LCHARG"] == ".FALSE."
    assert [run["tags"]["ENCUT"] for run in runs[2:]] == ["520", "520"]
    assert [(run["tags"]["ICHARG"], run["CHGCAR"]) for run in runs[2:]] == [("2", False), ("1", True)]


def test_a_diagnosed_failure_is_remedied_and_the_rerun_succeeds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace, job_id = _campaign(tmp_path / "remedy", monkeypatch, "vasp-relax", fail_once=True)

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    state = _job_state(payload)

    # One remedy was applied, and it was the first rung of the reviewed zpotrf
    # ladder: the lattice was scaled by five percent, and the stages after it
    # carried the scaled structure forward.
    assert state["remedies"] == 1
    assert state["classification"] == "completed"
    assert (workdir / "POSCAR").read_text(encoding="utf-8").splitlines()[1] == "1.05"
    assert (workdir / "fake-vasp-attempts").read_text(encoding="utf-8") == "5"

    # The ladder itself is job state, outside every workdir, so an isolated workdir
    # would find it too.
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text(encoding="utf-8"))
    assert history["attempts"] == {"zpotrf": 1}
    assert history["events"][0]["problem"] == "zpotrf"
    assert history["events"][0]["files"][0]["path"] == "POSCAR"
    assert not (workdir / ".httk-vasp").exists()


@pytest.mark.parametrize("recovery", ("decomposition", "bands", "never"))
def test_zhegv_manager_retries_each_rung_and_stops(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, recovery: str
) -> None:
    # Simulate acceptance of the edited inputs, not VASP's eigensolver physics.
    condition = (
        f"{recovery == 'never'!r} or tags.get('NPAR') != '1' or ({recovery == 'bands'!r} and int(tags['NBANDS']) < 8)"
    )
    workspace, job_id = _campaign(
        tmp_path / "zhegv",
        monkeypatch,
        "vasp-relax",
        fake={
            "fail_when": condition,
            "failure": "| EDDAV: Call to ZHEGV failed. Returncode = 42 2 64 |\n"
            "| I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! |",
        },
        parameters={"incar_tags": {"NPAR": 32, "NCORE": 1, "NBANDS": 6}},
    )
    kind, payload = _payload_of(workspace, job_id)
    assert kind == ("failed" if recovery == "never" else "succeeded")
    inputs = [(run["tags"]["NPAR"], run["tags"]["NBANDS"]) for run in _runs(payload / "run")]
    expected = [("32", "6"), ("1", "6")]
    if recovery != "decomposition":
        expected.append(("1", "8"))
    if recovery != "never":
        # The three stages after preprerelax start from NPAR = 32 and NBANDS = 6
        # again, and run with the remedies re-applied.
        expected += [expected[-1]] * 3
    assert inputs == expected
    history = json.loads((payload / ".httk-job" / "vasp-remedies.json").read_text())
    rungs = 1 if recovery == "decomposition" else 2
    assert history["attempts"] == {"edddav_zhegv": rungs}
    assert [event["step"] for event in history["events"]] == list(range(rungs))
    assert all(event["files"][0]["path"] == "INCAR" for event in history["events"])
    if recovery == "never":
        failure = _failure(workspace, job_id)
        assert failure["code"] == "vasp.failed"
        assert failure["details"]["problem"] == "edddav_zhegv"
        assert failure["details"]["give_up"]
        assert failure["details"]["step"] == 2


@pytest.mark.parametrize(
    ("label", "fail_once", "fake"),
    (
        ("plain", False, None),
        ("remedied", True, None),
        # Out of ionic steps in preprerelax and relax3, and a large step in relax2.
        ("rounds", False, {"energies": (-9.0, -10.0, -10.0, -10.5, -10.505, -10.506), "ionic": (0, 4)}),
        # Out of ionic steps with too few bands is not accepted: it goes to the ladder.
        ("too-few-bands", False, {"ionic": (0,), "too_few_bands": (0,)}),
    ),
)
def test_the_bash_runner_and_the_python_runner_publish_the_same_result(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    label: str,
    fail_once: bool,
    fake: dict[str, Any] | None,
    test_profile,
) -> None:
    """Parity: one workflow, two languages, the same stages, job state and files."""

    if not test_profile.extended and label != "plain":
        pytest.skip("only the plain parity case runs outside HTTK_TEST_PROFILE=extended")

    observed: dict[str, dict[str, Any]] = {}
    for directory in ("vasp-relax", "vasp-relax-bash"):
        workspace, job_id = _campaign(
            tmp_path / f"{label}-{directory}",
            monkeypatch,
            directory,
            fail_once=fail_once,
            fake=fake,
            # A null tag is VASP's default: it beats the template's EDIFF the runners merge in before it.
            parameters={"kpoint_density": 15.0, "incar_tags": {"NELM": 40, "EDIFF": None}},
        )
        kind, payload = _payload_of(workspace, job_id)
        state = _job_state(payload)
        state["energy"] = round(float(str(state["energy"])), 9)
        history = payload / ".httk-job" / "vasp-remedies.json"
        ladder = {} if not history.is_file() else json.loads(history.read_text(encoding="utf-8"))["attempts"]
        observed[directory] = {
            "kind": kind,
            "state": state,
            "ladder": ladder,
            "workdir": _files(payload / "run"),
            # Every stage's INCAR as VASP saw it, and what each archived stage kept.
            "runs": [run["INCAR"] for run in _runs(payload / "run")],
            "stages": {
                name: (payload / "run" / name).read_text(encoding="utf-8")
                for name in _files(payload / "run")
                if name.startswith("stages/") and name.rsplit("/", 1)[1] in ("INCAR", "KPOINTS", "CONTCAR")
            },
            "data": _files(payload / "data"),
            "inputs": [(payload / "run" / name).read_text(encoding="utf-8") for name in ("INCAR", "KPOINTS", "POSCAR")],
            "outputs": [
                (payload / "data" / "vasp" / name).read_text(encoding="utf-8")
                for name in ("CONTCAR", "OUTCAR", "OSZICAR")
            ],
        }

    assert observed["vasp-relax"]["kind"] == "succeeded"
    if fail_once or label == "too-few-bands":
        assert observed["vasp-relax"]["state"]["remedies"] == 1
    if label == "too-few-bands":
        assert observed["vasp-relax"]["ladder"]
    assert "NELM = 40" in observed["vasp-relax"]["inputs"][0]
    assert not any(re.search(r"(?m)^\s*EDIFF\s*=", run) for run in observed["vasp-relax"]["runs"])
    assert len(observed["vasp-relax"]["runs"]) == {"plain": 4, "remedied": 5, "rounds": 6, "too-few-bands": 5}[label]
    assert observed["vasp-relax"] == observed["vasp-relax-bash"]


def test_the_bash_vasp_runner_uses_the_workspace_command_setting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("HTTK_VASP_COMMAND", raising=False)
    workspace, job_id = _campaign(
        tmp_path / "settings-command",
        monkeypatch,
        "vasp-relax-bash",
        workspace_settings={"vasp.command": str(tmp_path / "settings-command" / "fake-vasp")},
        set_command_environment=False,
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    assert (payload / "run" / "fake-vasp-attempts").read_text(encoding="utf-8") == "4"


def test_a_job_with_no_configured_vasp_command_fails_by_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "nocommand", monkeypatch, "vasp-relax", command="")

    kind, _ = _payload_of(workspace, job_id)
    assert kind == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.command_missing"
    assert "HTTK_VASP_COMMAND" in failure["message"]


def test_a_vasp_command_that_cannot_start_fails_by_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    missing = tmp_path / "no-such-vasp"
    workspace, job_id = _campaign(tmp_path / "unstartable", monkeypatch, "vasp-relax", command=str(missing))

    assert _payload_of(workspace, job_id)[0] == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.failed"
    # run_vasp starts VASP through httk-workflow-vasp's launch wrapper, which reports the missing executable on
    # stderr and exits nonzero: a process failure without a remedy, not an OSError in the runner.
    assert failure["message"] == "VASP process_failure with no remaining remedy"
    _kind, payload = _payload_of(workspace, job_id)
    assert f"cannot execute {str(missing)!r}" in (payload / "run" / "vasp.err").read_text(encoding="utf-8")


def test_a_job_whose_structure_is_not_in_its_payload_fails_by_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Driven through the package directory, the manifest's declared required
    # input ("structure") is refused at submission time instead — a stronger,
    # earlier form of the same guarantee, covered by the ValueError case below.
    # Reaching the runner's own defensive "missing input" failure code needs the
    # bare runner file, whose resolution carries no declared-input metadata.
    workspace, job_id = _campaign(
        tmp_path / "nostructure",
        monkeypatch,
        "vasp-relax",
        files=("POTCAR",),
        bare_runner=True,
    )

    assert _payload_of(workspace, job_id)[0] == "failed"
    assert _failure(workspace, job_id)["code"] == "vasp.input_missing"


def test_a_package_directory_refuses_a_job_missing_its_required_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with pytest.raises(ValueError, match=r"workflow input 'structure' is required"):
        _campaign(tmp_path / "nostructure-dir", monkeypatch, "vasp-relax", files=("POTCAR",))


def test_the_static_runner_runs_its_template_single_point(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(
        tmp_path / "static",
        monkeypatch,
        "vasp-static",
        parameters={"poscar": "files/POSCAR", "data_prefix": "single-point"},
        extra_files={"INCAR.static": _template("vasp-static", "INCAR.static").replace("ENCUT=600", "ENCUT=520")},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    # The payload's INCAR.static replaced the template, and it is a single point.
    (run,) = _runs(payload / "run")
    assert (run["tags"]["NSW"], run["tags"]["IBRION"], run["tags"]["ENCUT"]) == ("0", "-1", "520")
    # Nothing to read a charge density from: the template's ICHARG = 1 became 2.
    assert run["tags"]["ICHARG"] == "2"
    assert float(run["tags"]["EDIFF"]) == pytest.approx(0.001 * 2 / 33)
    assert _files(payload / "data") == [f"single-point/{name}" for name in sorted(_COLLECTED)]


def test_the_chain_runner_relaxes_promotes_and_runs_statically(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace, job_id = _campaign(tmp_path / "chain", monkeypatch, "vasp-relax-static")

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"

    # The staged relaxation ran first, then the single point of its relaxed structure.
    runs = _runs(workdir)
    assert [(run["tags"]["NSW"], run["tags"]["IBRION"]) for run in runs] == [
        ("20", "2"),
        ("20", "2"),
        ("20", "1"),
        ("20", "1"),
        ("0", "-1"),
    ]
    # The single point reads the charge density of the final relax round.
    assert (runs[-1]["tags"]["ICHARG"], runs[-1]["CHGCAR"]) == ("1", True)
    assert runs[-1]["tags"]["NELM"] == "60"
    assert sorted(path.name for path in (workdir / "stages").iterdir()) == ["preprerelax", "prerelax", "relax1"]

    # The final relaxation round was archived before the single point overwrote the
    # workdir, and the relaxed structure became the structure of the single point.
    assert (workdir / "relax" / "OUTCAR").is_file()
    assert (workdir / "relax" / "INCAR").read_text(encoding="utf-8") == runs[3]["INCAR"]
    poscar = (workdir / "POSCAR").read_text(encoding="utf-8").splitlines()
    assert poscar[0] == "silicon"
    assert poscar[-1].startswith("0.51")

    state = _job_state(payload)
    assert state["relax_classification"] == "completed"
    assert float(str(state["relax_energy"])) == pytest.approx(-10.5)
    assert float(str(state["static_energy"])) == pytest.approx(-10.5)

    published = _files(payload / "data")
    assert published == sorted(
        [f"relax/{name}" for name in _COLLECTED] + [f"static/{name}" for name in _COLLECTED],
    )
    assert (workdir / "fake-vasp-attempts").read_text(encoding="utf-8") == "5"


def test_the_repeated_templates_stay_identical() -> None:
    # The runners are self-contained on purpose, so vasp-relax-static and vasp-relax-bash carry their own copies.
    for directory, name in (
        ("vasp-relax", "INCAR.preprerelax"),
        ("vasp-relax", "INCAR.prerelax"),
        ("vasp-relax", "INCAR.relax"),
        ("vasp-static", "INCAR.static"),
    ):
        assert _template(directory, name) == _template("vasp-relax-static", name), name
    for name in ("INCAR.preprerelax", "INCAR.prerelax", "INCAR.relax"):
        assert _template("vasp-relax", name) == _template("vasp-relax-bash", name), name


_GENERAL_INCAR = "ENCUT = 300\nISPIN = 2 ; ISYM = 2\nKSPACING=0.3 # kept as written\n"
_GENERAL_KPOINTS = "hand-written mesh\n0\nGamma\n3 3 3\n0 0 0\n"


def test_the_general_runner_runs_the_supplied_inputs_as_given(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("httk.atomistic")
    from httk.workflow import collect

    workspace, job_id = _campaign(
        tmp_path / "general",
        monkeypatch,
        "vasp-general",
        extra_files={
            "INCAR": _GENERAL_INCAR,
            "KPOINTS": _GENERAL_KPOINTS,
            "WAVECAR": "fake wavefunctions\n",
            "POTCAR.provenance.json": "{}\n",
        },
        publish_data=False,
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    workdir = payload / "run"
    # Nothing was derived: VASP ran on the job's INCAR and KPOINTS byte for byte,
    # next to every other file the job supplied.
    (run,) = _runs(workdir)
    assert run["INCAR"] == _GENERAL_INCAR
    assert (workdir / "INCAR").read_text(encoding="utf-8") == _GENERAL_INCAR
    assert (workdir / "KPOINTS").read_text(encoding="utf-8") == _GENERAL_KPOINTS
    assert (workdir / "WAVECAR").read_text(encoding="utf-8") == "fake wavefunctions\n"
    assert (workdir / "POTCAR").read_text(encoding="utf-8") == _POTCAR
    # The runner's own bookkeeping is not taken from the job's files.
    assert not (workdir / "POTCAR.provenance.json").exists()

    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {"final_structure", "total_energy"}
    assert item.outputs["total_energy"].value == pytest.approx(-10.5)
    # A calculation that writes no CONTCAR has no final structure, and collecting
    # reports the role as unfulfilled rather than failing.
    (workdir / "CONTCAR").unlink()
    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    assert set(item.outputs) == {"total_energy"}
    assert item.unfulfilled == ("final_structure",)


def test_the_general_runner_assembles_a_potcar_when_none_is_supplied(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    library = tmp_path / "library"
    (library / "Si").mkdir(parents=True)
    (library / "Si" / "POTCAR").write_text(_POTCAR, encoding="utf-8")
    workspace, job_id = _campaign(
        tmp_path / "general-potcar",
        monkeypatch,
        "vasp-general",
        files=("POSCAR",),
        extra_files={"INCAR": _GENERAL_INCAR},
        parameters={"pseudopotential_library": str(library)},
    )

    kind, payload = _payload_of(workspace, job_id)
    assert kind == "succeeded"
    assert (payload / "run" / "POTCAR").read_text(encoding="utf-8") == _POTCAR
    provenance = json.loads((payload / "run" / "POTCAR.provenance.json").read_text(encoding="utf-8"))
    assert provenance["potentials"][0]["variant"] == "Si"


@pytest.mark.parametrize("missing", ("INCAR", "POTCAR"))
def test_the_general_runner_names_a_missing_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, missing: str
) -> None:
    workspace, job_id = _campaign(
        tmp_path / f"general-no-{missing}",
        monkeypatch,
        "vasp-general",
        files=tuple(name for name in ("POSCAR", "POTCAR") if name != missing),
        extra_files={} if missing == "INCAR" else {"INCAR": _GENERAL_INCAR},
    )

    assert _payload_of(workspace, job_id)[0] == "failed"
    failure = _failure(workspace, job_id)
    assert failure["code"] == "vasp.input_missing"
    assert list(failure["details"]["missing"]) == [missing]


@pytest.mark.parametrize("maximum_remedies", (None, 1))
def test_the_general_runner_remedies_only_when_the_job_opts_in(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, maximum_remedies: int | None
) -> None:
    workspace, job_id = _campaign(
        tmp_path / f"general-remedies-{maximum_remedies}",
        monkeypatch,
        "vasp-general",
        fail_once=True,
        extra_files={"INCAR": _GENERAL_INCAR},
        parameters={} if maximum_remedies is None else {"maximum_remedies": maximum_remedies},
    )

    kind, payload = _payload_of(workspace, job_id)
    attempts = (payload / "run" / "fake-vasp-attempts").read_text(encoding="utf-8")
    if maximum_remedies is None:
        # Off by default: one run, and the failure says what would have helped.
        assert kind == "failed"
        assert attempts == "1"
        failure = _failure(workspace, job_id)
        assert failure["code"] == "vasp.failed"
        assert "remedies are off" in failure["message"]
        assert failure["details"]["problem"] == "zpotrf"
        assert failure["details"]["changes"]
        assert not (payload / ".httk-job" / "vasp-remedies.json").exists()
    else:
        assert kind == "succeeded"
        assert attempts == "2"
        assert _job_state(payload)["remedies"] == 1


@pytest.mark.parametrize("directory", ("vasp-relax", "vasp-relax-bash", "vasp-static", "vasp-relax-static"))
@pytest.mark.parametrize("publish_data", (False, True), ids=("default-workdir", "publish-data-opt-in"))
def test_vasp_cli_runs_and_collects_default_workdir_or_published_results(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    directory: str,
    publish_data: bool,
    test_profile,
) -> None:
    if not test_profile.extended and publish_data:
        pytest.skip("the publish-data-opt-in half of this matrix only runs under HTTK_TEST_PROFILE=extended")

    atomistic = cast(Any, pytest.importorskip("httk.atomistic"))
    atomistic_structures = cast(Any, pytest.importorskip("httk.atomistic.entries.structures"))
    store_module = cast(Any, pytest.importorskip("httk.store"))
    UnitcellStructureView = atomistic.UnitcellStructureView
    StructureEntry = atomistic_structures.StructureEntry
    Backend = store_module.Backend
    SqlStore = store_module.SqlStore
    from httk.core import DataRecord, TotalEnergyRecord
    from httk.core.cli import CLIContext
    from httk.core.storage import content_id
    from httk.workflow import collect
    from httk.workflow.packages import load_workflow_package
    from httk.workflow.postprocessing import run_postprocess_script
    from httk.workflow.workflow_cli import command

    from conftest import register_ws

    executable = _fake_vasp(tmp_path, fail_once=False, static_energy=-11.5)
    monkeypatch.setenv("HTTK_VASP_COMMAND", str(executable))
    workspace = Workspace.initialize(tmp_path / "workspace")
    context = CLIContext("httk", tmp_path)
    name = register_ws(workspace.root)
    structure = tmp_path / "POSCAR"
    structure.write_text(_POSCAR, encoding="utf-8")
    args = [
        "job",
        "new",
        "--install",
        "--workspace",
        name,
        "--workflow-dir",
        str(REPO_ROOT / directory),
        "--input",
        f"structure={structure}",
    ]
    if publish_data:
        args += ["--parameter", "publish_data=true"]
    assert command(args, context) == 0
    key, _submitted = capsys.readouterr().out.strip().split("\t")
    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=300.0)

    # The job directory moved as the job ran; it is found by its key.
    kind, payload = _payload_of(workspace, key)
    assert kind == "succeeded"
    assert (payload / "data").exists() == publish_data
    assert set(_COLLECTED) <= set(_files(payload / "run"))
    assert _job_state(payload)["classification"] == "completed"
    if not publish_data and directory in (
        "vasp-relax",
        "vasp-relax-bash",
        "vasp-static",
    ):
        # The result files occur exactly once in the whole payload.
        # (Earlier relaxation stages are kept apart under run/stages/.)
        for filename in _COLLECTED:
            found = [path for path in payload.rglob(filename) if "stages" not in path.relative_to(payload).parts]
            assert found == [payload / "run" / filename]
    if publish_data:
        prefixes = ("relax", "static") if directory == "vasp-relax-static" else ("vasp",)
        assert _files(payload / "data") == sorted(f"{prefix}/{file}" for prefix in prefixes for file in _COLLECTED)

    # No process-wide registered provider names these package-directory
    # workflows any more (the built-in ones are gone); collect them from the
    # workflow installed in the workspace instead.
    (item,) = collect(workspace, fail_fast=True, allow_job_collector=True)
    roles = {"total_energy"} if directory == "vasp-static" else {"relaxed_structure", "total_energy"}
    energy = -11.5 if directory in ("vasp-static", "vasp-relax-static") else -10.5
    assert set(item.outputs) == roles
    assert not item.unfulfilled
    collected_energy = item.outputs["total_energy"]
    assert isinstance(collected_energy, DataRecord)
    assert collected_energy.value == pytest.approx(energy)
    assert item.record.workdir == payload / "run"
    assert (item.record.data is None) == (not publish_data)
    if "relaxed_structure" in roles:
        relaxed = item.outputs["relaxed_structure"]
        assert isinstance(relaxed, UnitcellStructureView)
        assert float(cast(Any, relaxed).sites.reduced_coords[1][0]) == pytest.approx(0.51)

    for into in (False, True):
        args = ["collect", "--workspace", name, "--allow-job-collector"]
        if into:
            args += [
                "--into",
                str(tmp_path / "results.sqlite"),
                "--id-base",
                "httk.test",
                "--no-id-ledger",
            ]
        assert command(args, context) == 0
        report, summary = [json.loads(line) for line in capsys.readouterr().out.splitlines()]
        assert report["job_key"] == key
        assert set(report["outputs"]) == roles
        assert report["unfulfilled"] == []
        assert report["missing_collector"] is None
        if into:
            # The outputs and the staged structure, which collect records as the run's input: the
            # main structure plus its exact as-written ``~file_`` alternative.
            assert set(item.inputs) == {"initial_structure"}
            assert report["run"]["edges"]["inputs"] == 1
            assert len(report["stored"]["entries"]) == len(roles) + 2
            assert report["stored"]["run"]
        assert summary["format"] == "httk-workflow-collect-summary"

    with Backend.sqlite(tmp_path / "results.sqlite") as database:
        store = SqlStore(database)
        # Storing rewrites product_of links to public structure IDs, which also
        # changes the energy's content ID. Query the stored records instead: a
        # total energy is stored as the typed TotalEnergyRecord.
        searcher = store.searcher()
        variable = searcher.variable(TotalEnergyRecord)
        energies = [row.energy for row in searcher.results(energy=variable)]
        assert len(energies) == 1
        assert energies[0].total_energy == pytest.approx(energy)
        assert energies[0].id in report["stored"]["entries"]
        if "relaxed_structure" in roles:
            assert (
                store.fetch_entry(
                    StructureEntry,
                    content_id(item.outputs["relaxed_structure"]),
                    eager=True,
                )
                is not None
            )

    if "relaxed_structure" in roles:
        provider = load_workflow_package(REPO_ROOT / directory, register=False)
        for script in ("relaxation-report", "relaxation-plot"):
            result = run_postprocess_script(provider, script, item.record)
            assert result.returncode == 0, result.stderr
            if script == "relaxation-report":
                report = json.loads((result.output_dir / "relaxation_report.json").read_text())
                assert report["final_energy"] == pytest.approx(energy)
            else:
                assert f"{energy:.8f} eV" in (result.output_dir / "relaxation_energies.svg").read_text()
