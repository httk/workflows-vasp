#!/usr/bin/env bash

# The staged robust VASP relaxation, authored in Bash: the same workflow as the
# Python vasp-relax/run.py, step for step.
#
# It is the httk v1 vasp-robust-relax-formenrg template (version 1.0 by
# Christopher Tholander): preprerelax, prerelax, then relax rounds relax1, relax2,
# ... until a round lands within ten times EDIFFG of the round before it, each
# stage prepared from its INCAR template in templates/ (or a job's
# files/INCAR.<stage>). vasp-relax/run.py documents the stages and the rules in
# full; this file implements the same ones.
#
# The two runners implement one contract — the same workflow name, the same steps,
# the same job inputs and parameters, the same job state, the same failure codes, and the same
# files in the workdir and in the published data — because they call the same
# helpers through the Bash bridge. The job inputs and parameters are documented in
# this repository's README; this file and its templates/ are the entry of this
# package directory, reference it by git URI (e.g. git+https://github.com/httk/workflows-vasp#vasp-relax-bash)
# or with --workflow-dir; a job may run either runner and get the same result. To
# edit it, copy the package directory: this file alone, without its templates/,
# fails with vasp.template_missing.
set -euo pipefail
source "$HTTK_WORKFLOW_BASH_API"
source "$HTTK_WORKFLOW_VASP_BASH_API"
httk_workflow_runner vasp.relax-bash prepare stage run publish

_vasp_collect_default="INCAR KPOINTS OUTCAR CONTCAR OSZICAR vasprun.xml vasp-run-report.json POTCAR.provenance.json"
# What stages/<stage>/ keeps of every finished stage.
_vasp_stage_archive="POSCAR INCAR KPOINTS OUTCAR OSZICAR XDATCAR CONTCAR vasprun.xml vasp-run-report.json vasp.out vasp.err"
_vasp_templates=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)/templates

# Stage the structure and POTCAR of the payload, and go on to the first stage.
step_prepare() {
    local poscar potcar expected
    poscar=$HTTK_WORKFLOW_JOB_DIR/$(httk_workflow_parameter poscar files/POSCAR)
    if [ ! -f "$poscar" ]; then
        expected=$(_vasp_json_string poscar "$poscar") || return 0
        printf '{"expected": %s}\n' "$expected" >"$HTTK_WORKFLOW_CONTROL_DIR/vasp-input-missing.json"
        httk_workflow_fail vasp.input_missing \
            "the starting structure ${poscar##*/} is not in this payload" \
            --details "@$HTTK_WORKFLOW_CONTROL_DIR/vasp-input-missing.json"
        return
    fi
    cp "$poscar" POSCAR
    potcar=$HTTK_WORKFLOW_JOB_DIR/$(httk_workflow_parameter potcar files/POTCAR)
    if [ -f "$potcar" ]; then
        cp "$potcar" POTCAR
    fi
    httk_workflow_advance stage --state stage=preprerelax --state relax_round=0
}

# Prepare the current stage's inputs from its INCAR template, and go on to run
# VASP. The preparation options are written beside the attempt, not into the
# workdir, so the workdir holds exactly the calculation.
step_stage() {
    local name label template density library tags tag value margin options parallel_tag parallel_value expected centering
    name=$(httk_workflow_state_get stage)
    label=$name
    if [ "$name" = relax ]; then
        label=relax$(httk_workflow_state_get relax_round)
    fi
    template=$HTTK_WORKFLOW_JOB_DIR/files/INCAR.$name
    if [ ! -f "$template" ]; then
        template=$_vasp_templates/INCAR.$name
    fi
    if [ ! -f "$template" ]; then
        expected=$(_vasp_json_string template "$template") || return 0
        printf '{"expected": %s}\n' "$expected" >"$HTTK_WORKFLOW_CONTROL_DIR/vasp-template-missing.json"
        httk_workflow_fail vasp.template_missing \
            "the INCAR template $template is missing: reference this workflow as its package directory (a git URI ending in #vasp-relax-bash, or --workflow-dir), not as its runner file alone" \
            --details "@$HTTK_WORKFLOW_CONTROL_DIR/vasp-template-missing.json"
        return
    fi
    cp "$template" INCAR
    if [ "$name" = relax ]; then
        density=$(httk_workflow_parameter kpoint_density 20.0)
    else
        density=$(httk_workflow_parameter prerelax_kpoint_density 10.0)
    fi
    # Only the first stage assembles a POTCAR from the library; later stages find it.
    library=
    if [ ! -f POTCAR ]; then
        library=$(httk_workflow_setting vasp.pseudo_library "$(httk_workflow_parameter pseudopotential_library '')")
        if [ "$library" = null ]; then
            library=
        fi
    fi
    # v1: a relaxation converges its electrons tighter, so its ionic steps are more accurate.
    margin=33.0
    if [ "$(httk_calc "$(httk_vasp_get_tag NSW || echo 0) > 1")" = 1 ]; then
        margin=500.0
    fi
    # The template's own EDIFF/EDIFFG are kept rather than derived; the job's tags
    # win over both. The job's incar_tags is a compact JSON object, and its members
    # come after the template's in the merged object, so they win: a JSON object
    # with a repeated name takes its last value. That includes a null member, so a
    # job's {"EDIFF": null} still removes the template's EDIFF (VASP's default).
    tags=$(httk_workflow_parameter incar_tags '{}')
    for tag in EDIFFG EDIFF; do
        if value=$(httk_vasp_get_tag "$tag"); then
            value=$(_vasp_json_string "$tag" "$value") || return 0
            if [ "$tags" = '{}' ]; then
                tags="{\"$tag\": $value}"
            else
                tags="{\"$tag\": $value, ${tags#\{}"
            fi
        fi
    done
    # String values are encoded first, so a value this runner cannot encode fails the step here.
    centering=$(_vasp_json_string centering "$(httk_workflow_parameter centering Monkhorst-Pack)") || return 0
    if [ -n "$library" ]; then
        library=$(_vasp_json_string pseudopotential_library "$library") || return 0
    fi
    parallel_tag=$(httk_workflow_parameter parallel_tag '')
    if [ "$parallel_tag" = null ]; then
        parallel_tag=
    fi
    if [ -n "$parallel_tag" ]; then
        parallel_tag=$(_vasp_json_string parallel_tag "$parallel_tag") || return 0
    fi
    parallel_value=$(httk_workflow_parameter parallel_value 0)
    options=$HTTK_WORKFLOW_CONTROL_DIR/vasp-options.json
    {
        printf '{"kpoint_density": %s, "centering": %s, "accuracy_per_atom": %s, "ediff_margin": %s, "incar_tags": %s' \
            "$density" "$centering" "$(httk_workflow_parameter accuracy_per_atom 0.001)" "$margin" "$tags"
        if [ -n "$library" ]; then
            printf ', "pseudopotential_library": %s' "$library"
        fi
        if [ -n "$parallel_tag" ]; then
            printf ', "parallel_tag": %s, "parallel_value": %s' "$parallel_tag" "$parallel_value"
        fi
        printf '}\n'
    } >"$options"
    httk_vasp_prepare --directory . --options "$options" >/dev/null
    # Remedies are sticky across stages, as in v1: what an earlier stage needed, this one gets too.
    httk_vasp_remedy_reapply --directory . >/dev/null
    # The last run's charge density is only read on the grid it was written on:
    # ICHARG = 1 stays only with a CHGCAR from an accepted run with the same PREC,
    # ENCUT and ENAUG, as recorded in job state (not the workdir INCAR, which a
    # retried stage step has already replaced with its template).
    if [ "$(httk_vasp_get_tag ICHARG || true)" = 1 ]; then
        if [ ! -s CHGCAR ] || [ "$(httk_workflow_state_get chgcar_grid || true)" != "$(_vasp_grid INCAR)" ]; then
            httk_vasp_set_tag ICHARG 2
        fi
    fi
    httk_workflow_runlog_note "prepared the $HTTK_WORKFLOW_RUNNER_WORKFLOW stage $label"
    httk_workflow_advance run
}

# Run VASP once and publish what its classified result means: the end of the
# stage when it completed or only ran out of ionic steps, another attempt when a
# remedy was applied, and vasp.failed when the reviewed ladder has nothing left to try.
step_run() {
    local command timeout status classification energy applied maximum decision policy problem message amplitude
    command=$(httk_workflow_setting vasp.command "$(httk_workflow_parameter vasp_command '')")
    if [ -z "$command" ]; then
        httk_workflow_fail vasp.command_missing \
            "no VASP command is configured: set it with httk workspace settings set --key vasp.command --value '...' WORKSPACE, or set HTTK_VASP_COMMAND on the machine that runs this job, or give the job a vasp_command parameter"
        return
    fi
    timeout=$(httk_workflow_parameter timeout 86400)
    httk_vasp_preclean --directory . --keep WAVECAR --keep CHGCAR --keep CHG >/dev/null
    status=0
    # Deliberately unquoted: the resolved VASP command is one argv string, not one path.
    # shellcheck disable=SC2086
    httk_vasp_run --directory . --timeout "$timeout" --report vasp-run-report.json -- $command || status=$?
    case $status in
        0) classification=completed ;;
        20) classification=diagnosed_stop ;;
        21) classification=nonconverged ;;
        22) classification=process_failure ;;
        124) classification=timeout ;;
        *)
            httk_workflow_fail vasp.failed "httk_vasp_run could not run VASP at all (status $status)"
            return
            ;;
    esac
    energy=$(httk_vasp_energy OSZICAR || true)
    httk_workflow_runlog_note "VASP $classification"
    if [ "$classification" = completed ]; then
        _vasp_finish_stage "$classification" "$energy" 0
        return
    fi
    # Running out of NSW ionic steps is progress, not a failure, of a relaxation
    # stage: accepted when it is the only non-warning problem the run report names.
    if [ "$classification" = nonconverged ] &&
        [ "$(httk_vasp_report_problems vasp-run-report.json)" = ionic_nonconvergence ]; then
        _vasp_finish_stage "$classification" "$energy" 1
        return
    fi
    applied=$(httk_workflow_state_get remedies || echo 0)
    maximum=$(httk_workflow_parameter maximum_remedies 8)
    decision=$HTTK_WORKFLOW_CONTROL_DIR/vasp-remedy-decision.json
    policy=$(httk_workflow_parameter remedy_policy reviewed-v1)
    status=0
    problem=$(
        httk_vasp_remedy_plan vasp-run-report.json --directory . --policy "$policy" --output "$decision"
    ) || status=$?
    if [ "$status" -ne 0 ] && [ "$status" -ne 3 ]; then
        httk_workflow_fail vasp.failed "planning a VASP remedy failed with status $status"
        return
    fi
    if [ "$status" -eq 3 ] || [ "$applied" -ge "$maximum" ]; then
        if [ "$applied" -ge "$maximum" ]; then
            message="VASP $classification after $applied remedies"
        else
            message="VASP $classification with no remaining remedy"
        fi
        _vasp_merge_run_state "$classification" "$energy" ''
        httk_workflow_fail vasp.failed "$message" --details "@$decision"
        return
    fi
    httk_vasp_remedy_apply "$decision" --directory .
    amplitude=$(httk_workflow_parameter rattle_amplitude 0)
    if [ "$(httk_calc "$amplitude > 0")" = 1 ]; then
        # The entropy is the attempt itself, so the perturbation is reproducible
        # and no two attempts of this job rattle the same way.
        httk_vasp_rattle_poscar POSCAR --amplitude "$amplitude" \
            --entropy "$(httk_workflow_context job_key):$(httk_workflow_context attempt_ordinal)"
    fi
    _vasp_merge_run_state "$classification" "$energy" "$((applied + 1))"
    httk_workflow_runlog_note "applied a remedy for $problem"
    # The stage keeps its remedied inputs: the retry runs VASP again without re-preparing.
    httk_workflow_retry "applied the $policy remedy for $problem"
}

# End the relaxation when its last round converged, or archive the stage and go
# on to the next. Arguments: the classification, the energy (or empty), and 1
# when the run only ran out of ionic steps.
_vasp_finish_stage() {
    local classification=$1 energy=$2 ionic=$3 name rounds label state ediffg limit previous change reason next file lcharg
    local details=$HTTK_WORKFLOW_CONTROL_DIR/vasp-relax-not-converged.json
    name=$(httk_workflow_state_get stage)
    rounds=$(httk_workflow_state_get relax_round || echo 0)
    label=$name
    if [ "$name" = relax ]; then
        label=relax$rounds
    fi
    state=(classification="$classification")
    if [ -n "$energy" ]; then
        state+=(energy="$energy" "${label}_energy=$energy")
    fi
    # The grid of the CHGCAR this accepted run wrote, committed with its outcome so a
    # retried stage step still knows it; a run with LCHARG = .FALSE. wrote none.
    lcharg=$(httk_vasp_get_tag LCHARG || echo .TRUE.)
    case ${lcharg#.} in
        [Ff]*) state+=(chgcar_grid=null) ;;
        *) state+=("chgcar_grid=$(_vasp_grid INCAR)") ;;
    esac
    if [ "$name" = relax ] && [ "$rounds" -ge 2 ]; then
        # VASP's own default EDIFFG is ten times EDIFF.
        ediffg=$(httk_vasp_get_tag EDIFFG || echo "10 * $(httk_vasp_get_tag EDIFF || echo 1e-4)")
        previous=$(httk_workflow_state_get "relax$((rounds - 1))_energy" || true)
        if [ "$ionic" = 0 ] && [ -n "$energy" ] && [ -n "$previous" ] &&
            [ "$(httk_calc "abs(($energy) - ($previous)) <= 10 * abs($ediffg)")" = 1 ]; then
            _vasp_advance publish "${state[@]}"
            return
        fi
        if [ "$rounds" -gt "$(httk_workflow_parameter maximum_relax_rounds 5)" ]; then
            httk_workflow_state_merge "${state[@]}"
            limit=$(httk_calc "10 * abs($ediffg)")
            change=inf
            if [ -n "$energy" ] && [ -n "$previous" ]; then
                change=$(httk_calc "abs(($energy) - ($previous))")
            fi
            if [ "$ionic" = 1 ]; then
                reason="ran out of NSW ionic steps"
                printf '{"rounds": %s, "energy_change": null, "limit": %s}\n' "$rounds" "$limit" >"$details"
            else
                reason="changed the energy by $(printf %g "$change") eV (limit $(printf %g "$limit"))"
                printf '{"rounds": %s, "energy_change": %s, "limit": %s}\n' \
                    "$rounds" "${change/inf/Infinity}" "$limit" >"$details"
            fi
            httk_workflow_fail vasp.relax_not_converged \
                "the relaxation did not converge in rounds relax1 to relax$rounds: relax$rounds $reason" \
                --details "@$details"
            return
        fi
    fi
    if ! grep -qs '[^[:space:]]' CONTCAR; then
        httk_workflow_state_merge "${state[@]}"
        httk_workflow_fail vasp.no_relaxed_structure "the stage $label left no CONTCAR to continue from"
        return
    fi
    mkdir -p "stages/$label"
    for file in $_vasp_stage_archive; do
        if [ -f "$file" ]; then
            cp "$file" "stages/$label/$file"
        fi
    done
    # The reference POSCAR only lends its comment line, which carries the MAGMOM override.
    httk_vasp_promote_contcar
    case $name in
        preprerelax) next=prerelax ;;
        *) next=relax ;;
    esac
    if [ "$next" = relax ]; then
        state+=(stage=relax relax_round=$((rounds + 1)))
    else
        state+=(stage="$next" relax_round=0)
    fi
    if [ "$ionic" = 1 ]; then
        httk_workflow_runlog_note "finished the stage $label out of ionic steps"
    else
        httk_workflow_runlog_note "finished the stage $label"
    fi
    _vasp_advance stage "${state[@]}"
}

# The grid a CHGCAR of the INCAR $1 is on: its PREC, ENCUT and ENAUG,
# lower-cased and joined by ";", an absent tag empty.
_vasp_grid() {
    local tag value grid=
    for tag in PREC ENCUT ENAUG; do
        value=$(httk_vasp_get_tag "$tag" "$1" || true)
        grid+="${value,,};"
    done
    printf '%s' "${grid%;}"
}

# The value $2 as one JSON string, with backslashes, double quotes and tabs
# escaped. A newline is refused: the step fails with vasp.input_invalid naming
# $1, and this returns 1 so the caller stops.
_vasp_json_string() {
    local value=$2
    if [[ $value == *$'\n'* ]]; then
        httk_workflow_fail vasp.input_invalid "the $1 parameter contains a newline, which this runner cannot pass on"
        return 1
    fi
    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    value=${value//$'\t'/\\t}
    printf '"%s"' "$value"
}

# Advance to step $1 with the NAME=VALUE state assignments that follow.
_vasp_advance() {
    local step=$1 item arguments=()
    shift
    for item in "$@"; do
        arguments+=(--state "$item")
    done
    httk_workflow_advance "$step" "${arguments[@]}"
}

# Publish the collected files into data/ when the job opts in with publish_data,
# or leave them in the persistent workdir.
step_publish() {
    local prefix name publish published= stages=
    prefix=$(httk_workflow_parameter data_prefix vasp)
    publish=$(httk_workflow_parameter publish_data false)
    for name in $(httk_workflow_parameter collect "$_vasp_collect_default"); do
        if [ ! -f "$name" ]; then
            continue
        fi
        if [ "$publish" = true ]; then
            httk_workflow_put "$name" "$prefix/$name" >/dev/null
        fi
        if [ -n "$published" ]; then
            published="$published, $name"
        else
            published=$name
        fi
    done
    if [ "$publish" = true ]; then
        httk_workflow_runlog_note "published to data/$prefix: ${published:-nothing}"
    else
        httk_workflow_runlog_note "kept in the workdir: ${published:-nothing}"
    fi
    if [ -d stages ]; then
        stages=$(cd stages && LC_ALL=C ls | paste -sd, - | sed 's/,/, /g')
    fi
    httk_workflow_runlog_note "earlier stages kept in the workdir under stages/: ${stages:-none}"
    httk_workflow_succeed
}

# Record exactly the run state the Python runner records: the classification, the
# energy when the run produced one, and the remedy count when one was applied.
_vasp_merge_run_state() {
    local assignments=(classification="$1")
    if [ -n "$2" ]; then
        assignments+=(energy="$2")
    fi
    if [ -n "$3" ]; then
        assignments+=(remedies="$3")
    fi
    httk_workflow_state_merge "${assignments[@]}"
}

httk_workflow_main
