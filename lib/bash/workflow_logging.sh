# Compact phase messages for Snakemake jobs; detailed tool output stays in log/.
# Source this file, then call trm_begin RULE DESCRIPTION WORKDIR [LOG ...].

trm_emit() {
    local level="$1" message="$2" stamp line
    printf -v stamp '%(%Y-%m-%d %H:%M:%S%z)T' -1
    message="${message//$'\n'/ }"
    message="${message//$'\r'/ }"
    printf -v line '%s [TrEMOLO][%s][%s][%s] %s' "$stamp" "$trm_sample" "$trm_rule" "$level" "$message"
    printf '%s\n' "$line" >&"$trm_global_fd" || :
    # O_APPEND, one short write per message; a logging failure must not hide the
    # original command's status or invalidate otherwise valid scientific output.
    printf '%s\n' "$line" >> "$trm_progress_file" 2>/dev/null || :
}

trm_step() {
    trm_phase="$1"
    trm_emit STEP "$trm_phase"
}

trm_finish() {
    local code="$1" detail=""
    trap - EXIT
    if [[ "$code" -eq 0 ]]; then
        trm_emit DONE "$trm_description"
    else
        trm_emit ERROR "$trm_phase failed (exit $code). Detailed log: ${trm_logs[*]:-(none)}"
        detail=$(python3 "$trm_error_helper" --since "$trm_started" "${trm_logs[@]}" 2>/dev/null) || detail=""
        if [[ -n "$detail" ]]; then trm_emit ERROR "$detail"; fi
    fi
    return "$code"
}

trm_begin() {
    trm_rule="$1"
    trm_description="$2"
    local directory="$3"
    shift 3
    trm_logs=("$@")
    local index
    for index in "${!trm_logs[@]}"; do
        if [[ "${trm_logs[$index]}" != /* ]]; then trm_logs[$index]="$PWD/${trm_logs[$index]}"; fi
    done
    trm_sample="${directory%/}"
    trm_sample="${trm_sample##*/}"
    trm_phase="$trm_description"
    if [[ "$directory" != /* ]]; then directory="$PWD/$directory"; fi
    trm_progress_file="$directory/log/progress.log"
    trm_error_helper="${BASH_SOURCE[0]%/*}/../python/workflow/workflow_log_error.py"
    printf -v trm_started '%(%s)T' -1
    exec {trm_global_fd}>&2
    if [[ ! -d "$directory/log" ]]; then mkdir -p -- "$directory/log"; fi
    trap 'trm_finish "$?"' EXIT
    trm_emit START "$trm_description"
}
