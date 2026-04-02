#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REGRESSION_ROOT="${SCRIPT_DIR}/regressionTests"
CASE_DIR="${REGRESSION_ROOT}/main"

# ============================================================
# Plate2D regression test
# Uses a shortened run and checks the final peak temperature and
# deposited laser power to confirm the thermal/ray-tracing path.
# ============================================================

REG_END_TIME=0.001

T_MAX_MIN=4375
T_MAX_MAX=4390
Q0_MIN=318
Q0_MAX=320
Q1_MIN=426
Q1_MAX=429

ALLRUN_LOGFILE="log.Allrun"
SOLVER_LOGFILE="log.laserbeamFoam"
T_MINMAX_FILE="postProcessing/TMinMax/0/fieldMinMax.dat"

echo "============================================================"
echo "Plate2D regression test"
echo "Regression end time = ${REG_END_TIME}"
echo "Final max(T) in [${T_MAX_MIN}, ${T_MAX_MAX}] K"
echo "Final laser0 Q in [${Q0_MIN}, ${Q0_MAX}]"
echo "Final laser1 Q in [${Q1_MIN}, ${Q1_MAX}]"
echo "============================================================"
echo

prepare_case() {
    rm -rf "${CASE_DIR}"
    mkdir -p "${CASE_DIR}"

    for item in "${SCRIPT_DIR}"/*; do
        base_item=$(basename "${item}")
        if [[ "${base_item}" == "regressionTests" ]]; then
            continue
        fi
        cp -a "${item}" "${CASE_DIR}/"
    done
}

shorten_case() {
    local control_dict="${CASE_DIR}/system/controlDict"

    sed -i.bak "s/^endTime.*/endTime         ${REG_END_TIME};/" "${control_dict}"
    rm -f "${control_dict}.bak"

    cat >> "${control_dict}" <<'EOF'

functions
{
    TMinMax
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        fields          (T);
        writeControl    timeStep;
        writeInterval   1;
        write           true;
        log             false;
        mode            magnitude;
    }
}
EOF
}

extract_final_tmax() {
    awk -F '\t+' 'END {print $5}' "${CASE_DIR}/${T_MINMAX_FILE}"
}

extract_final_q_values() {
    awk '/Total Q deposited:/ {print $4}' "${CASE_DIR}/${SOLVER_LOGFILE}" | tail -n 2
}

CHECK_ONLY=false

for arg in "$@"; do
    case "$arg" in
        --check-only|--no-run)
            CHECK_ONLY=true
            ;;
        *)
            ;;
    esac
done

if [ "${CHECK_ONLY}" = false ]; then
    prepare_case
    shorten_case
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
    ( cd "${CASE_DIR}" && ./Allrun > "${ALLRUN_LOGFILE}" 2>&1 )
else
    echo "Running in check-only mode: skipping Allclean and Allrun"
fi

if [[ ! -f "${CASE_DIR}/${T_MINMAX_FILE}" ]]; then
    echo "FAIL: Could not find ${T_MINMAX_FILE}"
    exit 1
fi

if [[ ! -f "${CASE_DIR}/${SOLVER_LOGFILE}" ]]; then
    echo "FAIL: Could not find ${SOLVER_LOGFILE}"
    exit 1
fi

final_tmax=$(extract_final_tmax)

q_values=($(extract_final_q_values))
final_q0="${q_values[0]:-}"
final_q1="${q_values[1]:-}"

if [[ -z "${final_tmax}" || -z "${final_q0}" || -z "${final_q1}" ]]; then
    echo "FAIL: Could not extract one or more regression values"
    exit 1
fi

failures=0

if awk "BEGIN {exit !(${final_tmax} >= ${T_MAX_MIN} && ${final_tmax} <= ${T_MAX_MAX})}"; then
    printf "PASS: Final max(T) = %.6f K\n" "${final_tmax}"
else
    printf "FAIL: Final max(T) = %.6f K\n" "${final_tmax}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_q0} >= ${Q0_MIN} && ${final_q0} <= ${Q0_MAX})}"; then
    printf "PASS: Final laser0 Q = %.6f\n" "${final_q0}"
else
    printf "FAIL: Final laser0 Q = %.6f\n" "${final_q0}"
    failures=$((failures + 1))
fi

if awk "BEGIN {exit !(${final_q1} >= ${Q1_MIN} && ${final_q1} <= ${Q1_MAX})}"; then
    printf "PASS: Final laser1 Q = %.6f\n" "${final_q1}"
else
    printf "FAIL: Final laser1 Q = %.6f\n" "${final_q1}"
    failures=$((failures + 1))
fi

if [ "${CHECK_ONLY}" = false ]; then
    ( cd "${CASE_DIR}" && ./Allclean > /dev/null 2>&1 ) || true
fi

echo
if (( failures == 0 )); then
    echo "============================================================"
    echo "Regression test PASSED"
    echo "============================================================"
    exit 0
else
    echo "============================================================"
    echo "Regression test FAILED (${failures} checks)"
    echo "============================================================"
    exit 1
fi
