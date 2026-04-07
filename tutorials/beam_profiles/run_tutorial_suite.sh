#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP_ROOT="${TMPDIR:-/tmp}/beam_profiles_runs"
MODE="smoke"
FILTER=""

CASES=(
    "laserbeamFoam/bessel_beam_example"
    "laserbeamFoam/core_ring_example"
    "laserbeamFoam/csv_example"
    "laserbeamFoam/opa_cbc_example"
    "laserbeamFoam/supergaussian_example"
    "laserbeamFoam/cailabs_ringCore/core_ref"
    "laserbeamFoam/cailabs_ringCore/ring_core"
    "compressiblelaserbeamFoam/LPBF_small_vapour"
    "compressiblelaserbeamFoam/multiComponentLaserIrradiation"
)

is_smoke_supported() {
    local case="$1"

    case "$case" in
        compressiblelaserbeamFoam/LPBF_small_vapour)
            return 1
            ;;
        compressiblelaserbeamFoam/multiComponentLaserIrradiation)
            return 1
            ;;
        *)
            return 0
            ;;
    esac
}

is_enabled() {
    local case="$1"

    case "$case" in
        laserbeamFoam/cailabs_ringCore/core_ref)
            return 1
            ;;
        laserbeamFoam/cailabs_ringCore/ring_core)
            return 1
            ;;
        *)
            return 0
            ;;
    esac
}

usage() {
    cat <<'EOF'
Usage: run_tutorial_suite.sh [--mode smoke|full] [--filter substring] [--list]

Runs the beam-profile tutorial suite from copied cases under /tmp so the
tracked tutorial directories stay clean.

Options:
  --mode smoke   Shorten runtime-related settings in the copied case (default)
  --mode full    Run the copied case exactly as configured
  --filter STR   Run only cases whose relative path contains STR
  --list         Print the case list and exit
EOF
}

list_cases() {
    for case in "${CASES[@]}"; do
        echo "$case"
    done
}

rewrite_for_smoke() {
    local case_dir="$1"

    if command -v foamDictionary >/dev/null 2>&1; then
        local control_dict="$case_dir/system/controlDict"
        local decompose_dict="$case_dir/system/decomposeParDict"

        if [[ -f "$control_dict" ]]; then
            foamDictionary -entry startFrom -set startTime "$control_dict" >/dev/null 2>&1 || true
            foamDictionary -entry startTime -set 0 "$control_dict" >/dev/null 2>&1 || true
            foamDictionary -entry endTime -set 5e-7 "$control_dict" >/dev/null 2>&1 || true
            foamDictionary -entry deltaT -set 1e-7 "$control_dict" >/dev/null 2>&1 || true
            foamDictionary -entry maxDeltaT -set 5e-7 "$control_dict" >/dev/null 2>&1 || true
            foamDictionary -entry writeInterval -set 5e-7 "$control_dict" >/dev/null 2>&1 || true
        fi

        if [[ -f "$decompose_dict" ]]; then
            foamDictionary -entry numberOfSubdomains -set "${FOAM_RUN_NP:-4}" "$decompose_dict" >/dev/null 2>&1 || true
        fi
    fi
}

copy_case() {
    local rel="$1"
    local src="$SCRIPT_DIR/$rel"
    local dst="$TMP_ROOT/$rel"

    rm -rf "$dst"
    mkdir -p "$(dirname "$dst")"
    cp -a "$src" "$dst"

    if [[ "$MODE" == "smoke" ]]; then
        rewrite_for_smoke "$dst"
    fi

    echo "$dst"
}

run_case() {
    local rel="$1"
    local run_dir
    run_dir="$(copy_case "$rel")"

    echo "==> Running $rel ($MODE) in $run_dir"
    (
        cd "$run_dir"
        bash ./Allrun
    )
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)
            MODE="$2"
            shift 2
            ;;
        --filter)
            FILTER="$2"
            shift 2
            ;;
        --list)
            list_cases
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "$MODE" != "smoke" && "$MODE" != "full" ]]; then
    echo "Unsupported mode: $MODE" >&2
    exit 1
fi

mkdir -p "$TMP_ROOT"

failures=()
skipped=()
selected_count=0

for case in "${CASES[@]}"; do
    if [[ -n "$FILTER" && "$case" != *"$FILTER"* ]]; then
        continue
    fi

    selected_count=$((selected_count + 1))

    if ! is_enabled "$case"; then
        echo "==> Skipping $case: temporarily excluded while the cailabs_ringCore setup is still being refined"
        skipped+=("$case")
        continue
    fi

    if [[ "$MODE" == "smoke" ]] && ! is_smoke_supported "$case"; then
        echo "==> Skipping $case (smoke): known unstable as a shortened smoke run"
        skipped+=("$case")
        continue
    fi

    if ! run_case "$case"; then
        failures+=("$case")
        echo "==> FAILED $case" >&2
    fi
done

echo
echo "Suite summary"
echo "  passed: $(( selected_count - ${#failures[@]} - ${#skipped[@]} ))"
echo "  skipped: ${#skipped[@]}"
echo "  failed: ${#failures[@]}"

if [[ ${#skipped[@]} -gt 0 ]]; then
    echo "  skipped cases:"
    for case in "${skipped[@]}"; do
        echo "    $case"
    done
fi

if [[ ${#failures[@]} -gt 0 ]]; then
    echo "  failed cases:" >&2
    for case in "${failures[@]}"; do
        echo "    $case" >&2
    done
    exit 1
fi
