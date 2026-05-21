#!/usr/bin/env bash
# =============================================================================
# LocusLift.sh
# Production-grade genome build liftover for PLINK and VCF datasets.
#
# Supports:
#   - Input/Output : PLINK 1 (.bed/.bim/.fam), PLINK 2 (.pgen/.pvar/.psam), VCF (.vcf/.vcf.gz)
#   - Builds       : GRCh37/hg19 <-> GRCh38/hg38 (any supported chain)
#
# Core VCF flow:
#   1) Detect input assembly and validate chain direction
#   2) Extract VCF variants and build UCSC BED
#   3) Run UCSC liftOver
#   4) Build coordinate mapping for VCF rewrite
#   5) Rewrite VCF body with new coordinates
#   6) Sort and write lifted VCF (with index)
#
# Core PLINK flow:
#   1) Detect input assembly and validate chain direction
#   2) Extract coordinates in UCSC BED format (with automatic deduplication)
#   3) Run UCSC liftOver
#   4) Build chromosome/position update files
#   5) Apply updates to dataset and write final sorted dataset
# =============================================================================
set -Eeuo pipefail
IFS=$'\n\t'

SCRIPT_NAME="$(basename "$0")"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
VERSION="1.0"
PLINK2_DOC_URL="https://www.cog-genomics.org/plink/2.0/index"
BCFTOOLS_DOC_URL="https://samtools.github.io/bcftools/bcftools.html"

# -----------------------------------------------------------------------------
# CLI defaults
# -----------------------------------------------------------------------------
INPUT_PREFIX=""
INPUT_VCF=""
CHAIN_FILE=""
LIFTOVER_BIN=""
OUTPUT_PREFIX=""
OUTPUT_VCF=""
INPUT_FORMAT=""       # auto-detect when empty
OUTPUT_FORMAT=""      # defaults to input format
VARID_TEMPLATE=""     # auto-selected from chain direction when empty
KEEP_INTERMEDIATE=false
KEEP_ALL_CHR=false
FORCE_LIFTOVER=false
FORCE_ALL_LIFTOVER=false
MODE="plink"

# -----------------------------------------------------------------------------
# Visuals (minimal/premium; degrades gracefully for non-TTY)
# -----------------------------------------------------------------------------
RED=""
GREEN=""
YELLOW=""
CYAN=""
BOLD=""
DIM=""
RESET=""

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    RED=$'\033[0;31m'
    GREEN=$'\033[0;32m'
    YELLOW=$'\033[0;33m'
    CYAN=$'\033[0;36m'
    BOLD=$'\033[1m'
    DIM=$'\033[2m'
    RESET=$'\033[0m'
fi

GLYPH_INFO="i"
GLYPH_OK="OK"
GLYPH_WARN="!"
GLYPH_ERR="X"
GLYPH_STEP=">"
RULE_LINE="------------------------------------------------------------------------------------------"

if [[ "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" =~ UTF-8|utf8 ]]; then
    GLYPH_INFO="•"
    GLYPH_OK="✓"
    GLYPH_WARN="⚠"
    GLYPH_ERR="✖"
    GLYPH_STEP="→"
    RULE_LINE="──────────────────────────────────────────────────────────────────────────────────────────"
fi

# -----------------------------------------------------------------------------
# Logging helpers
# These functions format output for the user, providing visual cues like
# timestamps, colored glyphs (success, error, warning), and text wrapping
# to keep the terminal output clean and readable.
# -----------------------------------------------------------------------------
ts() {
    date '+%H:%M:%S'
}

line() {
    printf "%s\n" "${DIM}${RULE_LINE}${RESET}"
}

_wrap() {
    echo "$*" | fold -s -w 75
}

info() {
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %s %s\n" "$DIM" "$(ts)" "$RESET" "$GLYPH_INFO" "$l"
            first=false
        else
            printf "             %s\n" "$l"
        fi
    done < <(_wrap "$*")
}

step() {
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %s %b%s%b\n" "$DIM" "$(ts)" "$RESET" "$GLYPH_STEP" "$BOLD" "$l" "$RESET"
            first=false
        else
            printf "             %b%s%b\n" "$BOLD" "$l" "$RESET"
        fi
    done < <(_wrap "$*")
}

ok() {
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %b%s%b %s\n" "$DIM" "$(ts)" "$RESET" "$GREEN" "$GLYPH_OK" "$RESET" "$l"
            first=false
        else
            printf "             %s\n" "$l"
        fi
    done < <(_wrap "$*")
}

warn() {
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %b%s%b %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$l"
            first=false
        else
            printf "             %s\n" "$l"
        fi
    done < <(_wrap "$*")
}

die() {
    # Unset error trap to prevent recursive error handling
    trap - ERR
    set +e
    RUN_FAILED=true
    write_error_log "FATAL" "$*"
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %b%s%b %s\n" "$DIM" "$(ts)" "$RESET" "$RED" "$GLYPH_ERR" "$RESET" "$l" >&2
            first=false
        else
            printf "             %s\n" "$l" >&2
        fi
    done < <(_wrap "$*")
    
    local log_msg="See error log: ${ERROR_LOG_FILE:-./${SCRIPT_NAME%.sh}.error.log}"
    first=true
    while IFS= read -r l; do
        if $first; then
            printf "%b[%s]%b %b%s%b %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$l" >&2
            first=false
        else
            printf "             %s\n" "$l" >&2
        fi
    done < <(_wrap "$log_msg")
    
    if [[ -n "${FAILED_STEP_LOG_FILE:-}" && -f "$FAILED_STEP_LOG_FILE" ]]; then
        local fail_msg="Retained failure log: $FAILED_STEP_LOG_FILE"
        first=true
        while IFS= read -r l; do
            if $first; then
                printf "%b[%s]%b %b%s%b %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$l" >&2
                first=false
            else
                printf "             %s\n" "$l" >&2
            fi
        done < <(_wrap "$fail_msg")
    fi
    exit 1
}

print_kv() {
    printf "  %-16s %s\n" "$1" "$2"
}

log_line() {
    printf "%s\n" "${RULE_LINE}" >> "$RUN_SUMMARY_FILE"
}

log_kv() {
    local key="$1"
    local val="$2"
    local first=true
    local l
    while IFS= read -r l; do
        if $first; then
            printf "  %-16s %s\n" "$key" "$l" >> "$RUN_SUMMARY_FILE"
            first=false
        else
            printf "                   %s\n" "$l" >> "$RUN_SUMMARY_FILE"
        fi
    done < <(echo "$val" | fold -s -w 71)
}

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------
usage() {
    cat <<USAGE

${BOLD}LocusLift${RESET} ${DIM}v${VERSION}${RESET} ${DIM}- Genome LiftOver Pipeline (PLINK + VCF)${RESET}
${DIM}GitHub: https://github.com/BioInUmer/LocusLift${RESET}
${DIM}${RULE_LINE}${RESET}
${BOLD}Usage${RESET}
  PLINK  ./${SCRIPT_NAME} -i <prefix> -c <chain> -o <out_prefix> [opts]
  VCF    ./${SCRIPT_NAME} -I <in.vcf[.gz]> -c <chain> -O <out.vcf[.gz]> [opts]

${BOLD}Required${RESET}
  -c       UCSC chain file  (e.g. hg38ToHg19.over.chain.gz)
  -i / -o  PLINK input / output prefix
  -I / -O  VCF  input / output file

${BOLD}Options${RESET}
  -f FORMAT        Input  format: plink1|plink2  ${DIM}[auto]${RESET}
  -t FORMAT        Output format: plink1|plink2  ${DIM}[same as -f]${RESET}
  -l PATH          liftOver binary               ${DIM}[auto-search]${RESET}
  -v TMPL          Variant-ID template           ${DIM}[default: chr:pos:ref:alt]${RESET}
  -k               Keep intermediate files
  --keep-all-chr   Retain non-primary contigs    ${DIM}[default: 1-22/X/Y]${RESET}
  --force          Skip build-compatibility check for inconclusive inferences
  --force-all      Bypass all safety checks and force liftover despite explicit mismatches
  -h               Show this help

${BOLD}Examples${RESET}
  # PLINK hg19 → hg38
  ./${SCRIPT_NAME} -i data_hg19 -c hg19ToHg38.over.chain.gz -o data_hg38

  # VCF hg38→hg19
  ./${SCRIPT_NAME} -I cohort.hg38.vcf.gz -c hg38ToHg19.over.chain.gz -O cohort.hg19.vcf.gz
  
  # Format conversion (plink1→plink2), custom binary, keep all contigs
  ./${SCRIPT_NAME} -i data -c hg19ToHg38.over.chain.gz -o lifted \\
      -f plink1 -t plink2 -l ./liftOver --keep-all-chr

${DIM}Note: do not mix -i/-o (PLINK) with -I/-O (VCF). In VCF mode -f/-t/-v are unused.
For extra information, please refer to the README.
PLINK2 docs: ${PLINK2_DOC_URL}
BCFTOOLS docs: ${BCFTOOLS_DOC_URL}${RESET}
USAGE
    exit "${1:-0}"
}

# -----------------------------------------------------------------------------
# Runtime state used by traps
# These global variables track the progress and state of the pipeline,
# ensuring that error handlers have enough context to generate useful logs.
# -----------------------------------------------------------------------------
WORK_DIR=""
CURRENT_STEP=""
LAST_LOG_FILE=""
ERROR_LOG_FILE=""
ERROR_LOG_STARTED=false
FAILED_STEP_LOG_FILE=""
RUN_SUMMARY_FILE=""
RUN_STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S %Z')"
CHAIN_BUILD_GUARD_ENABLED=false
CHAIN_COMPATIBILITY_STATUS="N/A"
DETECTED_INPUT_BUILD="N/A"
DETECTED_INPUT_BUILD_REASON="N/A"
INPUT_PRIMARY_CHR_COUNT="N/A"
INPUT_PRIMARY_CHR_LIST="N/A"
ID_REWRITTEN_COUNT=0
ID_PRESERVED_COUNT=0
RUN_FAILED=false
PLINK2_ACTIVE_OUT_PREFIX=""
PLINK2_ACTIVE_STEP_LOG=""

init_error_log_path() {
    if [[ -n "${LOG_PREFIX:-}" ]]; then
        ERROR_LOG_FILE="${LOG_PREFIX}.error.log"
    elif [[ -n "$OUTPUT_PREFIX" ]]; then
        ERROR_LOG_FILE="${OUTPUT_PREFIX}.error.log"
    else
        ERROR_LOG_FILE="./${SCRIPT_NAME%.sh}.error.log"
    fi
}

init_failure_log_path() {
    if [[ -n "${LOG_PREFIX:-}" ]]; then
        FAILED_STEP_LOG_FILE="${LOG_PREFIX}.failed-step.log"
    elif [[ -n "$OUTPUT_PREFIX" ]]; then
        FAILED_STEP_LOG_FILE="${OUTPUT_PREFIX}.failed-step.log"
    else
        FAILED_STEP_LOG_FILE="./${SCRIPT_NAME%.sh}.failed-step.log"
    fi
}

init_run_summary_path() {
    if [[ -n "${LOG_PREFIX:-}" ]]; then
        RUN_SUMMARY_FILE="${LOG_PREFIX}.log"
    elif [[ -n "$OUTPUT_PREFIX" ]]; then
        RUN_SUMMARY_FILE="${OUTPUT_PREFIX}.log"
    else
        RUN_SUMMARY_FILE="./${SCRIPT_NAME%.sh}.log"
    fi
}

write_success_log_header() {
    : > "$RUN_SUMMARY_FILE"
    log_line
    printf "LocusLift (v%s) - Genome LiftOver Pipeline (PLINK/VCF)\n" "$VERSION" >> "$RUN_SUMMARY_FILE"
    log_line
    log_kv "Timestamp start" "$RUN_STARTED_AT"
    log_kv "Timestamp end" "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    log_kv "Status" "SUCCESS"
}

persist_failure_log() {
    local source_log="${1:-}"
    if [[ -z "$source_log" || ! -f "$source_log" ]]; then
        return 0
    fi
    if [[ -z "$FAILED_STEP_LOG_FILE" ]]; then
        init_failure_log_path
    fi
    cp "$source_log" "$FAILED_STEP_LOG_FILE" 2>/dev/null || return 0
    return 0
}

purge_work_logs() {
    [[ -n "$WORK_DIR" && -d "$WORK_DIR" ]] || return 0
    local f
    for f in "$WORK_DIR"/*.log; do
        [[ -e "$f" ]] || continue
        rm -f "$f"
    done
    return 0
}

write_error_log() {
    local error_kind="$1"
    local error_message="$2"

    # Initialize error log path if not already set
    if [[ -z "$ERROR_LOG_FILE" ]]; then
        init_error_log_path
    fi

    local error_dir
    error_dir="$(dirname "$ERROR_LOG_FILE")"
    if ! mkdir -p "$error_dir" 2>/dev/null; then
        ERROR_LOG_FILE="./${SCRIPT_NAME%.sh}.error.log"
    fi

    if [[ "$ERROR_LOG_STARTED" == false ]]; then
        if {
            echo "========================================"
            echo "LocusLift Error Log - Genome LiftOver Pipeline"
            echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S %Z')"
            echo "Script: ${SCRIPT_NAME} (v${VERSION})"
            echo "Input prefix: ${INPUT_PREFIX:-N/A}"
            echo "Output prefix: ${OUTPUT_PREFIX:-N/A}"
            echo "Mode: ${MODE:-N/A}"
            echo "Current step: ${CURRENT_STEP:-N/A}"
            if [[ "${DETECTED_INPUT_BUILD:-N/A}" != "N/A" ]]; then
                echo "Detected build: ${DETECTED_INPUT_BUILD:-N/A}"
                echo "Build evidence: ${DETECTED_INPUT_BUILD_REASON:-N/A}"
                echo "Chain source build: ${SOURCE_BUILD:-N/A}"
                echo "Chain target build: ${TARGET_BUILD:-N/A}"
                echo "Assembly check: ${CHAIN_COMPATIBILITY_STATUS:-N/A}"
            fi
            echo "========================================"
        } > "$ERROR_LOG_FILE" 2>/dev/null; then
            ERROR_LOG_STARTED=true
        else
            ERROR_LOG_FILE="./${SCRIPT_NAME%.sh}.error.log"
            if {
                echo "========================================"
                echo "LocusLift Error Log - Genome LiftOver Pipeline"
                echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S %Z')"
                echo "Script: ${SCRIPT_NAME} (v${VERSION})"
                echo "Input prefix: ${INPUT_PREFIX:-N/A}"
                echo "Output prefix: ${OUTPUT_PREFIX:-N/A}"
                echo "Mode: ${MODE:-N/A}"
                echo "Current step: ${CURRENT_STEP:-N/A}"
                if [[ "${DETECTED_INPUT_BUILD:-N/A}" != "N/A" ]]; then
                    echo "Detected build: ${DETECTED_INPUT_BUILD:-N/A}"
                    echo "Build evidence: ${DETECTED_INPUT_BUILD_REASON:-N/A}"
                    echo "Chain source build: ${SOURCE_BUILD:-N/A}"
                    echo "Chain target build: ${TARGET_BUILD:-N/A}"
                    echo "Assembly check: ${CHAIN_COMPATIBILITY_STATUS:-N/A}"
                fi
                echo "========================================"
            } > "$ERROR_LOG_FILE" 2>/dev/null; then
                ERROR_LOG_STARTED=true
            fi
        fi
    fi

    if [[ "$ERROR_LOG_STARTED" == true ]]; then
        {
            echo
            echo "[$(date '+%Y-%m-%d %H:%M:%S %Z')] ${error_kind}"
            echo "Message: ${error_message}"
            [[ -n "${CURRENT_STEP:-}" ]] && echo "Step: ${CURRENT_STEP}"
            if [[ "${DETECTED_INPUT_BUILD:-N/A}" != "N/A" ]]; then
                echo "Detected build: ${DETECTED_INPUT_BUILD:-N/A}"
                echo "Build evidence: ${DETECTED_INPUT_BUILD_REASON:-N/A}"
                echo "Chain source build: ${SOURCE_BUILD:-N/A}"
                echo "Chain target build: ${TARGET_BUILD:-N/A}"
                echo "Assembly check: ${CHAIN_COMPATIBILITY_STATUS:-N/A}"
            fi
            if [[ "${MODE:-}" == "vcf" ]]; then
                echo "ID rewrite counts: rewritten=${ID_REWRITTEN_COUNT:-0}, preserved=${ID_PRESERVED_COUNT:-0}"
            fi
            if [[ -n "${LAST_LOG_FILE:-}" && -f "$LAST_LOG_FILE" ]]; then
                echo "PLINK log: ${LAST_LOG_FILE}"
                echo "--- PLINK log tail (last 40 lines) ---"
                tail -n 40 "$LAST_LOG_FILE" 2>/dev/null || true
                echo "--- End PLINK log tail ---"
            fi
        } >> "$ERROR_LOG_FILE" 2>/dev/null
    fi

    return 0
}

preflight_liftover_binary() {
    # Checks if the provided UCSC liftOver binary is executable and compatible
    # with the current OS and CPU architecture.
    local os arch bin_desc probe_output probe_status
    os="$(uname -s 2>/dev/null || echo unknown)"
    arch="$(uname -m 2>/dev/null || echo unknown)"

    if command -v file >/dev/null 2>&1; then
        bin_desc="$(file -b "$LIFTOVER_BIN" 2>/dev/null || true)"
    else
        bin_desc=""
    fi

    # Temporarily disable error trap for probing the binary
    trap - ERR
    set +e
    probe_output="$("$LIFTOVER_BIN" 2>&1)"
    probe_status=$?
    set -e
    trap 'on_error "$?" "$LINENO" "$BASH_COMMAND"' ERR

    if [[ $probe_status -eq 126 || $probe_status -eq 127 ]] || \
       [[ "$probe_output" == *"cannot execute binary file"* ]] || \
       [[ "$probe_output" == *"Exec format error"* ]] || \
       [[ "$probe_output" == *"Bad CPU type in executable"* ]] || \
       [[ "$probe_output" == *"No such file or directory"* ]]; then
        write_error_log "PRECHECK" "liftOver preflight failed (os=${os}, arch=${arch}, file='${bin_desc:-unknown}', status=${probe_status}). Output: ${probe_output}"
        die "liftOver binary appears incompatible on this machine (OS=${os}, ARCH=${arch}). Use the matching UCSC build (e.g., macOSX.arm64 on Apple Silicon) or provide a working binary with -l."
    fi

    if [[ -n "$bin_desc" ]]; then
        if [[ "$os" == "Darwin" && "$arch" == "arm64" && "$bin_desc" == *"x86_64"* ]]; then
            warn "Using x86_64 liftOver on Apple Silicon. It may require Rosetta 2; prefer UCSC macOSX.arm64 binary."
        elif [[ "$os" == "Darwin" && "$arch" == "x86_64" && "$bin_desc" == *"arm64"* ]]; then
            warn "Using arm64 liftOver on Intel macOS; use UCSC macOSX.x86_64 binary."
        elif [[ "$os" == "Linux" && "$arch" == "aarch64" && "$bin_desc" == *"x86-64"* ]]; then
            warn "Using x86_64 liftOver on Linux aarch64; use UCSC linux.aarch64 binary."
        fi
    fi

    return 0
}

cleanup() {
    # Executes on script exit. Clears temporary directories and work logs
    # unless the user requested to keep intermediate files or if a crash occurred.
    if [[ -n "$WORK_DIR" && -d "$WORK_DIR" ]]; then
        if [[ "$RUN_FAILED" == true ]]; then
            purge_work_logs
        fi

        if [[ "$KEEP_INTERMEDIATE" == true ]]; then
            if [[ "$RUN_FAILED" == false ]]; then
                purge_work_logs
            fi
            info "Intermediate files kept: ${WORK_DIR}"
        else
            rm -rf "$WORK_DIR"
        fi
    fi
}

on_error() {
    # General error handler triggered via 'ERR' trap.
    # Collects the exit code, failing line, and failing command to log
    # comprehensive debugging information.
    local exit_code="${1:-$?}"
    local failed_line="${2:-unknown}"
    local failed_cmd="${3:-unknown}"
    local is_plink2_failure=false
    local failure_log_to_keep=""
    trap - ERR
    set +e
    if [[ "$exit_code" == "0" ]]; then
        set -e
        trap 'on_error "$?" "$LINENO" "$BASH_COMMAND"' ERR
        return
    fi
    RUN_FAILED=true
    if [[ "$failed_cmd" == *"plink2"* || -n "${PLINK2_ACTIVE_OUT_PREFIX:-}" ]]; then
        is_plink2_failure=true
        if [[ -n "${PLINK2_ACTIVE_OUT_PREFIX:-}" && -f "${PLINK2_ACTIVE_OUT_PREFIX}.log" ]]; then
            failure_log_to_keep="${PLINK2_ACTIVE_OUT_PREFIX}.log"
        elif [[ -n "${PLINK2_ACTIVE_STEP_LOG:-}" && -f "${PLINK2_ACTIVE_STEP_LOG}" ]]; then
            failure_log_to_keep="${PLINK2_ACTIVE_STEP_LOG}"
        elif [[ -n "${LAST_LOG_FILE:-}" && -f "$LAST_LOG_FILE" ]]; then
            failure_log_to_keep="$LAST_LOG_FILE"
        fi
        persist_failure_log "$failure_log_to_keep"
    fi
    write_error_log "UNHANDLED" "Command failed at line ${failed_line}: ${failed_cmd} (exit ${exit_code})"

    echo
    printf "%b[%s]%b %b%s%b Execution failed (exit code: %s)\n" "$DIM" "$(ts)" "$RESET" "$RED" "$GLYPH_ERR" "$RESET" "$exit_code" >&2
    [[ -n "$CURRENT_STEP" ]] && printf "%b[%s]%b %b%s%b Failed step: %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$CURRENT_STEP" >&2
    if [[ -n "${FAILED_STEP_LOG_FILE:-}" && -f "$FAILED_STEP_LOG_FILE" ]]; then
        printf "%b[%s]%b %b%s%b Retained failure log: %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$FAILED_STEP_LOG_FILE" >&2
    elif [[ -n "$LAST_LOG_FILE" && -f "$LAST_LOG_FILE" ]]; then
        printf "%b[%s]%b %b%s%b Failing step log (temporary): %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$LAST_LOG_FILE" >&2
    fi
    printf "%b[%s]%b %b%s%b Error log written to: %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$ERROR_LOG_FILE" >&2
    if [[ "$is_plink2_failure" == true ]]; then
        printf "%b[%s]%b %b%s%b PLINK2 docs: %s\n" "$DIM" "$(ts)" "$RESET" "$YELLOW" "$GLYPH_WARN" "$RESET" "$PLINK2_DOC_URL" >&2
    fi
    exit "$exit_code"
}

trap cleanup EXIT
trap 'on_error "$?" "$LINENO" "$BASH_COMMAND"' ERR

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
if [[ $# -eq 0 ]]; then
    usage 1
fi

# Handle supported long options before getopts.
PARSED_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --keep-all-chr)
            KEEP_ALL_CHR=true
            shift
            ;;
        --force)
            FORCE_LIFTOVER=true
            shift
            ;;
        --force-all)
            FORCE_ALL_LIFTOVER=true
            shift
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                PARSED_ARGS+=("$1")
                shift
            done
            ;;
        --*)
            die "Unknown option: $1"
            ;;
        *)
            PARSED_ARGS+=("$1")
            shift
            ;;
    esac
done
set -- "${PARSED_ARGS[@]}"

while getopts ":i:I:c:l:o:O:f:t:v:kh" opt; do
    case "$opt" in
        i) INPUT_PREFIX="$OPTARG" ;;
        I) INPUT_VCF="$OPTARG" ;;
        c) CHAIN_FILE="$OPTARG" ;;
        l) LIFTOVER_BIN="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        O) OUTPUT_VCF="$OPTARG" ;;
        f) INPUT_FORMAT="$OPTARG" ;;
        t) OUTPUT_FORMAT="$OPTARG" ;;
        v) VARID_TEMPLATE="$OPTARG" ;;
        k) KEEP_INTERMEDIATE=true ;;
        h) usage 0 ;;
        :) die "Option -${OPTARG} requires an argument." ;;
        \?) usage 1 ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate arguments and dependencies
# -----------------------------------------------------------------------------
is_vcf_path() {
    local p="$1"
    [[ "$p" == *.vcf || "$p" == *.vcf.gz ]]
}

is_plink_component_path() {
    local p="$1"
    [[ "$p" == *.bed || "$p" == *.bim || "$p" == *.fam || "$p" == *.pgen || "$p" == *.pvar || "$p" == *.psam ]]
}

normalize_plink_prefix_arg() {
    local p="$1"
    if is_plink_component_path "$p"; then
        printf '%s\n' "${p%.*}"
        return 0
    fi
    printf '%s\n' "$p"
}

filter_primary_target_contigs() {
    # Removes any variants mapped to non-primary contigs (e.g., alternative loci or
    # unplaced contigs) from the lifted BED file, keeping only standard chromosomes (1-22, X, Y).
    local in_bed="$1"
    local out_bed="$2"
    awk '
        function canon_chr(c, x) {
            x = toupper(c)
            sub(/^CHR/, "", x)
            if (x == "23") x = "X"
            else if (x == "24") x = "Y"
            else if (x == "M") x = "MT"
            return x
        }
        function is_primary(c) {
            return (c ~ /^([1-9]|1[0-9]|2[0-2])$/ || c == "X" || c == "Y")
        }
        {
            c = canon_chr($1)
            if (is_primary(c)) print
        }
    ' "$in_bed" > "$out_bed"
}

detect_vcf_build() {
    # Attempts to auto-detect the human genome build (GRCh37 vs GRCh38) from a VCF file.
    # It analyzes two streams of evidence:
    # 1. Header Analysis: Checks for contig lengths and known sequence tokens (e.g., NC_000001.10).
    # 2. Coordinate Analysis: Uses terminal positional coordinates on primary chromosomes 
    #    to determine if the maximum variant coordinate exceeds the known length of GRCh37 or GRCh38.
    local input_vcf="$1"
    local header_out="$2"
    local log_out="$3"
    local expected_build="$4"
    local header_result header_build header_reason
    local pos_result pos_build pos_reason
    local pos_primary_count pos_primary_list

    # Extract header lines and check for exact reference genome markers and specific contig lengths
    bcftools view -h "$input_vcf" > "$header_out" 2> "$log_out"

    header_result="$(awk '
        BEGIN {
            strong37 = 0
            strong38 = 0
            token37 = 0
            token38 = 0
            set_len("1", 249250621, 248956422)
            set_len("2", 243199373, 242193529)
            set_len("3", 198022430, 198295559)
            set_len("4", 191154276, 190214555)
            set_len("5", 180915260, 181538259)
            set_len("6", 171115067, 170805979)
            set_len("7", 159138663, 159345973)
            set_len("8", 146364022, 145138636)
            set_len("9", 141213431, 138394717)
            set_len("10", 135534747, 133797422)
            set_len("11", 135006516, 135086622)
            set_len("12", 133851895, 133275309)
            set_len("13", 115169878, 114364328)
            set_len("14", 107349540, 107043718)
            set_len("15", 102531392, 101991189)
            set_len("16", 90354753, 90338345)
            set_len("17", 81195210, 83257441)
            set_len("18", 78077248, 80373285)
            set_len("19", 59128983, 58617616)
            set_len("20", 63025520, 64444167)
            set_len("21", 48129895, 46709983)
            set_len("22", 51304566, 50818468)
            set_len("X", 155270560, 156040895)
            set_len("Y", 59373566, 57227415)
        }
        function set_len(chr, l37, l38) {
            len37[chr] = l37
            len38[chr] = l38
        }
        function canon_chr(c, x) {
            x = toupper(c)
            sub(/^CHR/, "", x)
            if (x == "M") x = "MT"
            return x
        }
        function score_contig(id, len, cid) {
            cid = canon_chr(id)
            if (!(cid in len37) || !(cid in len38)) return
            if (len == len37[cid]) strong37++
            else if (len == len38[cid]) strong38++
        }
        /^##contig=<.*>/ {
            id = ""
            len = ""
            if (match($0, /ID=[^,>]+/)) id = substr($0, RSTART + 3, RLENGTH - 3)
            if (match($0, /length=[0-9]+/)) len = substr($0, RSTART + 7, RLENGTH - 7)
            if (id != "" && len != "") score_contig(id, len)
        }
        {
            lower = tolower($0)
            if (lower ~ /grch37|hg19|b37/) token37++
            if (lower ~ /grch38|hg38/) token38++
            if ($0 ~ /NC_000001\.10|CM000663\.1/) strong37++
            if ($0 ~ /NC_000001\.11|CM000663\.2/) strong38++
        }
        END {
            build = "unknown"
            reason = "No GRCh37/GRCh38 markers found in VCF header"
            if (strong37 > 0 && strong38 > 0) {
                reason = "Conflicting GRCh37/GRCh38 contig markers in VCF header"
            } else if (strong37 > 0) {
                build = "GRCh37"
                reason = "Detected GRCh37/hg19 contig markers in VCF header"
            } else if (strong38 > 0) {
                build = "GRCh38"
                reason = "Detected GRCh38/hg38 contig markers in VCF header"
            } else if (token37 > 0 && token38 > 0) {
                reason = "Conflicting GRCh37/GRCh38 tokens in VCF header"
            } else if (token37 > 0) {
                build = "GRCh37"
                reason = "Detected GRCh37/hg19 tokens in VCF header"
            } else if (token38 > 0) {
                build = "GRCh38"
                reason = "Detected GRCh38/hg38 tokens in VCF header"
            }
            print build "\t" reason
        }
    ' "$header_out")"

    header_build="${header_result%%$'\t'*}"
    header_reason="${header_result#*$'\t'}"

    # Extract all positions from VCF to scan for terminal positions.
    # We find the largest coordinate (max_pos) for each primary chromosome.
    # These max coordinates are compared against standard GRCh37/GRCh38 chromosome lengths.
    pos_result="$(bcftools query -f '%CHROM\t%POS\n' "$input_vcf" 2>>"$log_out" | awk -v expected_build="$expected_build" '
        BEGIN {
            set_len("1", 249250621, 248956422)
            set_len("2", 243199373, 242193529)
            set_len("3", 198022430, 198295559)
            set_len("4", 191154276, 190214555)
            set_len("5", 180915260, 181538259)
            set_len("6", 171115067, 170805979)
            set_len("7", 159138663, 159345973)
            set_len("8", 146364022, 145138636)
            set_len("9", 141213431, 138394717)
            set_len("10", 135534747, 133797422)
            set_len("11", 135006516, 135086622)
            set_len("12", 133851895, 133275309)
            set_len("13", 115169878, 114364328)
            set_len("14", 107349540, 107043718)
            set_len("15", 102531392, 101991189)
            set_len("16", 90354753, 90338345)
            set_len("17", 81195210, 83257441)
            set_len("18", 78077248, 80373285)
            set_len("19", 59128983, 58617616)
            set_len("20", 63025520, 64444167)
            set_len("21", 48129895, 46709983)
            set_len("22", 51304566, 50818468)
            set_len("X", 155270560, 156040895)
            set_len("Y", 59373566, 57227415)
        }
        function set_len(chr, l37, l38) {
            len37[chr] = l37
            len38[chr] = l38
        }
        function canon_chr(c, x, ux) {
            x = c
            sub(/^chr/, "", x)
            sub(/^CHR/, "", x)
            ux = toupper(x)
            if (ux == "23") x = "X"
            else if (ux == "24") x = "Y"
            else if (ux == "M" || ux == "MT" || ux == "26") x = "MT"
            else if (ux == "X" || ux == "Y" || ux ~ /^([1-9]|1[0-9]|2[0-2])$/) x = ux
            else x = ux
            return x
        }
        function is_primary(c, uc) {
            uc = toupper(c)
            return (uc ~ /^([1-9]|1[0-9]|2[0-2])$/ || uc == "X" || uc == "Y")
        }
        function remember_primary(chr, c) {
            c = canon_chr(chr)
            if (!is_primary(c)) return
            if (!(c in primary_seen)) {
                primary_seen[c] = 1
                primary_count++
                primary_order[++primary_n] = c
            }
        }
        function capture(chr, pos, c, p) {
            if (pos !~ /^[0-9]+$/) return
            c = canon_chr(chr)
            p = pos + 0
            if (is_primary(c) && p > max_pos[c]) {
                max_pos[c] = p
            }
        }
        {
            remember_primary($1)
            capture($1, $2)
        }
        END {
            evidence37 = 0
            evidence38 = 0
            invalid = 0
            detail = ""
            primary_list = ""
            for (i = 1; i <= primary_n; i++) {
                c = primary_order[i]
                p = max_pos[c] + 0
                if (p <= 0) continue

                l37 = len37[c] + 0
                l38 = len38[c] + 0
                if (l37 <= 0 || l38 <= 0) continue

                max_len = (l37 > l38 ? l37 : l38)
                min_len = (l37 < l38 ? l37 : l38)

                if (p > max_len) {
                    invalid++
                    detail = detail "chr" c " exceeds GRCh37/GRCh38 lengths; "
                    continue
                }
                if (l37 > l38 && p > min_len) {
                    evidence37++
                    detail = detail "chr" c " tail supports GRCh37; "
                } else if (l38 > l37 && p > min_len) {
                    evidence38++
                    detail = detail "chr" c " tail supports GRCh38; "
                } else {
                    diff = max_len - min_len
                    if (diff >= 50000) {
                        dist37 = l37 - p
                        dist38 = l38 - p
                        if (dist37 >= 0 && dist37 <= 500000 && dist38 - dist37 >= 50000) {
                            evidence37++
                            detail = detail "chr" c " proximity to GRCh37 end; "
                        } else if (dist38 >= 0 && dist38 <= 500000 && dist37 - dist38 >= 50000) {
                            evidence38++
                            detail = detail "chr" c " proximity to GRCh38 end; "
                        }
                    }
                }
            }

            for (i = 1; i <= primary_n; i++) {
                if (i > 1) primary_list = primary_list ","
                primary_list = primary_list primary_order[i]
            }
            if (primary_list == "") primary_list = "none"

            if (invalid > 0) {
                print "unknown\tCoordinates exceed GRCh37/GRCh38 primary chromosome lengths (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence37 > 0 && evidence38 > 0) {
                print "unknown\tConflicting GRCh37/GRCh38 positional evidence (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence37 > 0) {
                if (evidence37 > 1) detail = evidence37 " chromosomes support GRCh37"
                else sub(/; $/, "", detail)
                print "GRCh37\tDetected GRCh37 positional evidence from chromosome-tail coordinates (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence38 > 0) {
                if (evidence38 > 1) detail = evidence38 " chromosomes support GRCh38"
                else sub(/; $/, "", detail)
                print "GRCh38\tDetected GRCh38 positional evidence from chromosome-tail coordinates (" detail ")\t" primary_count "\t" primary_list
            } else {
                if (expected_build != "unknown" && expected_build != "") {
                    print expected_build "\tCoordinates fit within both GRCh37 and GRCh38; assuming expected source build (" expected_build ")\t" primary_count "\t" primary_list
                } else {
                    print "unknown\tNo build-discriminating primary chromosome-tail coordinates found\t" primary_count "\t" primary_list
                }
            }
        }
    ')"

    IFS=$'\t' read -r pos_build pos_reason pos_primary_count pos_primary_list <<< "$pos_result"
    [[ -n "$pos_build" ]] || pos_build="unknown"
    [[ -n "$pos_reason" ]] || pos_reason="Unable to parse coordinate evidence"
    [[ "$pos_primary_count" =~ ^[0-9]+$ ]] || pos_primary_count=0
    [[ -n "$pos_primary_list" ]] || pos_primary_list="none"

    if [[ "$header_build" != "unknown" && "$pos_build" != "unknown" && "$header_build" != "$pos_build" ]]; then
        printf "unknown\tConflicting build evidence between VCF header and variant coordinates (header: %s; coordinates: %s)\t%s\t%s\n" "$header_reason" "$pos_reason" "$pos_primary_count" "$pos_primary_list"
    elif [[ "$header_build" != "unknown" ]]; then
        if [[ "$pos_build" == "unknown" ]]; then
            printf "%s\t%s; coordinate-check inconclusive: %s\t%s\t%s\n" "$header_build" "$header_reason" "$pos_reason" "$pos_primary_count" "$pos_primary_list"
        else
            printf "%s\t%s; coordinate-check agrees: %s\t%s\t%s\n" "$header_build" "$header_reason" "$pos_reason" "$pos_primary_count" "$pos_primary_list"
        fi
    elif [[ "$pos_build" != "unknown" ]]; then
        printf "%s\t%s; header-check inconclusive: %s\t%s\t%s\n" "$pos_build" "$pos_reason" "$header_reason" "$pos_primary_count" "$pos_primary_list"
    else
        printf "unknown\tCould not infer from VCF header or variant coordinates (%s; %s)\t%s\t%s\n" "$header_reason" "$pos_reason" "$pos_primary_count" "$pos_primary_list"
    fi
}

detect_plink_build() {
    # Attempts to auto-detect the human genome build (GRCh37 vs GRCh38) from PLINK variant files.
    # Since PLINK files typically lack detailed headers, this relies solely on
    # analyzing terminal variant coordinates against known chromosome lengths.
    local input_format="$1"
    local variant_file="$2"
    local expected_build="$3"

    awk -v fmt="$input_format" -v expected_build="$expected_build" '
        BEGIN {
            set_len("1", 249250621, 248956422)
            set_len("2", 243199373, 242193529)
            set_len("3", 198022430, 198295559)
            set_len("4", 191154276, 190214555)
            set_len("5", 180915260, 181538259)
            set_len("6", 171115067, 170805979)
            set_len("7", 159138663, 159345973)
            set_len("8", 146364022, 145138636)
            set_len("9", 141213431, 138394717)
            set_len("10", 135534747, 133797422)
            set_len("11", 135006516, 135086622)
            set_len("12", 133851895, 133275309)
            set_len("13", 115169878, 114364328)
            set_len("14", 107349540, 107043718)
            set_len("15", 102531392, 101991189)
            set_len("16", 90354753, 90338345)
            set_len("17", 81195210, 83257441)
            set_len("18", 78077248, 80373285)
            set_len("19", 59128983, 58617616)
            set_len("20", 63025520, 64444167)
            set_len("21", 48129895, 46709983)
            set_len("22", 51304566, 50818468)
            set_len("X", 155270560, 156040895)
            set_len("Y", 59373566, 57227415)
        }
        function set_len(chr, l37, l38) {
            len37[chr] = l37
            len38[chr] = l38
        }
        function canon_chr(c, x) {
            x = toupper(c)
            sub(/^CHR/, "", x)
            if (x == "23") x = "X"
            else if (x == "24") x = "Y"
            else if (x == "M") x = "MT"
            return x
        }
        function is_primary(c, uc) {
            uc = toupper(c)
            return (uc ~ /^([1-9]|1[0-9]|2[0-2])$/ || uc == "X" || uc == "Y")
        }
        function remember_primary(chr, c) {
            c = canon_chr(chr)
            if (!is_primary(c)) return
            if (!(c in primary_seen)) {
                primary_seen[c] = 1
                primary_count++
                primary_order[++primary_n] = c
            }
        }
        function capture(chr, pos, c, p) {
            if (pos !~ /^[0-9]+$/) return
            c = canon_chr(chr)
            p = pos + 0
            if (is_primary(c) && p > max_pos[c]) {
                max_pos[c] = p
            }
        }
        {
            if (fmt == "plink2") {
                if ($0 ~ /^#/) next
                remember_primary($1)
                capture($1, $2)
            } else {
                remember_primary($1)
                capture($1, $4)
            }
        }
        END {
            evidence37 = 0
            evidence38 = 0
            invalid = 0
            detail = ""
            primary_list = ""
            for (i = 1; i <= primary_n; i++) {
                c = primary_order[i]
                p = max_pos[c] + 0
                if (p <= 0) continue

                l37 = len37[c] + 0
                l38 = len38[c] + 0
                if (l37 <= 0 || l38 <= 0) continue

                max_len = (l37 > l38 ? l37 : l38)
                min_len = (l37 < l38 ? l37 : l38)

                if (p > max_len) {
                    invalid++
                    detail = detail "chr" c " exceeds GRCh37/GRCh38 lengths; "
                    continue
                }
                if (l37 > l38 && p > min_len) {
                    evidence37++
                    detail = detail "chr" c " tail supports GRCh37; "
                } else if (l38 > l37 && p > min_len) {
                    evidence38++
                    detail = detail "chr" c " tail supports GRCh38; "
                } else {
                    diff = max_len - min_len
                    if (diff >= 50000) {
                        dist37 = l37 - p
                        dist38 = l38 - p
                        if (dist37 >= 0 && dist37 <= 500000 && dist38 - dist37 >= 50000) {
                            evidence37++
                            detail = detail "chr" c " proximity to GRCh37 end; "
                        } else if (dist38 >= 0 && dist38 <= 500000 && dist37 - dist38 >= 50000) {
                            evidence38++
                            detail = detail "chr" c " proximity to GRCh38 end; "
                        }
                    }
                }
            }

            for (i = 1; i <= primary_n; i++) {
                if (i > 1) primary_list = primary_list ","
                primary_list = primary_list primary_order[i]
            }
            if (primary_list == "") primary_list = "none"

            if (invalid > 0) {
                print "unknown\tCoordinates exceed GRCh37/GRCh38 primary chromosome lengths (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence37 > 0 && evidence38 > 0) {
                print "unknown\tConflicting GRCh37/GRCh38 positional evidence (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence37 > 0) {
                if (evidence37 > 1) detail = evidence37 " chromosomes support GRCh37"
                else sub(/; $/, "", detail)
                print "GRCh37\tDetected GRCh37 positional evidence from chromosome-tail coordinates (" detail ")\t" primary_count "\t" primary_list
            } else if (evidence38 > 0) {
                if (evidence38 > 1) detail = evidence38 " chromosomes support GRCh38"
                else sub(/; $/, "", detail)
                print "GRCh38\tDetected GRCh38 positional evidence from chromosome-tail coordinates (" detail ")\t" primary_count "\t" primary_list
            } else {
                if (expected_build != "unknown" && expected_build != "") {
                    print expected_build "\tCoordinates fit within both GRCh37 and GRCh38; assuming expected source build (" expected_build ")\t" primary_count "\t" primary_list
                } else {
                    print "unknown\tNo build-discriminating primary chromosome-tail coordinates found\t" primary_count "\t" primary_list
                }
            }
        }
    ' "$variant_file"
}

HAS_PLINK_FLAGS=false
HAS_VCF_FLAGS=false

# Validate that the user did not mix VCF and PLINK CLI arguments
if [[ -n "$INPUT_PREFIX" || -n "$OUTPUT_PREFIX" ]]; then
    HAS_PLINK_FLAGS=true
fi
if [[ -n "$INPUT_VCF" || -n "$OUTPUT_VCF" ]]; then
    HAS_VCF_FLAGS=true
fi

if [[ "$HAS_PLINK_FLAGS" == true && "$HAS_VCF_FLAGS" == true ]]; then
    die "Do not mix PLINK flags (-i/-o) with VCF flags (-I/-O)."
elif [[ "$HAS_VCF_FLAGS" == true ]]; then
    MODE="vcf"
elif [[ "$HAS_PLINK_FLAGS" == true ]]; then
    MODE="plink"
else
    die "Missing input/output flags. Use either -i/-o (PLINK) or -I/-O (VCF)."
fi

[[ -z "$CHAIN_FILE" ]] && die "Missing required argument: -c"

if [[ "$MODE" == "plink" ]]; then
    [[ -z "$INPUT_PREFIX" ]] && die "Missing required argument: -i"
    [[ -z "$OUTPUT_PREFIX" ]] && die "Missing required argument: -o"
    [[ -z "$INPUT_VCF" ]] || die "Do not use -I with PLINK mode. Use -i."
    [[ -z "$OUTPUT_VCF" ]] || die "Do not use -O with PLINK mode. Use -o."

    is_vcf_path "$INPUT_PREFIX" && die "Input '$INPUT_PREFIX' looks like VCF. Use VCF flags -I/-O instead."
    is_vcf_path "$OUTPUT_PREFIX" && die "Output '$OUTPUT_PREFIX' looks like VCF. Use VCF flags -I/-O instead."

    INPUT_PREFIX="$(normalize_plink_prefix_arg "$INPUT_PREFIX")"
    OUTPUT_PREFIX="$(normalize_plink_prefix_arg "$OUTPUT_PREFIX")"
else
    [[ -z "$INPUT_PREFIX" ]] || die "Do not use -i with VCF mode. Use -I."
    [[ -z "$OUTPUT_PREFIX" ]] || die "Do not use -o with VCF mode. Use -O."
    [[ -z "$INPUT_VCF" ]] && die "Missing required argument: -I (input VCF path)"
    [[ -z "$OUTPUT_VCF" ]] && die "Missing required argument: -O (output VCF path)"
    is_vcf_path "$INPUT_VCF" || die "Input VCF must end with .vcf or .vcf.gz: $INPUT_VCF"
    is_vcf_path "$OUTPUT_VCF" || die "Output VCF must end with .vcf or .vcf.gz: $OUTPUT_VCF"
    [[ -f "$INPUT_VCF" ]] || die "Input VCF not found: $INPUT_VCF"

    OUTPUT_PREFIX="$OUTPUT_VCF"
    OUTPUT_PREFIX="${OUTPUT_PREFIX%.vcf.gz}"
    OUTPUT_PREFIX="${OUTPUT_PREFIX%.vcf}"
fi

# -- Directory Organization --
BASE_OUT_DIR="$OUTPUT_PREFIX"
BASE_NAME="$(basename "$OUTPUT_PREFIX")"
LOGS_DIR="${BASE_OUT_DIR}/logs"
LIFTED_DIR="${BASE_OUT_DIR}/lifted-data"
mkdir -p "$LOGS_DIR" "$LIFTED_DIR"

OUTPUT_PREFIX="${LIFTED_DIR}/${BASE_NAME}"
LOG_PREFIX="${LOGS_DIR}/LocusLift"

if [[ "$MODE" == "vcf" ]]; then
    if [[ "$OUTPUT_VCF" == *.vcf.gz ]]; then
        OUTPUT_VCF="${OUTPUT_PREFIX}.vcf.gz"
    else
        OUTPUT_VCF="${OUTPUT_PREFIX}.vcf"
    fi
fi

init_error_log_path
init_failure_log_path
init_run_summary_path
rm -f "$ERROR_LOG_FILE" "$FAILED_STEP_LOG_FILE" "$RUN_SUMMARY_FILE"

if [[ "$MODE" == "plink" && "$INPUT_PREFIX" == "$OUTPUT_PREFIX" ]]; then
    die "Input and output prefixes must differ to avoid in-place overwrite."
fi

[[ -f "$CHAIN_FILE" ]] || die "Chain file not found: $CHAIN_FILE"

if [[ -n "$LIFTOVER_BIN" ]]; then
    # If user passed a bare command name (no slash), resolve via PATH.
    if [[ "$LIFTOVER_BIN" != */* ]]; then
        if command -v "$LIFTOVER_BIN" >/dev/null 2>&1; then
            LIFTOVER_BIN="$(command -v "$LIFTOVER_BIN")"
        elif [[ -x "./${LIFTOVER_BIN}" ]]; then
            LIFTOVER_BIN="./${LIFTOVER_BIN}"
        elif [[ -x "${SCRIPT_DIR}/${LIFTOVER_BIN}" ]]; then
            LIFTOVER_BIN="${SCRIPT_DIR}/${LIFTOVER_BIN}"
        else
            die "Provided -l '${LIFTOVER_BIN}' was not found as a command or executable file."
        fi
    fi
elif [[ -z "$LIFTOVER_BIN" ]]; then
    if [[ -x "./liftOver" ]]; then
        LIFTOVER_BIN="./liftOver"
    elif [[ -x "${SCRIPT_DIR}/liftOver" ]]; then
        LIFTOVER_BIN="${SCRIPT_DIR}/liftOver"
    elif command -v liftOver >/dev/null 2>&1; then
        LIFTOVER_BIN="$(command -v liftOver)"
    elif [[ -f "./liftOver" || -f "${SCRIPT_DIR}/liftOver" ]]; then
        die "Found liftOver in current/script directory but it is not executable. Run: chmod +x ./liftOver"
    else
        die "liftOver binary not provided (-l) and not found in current directory, script directory, or PATH."
    fi
fi

[[ -x "$LIFTOVER_BIN" ]] || die "liftOver binary is not executable: $LIFTOVER_BIN"
command -v awk >/dev/null 2>&1 || die "awk not found in PATH"
if [[ "$MODE" == "plink" ]]; then
    command -v plink2 >/dev/null 2>&1 || die "plink2 not found in PATH"
else
    command -v bcftools >/dev/null 2>&1 || die "bcftools not found in PATH (required for VCF mode)"
fi
preflight_liftover_binary

if [[ "$MODE" == "plink" ]]; then
    # -----------------------------------------------------------------------------
    # Input format detection (PLINK mode)
    # Automatically infer PLINK 1 vs PLINK 2 by checking file extensions.
    # -----------------------------------------------------------------------------
    detect_input_format() {
        if [[ -f "${INPUT_PREFIX}.pgen" && -f "${INPUT_PREFIX}.pvar" && -f "${INPUT_PREFIX}.psam" ]]; then
            echo "plink2"
        elif [[ -f "${INPUT_PREFIX}.bed" && -f "${INPUT_PREFIX}.bim" && -f "${INPUT_PREFIX}.fam" ]]; then
            echo "plink1"
        else
            echo ""
        fi
    }

    if [[ -z "$INPUT_FORMAT" ]]; then
        INPUT_FORMAT="$(detect_input_format)"
        [[ -n "$INPUT_FORMAT" ]] || die "Could not auto-detect input format for prefix: ${INPUT_PREFIX}"
    fi

    case "$INPUT_FORMAT" in
        plink1)
            [[ -f "${INPUT_PREFIX}.bed" ]] || die "Missing file: ${INPUT_PREFIX}.bed"
            [[ -f "${INPUT_PREFIX}.bim" ]] || die "Missing file: ${INPUT_PREFIX}.bim"
            [[ -f "${INPUT_PREFIX}.fam" ]] || die "Missing file: ${INPUT_PREFIX}.fam"
            VARIANT_FILE="${INPUT_PREFIX}.bim"
            PLINK_INPUT_ARGS=(--bfile "$INPUT_PREFIX")
            ;;
        plink2)
            [[ -f "${INPUT_PREFIX}.pgen" ]] || die "Missing file: ${INPUT_PREFIX}.pgen"
            [[ -f "${INPUT_PREFIX}.pvar" ]] || die "Missing file: ${INPUT_PREFIX}.pvar"
            [[ -f "${INPUT_PREFIX}.psam" ]] || die "Missing file: ${INPUT_PREFIX}.psam"
            VARIANT_FILE="${INPUT_PREFIX}.pvar"
            PLINK_INPUT_ARGS=(--pfile "$INPUT_PREFIX")
            ;;
        *)
            die "Invalid input format: ${INPUT_FORMAT}. Use 'plink1' or 'plink2'."
            ;;
    esac

    if [[ -z "$OUTPUT_FORMAT" ]]; then
        OUTPUT_FORMAT="$INPUT_FORMAT"
    fi

    case "$OUTPUT_FORMAT" in
        plink1) PLINK_OUTPUT_ARGS=(--make-bed) ;;
        plink2) PLINK_OUTPUT_ARGS=(--make-pgen) ;;
        *) die "Invalid output format: ${OUTPUT_FORMAT}. Use 'plink1' or 'plink2'." ;;
    esac

    build_ucsc_bed() {
        # Converts PLINK 1 (.bim) or PLINK 2 (.pvar) coordinate data to UCSC BED format (0-based start, 1-based end).
        # This standardized format is required by the liftOver binary.
        local in_format="$1"
        local variant_file="$2"
        local output_bed="$3"

        if [[ "$in_format" == "plink1" ]]; then
            # .bim columns: CHR ID CM POS A1 A2
            awk -v OFS='\t' -v keep_all="$KEEP_ALL_CHR" '
                function canon_chr(c, x) {
                    x = c
                    sub(/^chr/, "", x)
                    sub(/^CHR/, "", x)
                    ux = toupper(x)
                    if (ux == "23") x = "X"
                    else if (ux == "24") x = "Y"
                    else if (ux == "26" || ux == "M" || ux == "MT") x = "MT"
                    else if (ux == "25" || ux == "XY") x = "XY"
                    else if (ux == "X" || ux == "Y" || ux ~ /^([1-9]|1[0-9]|2[0-2])$/) x = ux
                    return x
                }
                function is_primary(c) {
                    uc = toupper(c)
                    return (uc ~ /^([1-9]|1[0-9]|2[0-2])$/ || uc == "X" || uc == "Y")
                }
                function norm_chr(c, x) {
                    x = canon_chr(c)
                    if (x == "MT") return "chrM"
                    return "chr" x
                }
                $4 ~ /^[0-9]+$/ {
                    cc = canon_chr($1)
                    if (keep_all != "true" && !is_primary(cc)) next
                    chr = norm_chr($1)
                    start = ($4 > 0 ? $4 - 1 : 0)
                    print chr, start, $4, $2
                }
            ' "$variant_file" > "$output_bed"
        else
            # .pvar columns: #CHROM POS ID REF ALT ... (skip header lines)
            awk -v OFS='\t' -v keep_all="$KEEP_ALL_CHR" '
                function canon_chr(c, x) {
                    x = c
                    sub(/^chr/, "", x)
                    sub(/^CHR/, "", x)
                    ux = toupper(x)
                    if (ux == "23") x = "X"
                    else if (ux == "24") x = "Y"
                    else if (ux == "26" || ux == "M" || ux == "MT") x = "MT"
                    else if (ux == "25" || ux == "XY") x = "XY"
                    else if (ux == "X" || ux == "Y" || ux ~ /^([1-9]|1[0-9]|2[0-2])$/) x = ux
                    return x
                }
                function is_primary(c) {
                    uc = toupper(c)
                    return (uc ~ /^([1-9]|1[0-9]|2[0-2])$/ || uc == "X" || uc == "Y")
                }
                function norm_chr(c, x) {
                    x = canon_chr(c)
                    if (x == "MT") return "chrM"
                    return "chr" x
                }
                /^#/ { next }
                $2 ~ /^[0-9]+$/ {
                    cc = canon_chr($1)
                    if (keep_all != "true" && !is_primary(cc)) next
                    chr = norm_chr($1)
                    start = ($2 > 0 ? $2 - 1 : 0)
                    print chr, start, $2, $3
                }
            ' "$variant_file" > "$output_bed"
        fi
    }
else
    INPUT_FORMAT="vcf"
    OUTPUT_FORMAT="vcf"
fi

# -----------------------------------------------------------------------------
# Infer direction from chain filename; select sensible var-id template default
# e.g., 'hg19ToHg38.over.chain.gz' sets SOURCE=GRCh37 and TARGET=GRCh38.
# -----------------------------------------------------------------------------
CHAIN_BASENAME="$(basename "$CHAIN_FILE")"
CHAIN_LOWER="$(printf '%s' "$CHAIN_BASENAME" | tr '[:upper:]' '[:lower:]')"

SOURCE_BUILD="source"
TARGET_BUILD="target"
CHAIN_BUILD_GUARD_ENABLED=false

if [[ "$CHAIN_LOWER" =~ (hg38|grch38).*(hg19|grch37|b37) ]]; then
    SOURCE_BUILD="GRCh38"
    TARGET_BUILD="GRCh37"
    CHAIN_BUILD_GUARD_ENABLED=true
    if [[ "$MODE" == "plink" ]]; then
        [[ -z "$VARID_TEMPLATE" ]] && VARID_TEMPLATE='@:#:$r:$a'
    fi
elif [[ "$CHAIN_LOWER" =~ (hg19|grch37|b37).*(hg38|grch38) ]]; then
    SOURCE_BUILD="GRCh37"
    TARGET_BUILD="GRCh38"
    CHAIN_BUILD_GUARD_ENABLED=true
    if [[ "$MODE" == "plink" ]]; then
        [[ -z "$VARID_TEMPLATE" ]] && VARID_TEMPLATE='chr@:#:$r:$a'
    fi
else
    warn "Could not infer direction from chain filename; using generic labels."
    if [[ "$MODE" == "plink" ]]; then
        [[ -z "$VARID_TEMPLATE" ]] && VARID_TEMPLATE='@:#:$r:$a'
    fi
fi

# -----------------------------------------------------------------------------
# Prepare workspace
# -----------------------------------------------------------------------------
OUTPUT_DIR="${BASE_OUT_DIR:-$(dirname "$OUTPUT_PREFIX")}"
mkdir -p "$OUTPUT_DIR"
WORK_DIR="$(mktemp -d "${OUTPUT_DIR}/.plink-liftover.XXXXXX")"

# -----------------------------------------------------------------------------
# Header
# -----------------------------------------------------------------------------
echo
line
printf "%bLocusLift%b %s %b- Genome LiftOver Pipeline (PLINK/VCF)%b\n" "$BOLD" "$RESET" "$DIM(v${VERSION})$RESET" "$DIM" "$RESET"
line
if [[ "$MODE" == "plink" ]]; then
    print_kv "Input prefix" "$INPUT_PREFIX"
    print_kv "Input format" "$INPUT_FORMAT"
    print_kv "Output prefix" "$OUTPUT_PREFIX"
    print_kv "Output format" "$OUTPUT_FORMAT"
    print_kv "Var-ID template" "$VARID_TEMPLATE"
else
    print_kv "Input VCF" "$INPUT_VCF"
    print_kv "Input format" "$INPUT_FORMAT"
    print_kv "Output VCF" "$OUTPUT_VCF"
    print_kv "Output format" "$OUTPUT_FORMAT"
fi
print_kv "Chain file" "$CHAIN_FILE"
print_kv "Direction" "${SOURCE_BUILD} -> ${TARGET_BUILD}"
print_kv "liftOver bin" "$LIFTOVER_BIN"
print_kv "Temp dir" "$WORK_DIR"
line

if [[ "$KEEP_ALL_CHR" == false ]]; then
    warn "Non-primary variants (XY/MT etc.) are dropped by default for reliability. Use --keep-all-chr to retain them."
fi
if [[ "$FORCE_LIFTOVER" == true ]]; then
    warn "Force mode enabled (--force): liftover will continue when input build cannot be confidently inferred (at your own risk)."
fi
if [[ "$FORCE_ALL_LIFTOVER" == true ]]; then
    warn "Force all mode enabled (--force-all): liftover will bypass ALL safety checks, including proven mismatches (AT YOUR OWN RISK)."
fi

if [[ "$MODE" == "vcf" ]]; then
    # =============================================================================
    # VCF MODE PIPELINE
    # 1) Assembly Check -> 2) Query VCF -> 3) liftOver -> 4) Coordinate Map -> 
    # 5) Rewrite VCF body -> 6) Sort and write
    # =============================================================================

    # -----------------------------------------------------------------------------
    # STEP 1: Detect assembly build and validate chain direction
    # Ensures that the input VCF actually matches the source build specified by the chain file.
    # -----------------------------------------------------------------------------
    CURRENT_STEP="1/6 assembly check"
    step "1/6 Detect input assembly and validate chain direction"

    VCF_HEADER_CHECK="${WORK_DIR}/vcf-header-check.txt"
    STEP1_LOG="${WORK_DIR}/vcf-step1-assembly-check.log"
    LAST_LOG_FILE="$STEP1_LOG"

    IFS=$'\t' read -r DETECTED_INPUT_BUILD DETECTED_INPUT_BUILD_REASON INPUT_PRIMARY_CHR_COUNT INPUT_PRIMARY_CHR_LIST \
        < <(detect_vcf_build "$INPUT_VCF" "$VCF_HEADER_CHECK" "$STEP1_LOG" "$SOURCE_BUILD")

    [[ -n "$DETECTED_INPUT_BUILD" ]] || DETECTED_INPUT_BUILD="unknown"
    [[ -n "$DETECTED_INPUT_BUILD_REASON" ]] || DETECTED_INPUT_BUILD_REASON="Unable to parse build markers from VCF header"
    [[ "$INPUT_PRIMARY_CHR_COUNT" =~ ^[0-9]+$ ]] || INPUT_PRIMARY_CHR_COUNT=0
    [[ -n "$INPUT_PRIMARY_CHR_LIST" ]] || INPUT_PRIMARY_CHR_LIST="none"

    {
        echo "Detected build: ${DETECTED_INPUT_BUILD}"
        echo "Build evidence: ${DETECTED_INPUT_BUILD_REASON}"
        echo "Observed primary chromosome count: ${INPUT_PRIMARY_CHR_COUNT}"
        echo "Observed primary chromosomes: ${INPUT_PRIMARY_CHR_LIST}"
        echo "Chain source build: ${SOURCE_BUILD}"
        echo "Chain target build: ${TARGET_BUILD}"
    } >> "$STEP1_LOG"

    if [[ "$CHAIN_BUILD_GUARD_ENABLED" == true ]]; then
        if [[ "$DETECTED_INPUT_BUILD" == "unknown" ]]; then
            CHAIN_COMPATIBILITY_STATUS="FAILED"
            die "Could not confidently infer input build from input file (${DETECTED_INPUT_BUILD_REASON}). Verify input build metadata and requested target build. If you are sure of the builds, you can force liftover with --force (at your own risk)."
        else
            if [[ "$DETECTED_INPUT_BUILD" == "$TARGET_BUILD" && "$SOURCE_BUILD" != "$TARGET_BUILD" ]]; then
                if [[ "$FORCE_ALL_LIFTOVER" == true ]]; then
                    warn "Input VCF appears to be ${DETECTED_INPUT_BUILD}, while chain maps ${SOURCE_BUILD} -> ${TARGET_BUILD}. Continuing because --force-all is set."
                else
                    CHAIN_COMPATIBILITY_STATUS="FAILED"
                    die "Input VCF appears to be ${DETECTED_INPUT_BUILD}, while chain maps ${SOURCE_BUILD} -> ${TARGET_BUILD}. Verify the input build and the build you want to lift over to. If you are absolutely sure of your actions, use --force-all to bypass."
                fi
            fi

            if [[ "$DETECTED_INPUT_BUILD" != "$SOURCE_BUILD" ]]; then
                if [[ "$FORCE_ALL_LIFTOVER" == true ]]; then
                    warn "Input VCF appears to be ${DETECTED_INPUT_BUILD}, but chain source build is ${SOURCE_BUILD} (${SOURCE_BUILD} -> ${TARGET_BUILD}). Continuing because --force-all is set."
                else
                    CHAIN_COMPATIBILITY_STATUS="FAILED"
                    die "Input VCF appears to be ${DETECTED_INPUT_BUILD}, but chain source build is ${SOURCE_BUILD} (${SOURCE_BUILD} -> ${TARGET_BUILD}). Verify input build and desired target build. If you are absolutely sure of your actions, use --force-all to bypass."
                fi
            fi

            if [[ "$INPUT_PRIMARY_CHR_COUNT" == "1" && "$INPUT_PRIMARY_CHR_LIST" == "16" ]]; then
                warn "Only chr16 is present. The ~16kb size difference between GRCh37 and GRCh38 for this chromosome is too small for reliable build inference."
                warn "Assuming build ${DETECTED_INPUT_BUILD} entirely based on chain file ${CHAIN_FILE}."
                if [[ "$FORCE_LIFTOVER" != true ]]; then
                    read -p "Do you want to continue anyways? (y/n) " -r
                    if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
                        die "Liftover aborted by user."
                    fi
                else
                    warn "Continuing automatically because --force is set."
                fi
            fi

            CHAIN_COMPATIBILITY_STATUS="PASS (${DETECTED_INPUT_BUILD} matches chain source ${SOURCE_BUILD})"
            ok "Assembly check passed: ${CHAIN_COMPATIBILITY_STATUS}"
        fi
    else
        CHAIN_COMPATIBILITY_STATUS="SKIPPED (chain not recognized as GRCh37/GRCh38)"
        warn "Assembly check skipped: chain direction could not be inferred as GRCh37/GRCh38."
        if [[ "$DETECTED_INPUT_BUILD" == "unknown" ]]; then
            warn "Input build could not be inferred: ${DETECTED_INPUT_BUILD_REASON}"
        else
            info "Detected input build: ${DETECTED_INPUT_BUILD} (${DETECTED_INPUT_BUILD_REASON})"
        fi
    fi

    # -----------------------------------------------------------------------------
    # STEP 2: Query variants from input VCF and build UCSC BED
    # Extracts all coordinate info (CHR, POS, ID, REF, ALT) to generate the 0-based BED.
    # -----------------------------------------------------------------------------
    CURRENT_STEP="2/6 query VCF and build BED"
    step "2/6 Query VCF variants and build UCSC BED"

    VCF_QUERY="${WORK_DIR}/vcf-query.tsv"
    UCSC_BED="${WORK_DIR}/source.bed"
    STEP2_LOG="${WORK_DIR}/vcf-step2-query.log"

    LAST_LOG_FILE="$STEP2_LOG"
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$INPUT_VCF" \
        > "$VCF_QUERY" 2> "$STEP2_LOG"

    TOTAL_VARIANTS_RAW="$(wc -l < "$VCF_QUERY" | tr -d ' ')"
    (( TOTAL_VARIANTS_RAW > 0 )) || die "No variant records found in VCF input."

    VCF_HAS_CHR_PREFIX="$(awk 'NR==1 {if ($1 ~ /^chr/) print 1; else print 0; exit}' "$VCF_QUERY")"

    awk -v OFS='\t' -v keep_all="$KEEP_ALL_CHR" '
        function canon_chr(c, x) {
            x = c
            sub(/^chr/, "", x)
            sub(/^CHR/, "", x)
            ux = toupper(x)
            if (ux == "23") x = "X"
            else if (ux == "24") x = "Y"
            else if (ux == "26" || ux == "M" || ux == "MT") x = "MT"
            else if (ux == "25" || ux == "XY") x = "XY"
            else if (ux == "X" || ux == "Y" || ux ~ /^([1-9]|1[0-9]|2[0-2])$/) x = ux
            return x
        }
        function is_primary(c) {
            uc = toupper(c)
            return (uc ~ /^([1-9]|1[0-9]|2[0-2])$/ || uc == "X" || uc == "Y")
        }
        function norm_chr(c, x) {
            x = canon_chr(c)
            if (x == "MT") return "chrM"
            return "chr" x
        }
        $2 ~ /^[0-9]+$/ {
            cc = canon_chr($1)
            if (keep_all != "true" && !is_primary(cc)) next
            chr = norm_chr($1)
            start = ($2 > 0 ? $2 - 1 : 0)
            print chr, start, $2, NR
        }
    ' "$VCF_QUERY" > "$UCSC_BED"

    TOTAL_VARIANTS="$(wc -l < "$UCSC_BED" | tr -d ' ')"
    (( TOTAL_VARIANTS > 0 )) || die "No variants left after chromosome filtering. Use --keep-all-chr to keep XY/MT/non-primary contigs."
    SOURCE_FILTER_DROPPED="$(( TOTAL_VARIANTS_RAW - TOTAL_VARIANTS ))"
    if [[ "$KEEP_ALL_CHR" == false && "$SOURCE_FILTER_DROPPED" -gt 0 ]]; then
        warn "Dropped ${SOURCE_FILTER_DROPPED} non-primary/XY/MT input variants before liftover."
    fi

    ok "Prepared ${TOTAL_VARIANTS} VCF coordinates"

    # -----------------------------------------------------------------------------
    # STEP 3: liftOver
    # Invokes the UCSC liftOver binary on the generated BED file.
    # Unmapped variants are output to an unmapped file.
    # -----------------------------------------------------------------------------
    CURRENT_STEP="3/6 liftover"
    step "3/6 Run UCSC liftOver"

    LIFTED_BED="${WORK_DIR}/lifted.bed"
    UNMAPPED_RAW="${WORK_DIR}/unmapped.raw"
    STEP3_LOG="${WORK_DIR}/vcf-step3-liftover.log"
    LAST_LOG_FILE="$STEP3_LOG"

    "$LIFTOVER_BIN" "$UCSC_BED" "$CHAIN_FILE" "$LIFTED_BED" "$UNMAPPED_RAW" \
        > "$STEP3_LOG" 2>&1

    if [[ "$KEEP_ALL_CHR" == false ]]; then
        FILTERED_LIFTED_BED="${WORK_DIR}/lifted.primary-target.bed"
        LIFTED_COUNT_RAW="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
        filter_primary_target_contigs "$LIFTED_BED" "$FILTERED_LIFTED_BED"
        mv "$FILTERED_LIFTED_BED" "$LIFTED_BED"
        LIFTED_COUNT_FILTERED="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
        TARGET_FILTER_DROPPED="$(( LIFTED_COUNT_RAW - LIFTED_COUNT_FILTERED ))"
        if (( TARGET_FILTER_DROPPED > 0 )); then
            warn "Dropped ${TARGET_FILTER_DROPPED} lifted variants on non-primary target contigs before rewrite."
        fi
    fi

    LIFTED_COUNT="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
    UNMAPPED_COUNT="$(awk 'NF && $1 !~ /^#/' "$UNMAPPED_RAW" | wc -l | tr -d ' ')"
    UNMAPPED_PCT="$(awk -v u="$UNMAPPED_COUNT" -v t="$TOTAL_VARIANTS" 'BEGIN {if (t > 0) printf "%.3f", (u/t)*100; else print "0.000"}')"

    (( LIFTED_COUNT > 0 )) || die "liftOver mapped 0 variants after target chromosome filtering. Check chain/build compatibility or use --keep-all-chr."
    ok "Lifted ${LIFTED_COUNT}; unmapped ${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"

    UNMAPPED_VCF_BED="${OUTPUT_PREFIX}.unmapped.bed"
    if (( UNMAPPED_COUNT > 0 )); then
        awk 'NF && $1 !~ /^#/' "$UNMAPPED_RAW" > "$UNMAPPED_VCF_BED"
        warn "Unmapped variants written: ${UNMAPPED_VCF_BED}"
    fi

    # -----------------------------------------------------------------------------
    # STEP 4: Build lifted coordinate map
    # Creates a mapping file linking original BED row indices to their newly lifted coordinates.
    # -----------------------------------------------------------------------------
    CURRENT_STEP="4/6 build coordinate map"
    step "4/6 Build coordinate mapping for VCF rewrite"

    COORD_MAP="${WORK_DIR}/vcf-lift-map.tsv"
    awk -v OFS='\t' -v keep_chr="$VCF_HAS_CHR_PREFIX" '
        function to_vcf_chr(c) {
            if (keep_chr == 1) {
                if (c ~ /^chr/) return c
                if (c == "M" || c == "MT") return "chrM"
                return "chr" c
            }
            sub(/^chr/, "", c)
            if (c == "M") c = "MT"
            return c
        }
        { print $4, to_vcf_chr($1), $3 }
    ' "$LIFTED_BED" > "$COORD_MAP"

    # -----------------------------------------------------------------------------
    # STEP 5: Rewrite VCF body with lifted coordinates
    # Uses AWK to iterate over the VCF body and apply the updated positions from the coord map.
    # Also rewrites variant IDs if they match the old position format.
    # -----------------------------------------------------------------------------
    CURRENT_STEP="5/6 rewrite VCF"
    step "5/6 Rewrite VCF with lifted coordinates"

    VCF_HEADER="${WORK_DIR}/vcf-header.txt"
    VCF_BODY="${WORK_DIR}/vcf-body.txt"
    VCF_BODY_LIFTED="${WORK_DIR}/vcf-body-lifted.txt"
    ID_REWRITE_STATS="${WORK_DIR}/vcf-id-rewrite-stats.tsv"
    VCF_UNSORTED="${WORK_DIR}/vcf-unsorted.vcf"
    STEP5_LOG="${WORK_DIR}/vcf-step5-rewrite.log"

    LAST_LOG_FILE="$STEP5_LOG"
    bcftools view -h "$INPUT_VCF" > "$VCF_HEADER" 2> "$STEP5_LOG"
    bcftools view -H "$INPUT_VCF" > "$VCF_BODY" 2>> "$STEP5_LOG"

    awk -v OFS='\t' -v stats_file="$ID_REWRITE_STATS" '
        function canon_chr(c, x) {
            x = toupper(c)
            sub(/^CHR/, "", x)
            if (x == "M") x = "MT"
            return x
        }
        function should_rewrite_id(id, chr, pos, ref, alt, n, c, p, r, a, arr) {
            if (id == "." || id == "") return 0

            n = split(id, arr, ":")
            if (n != 4) n = split(id, arr, "_")
            if (n != 4) return 0

            c = arr[1]
            p = arr[2]
            r = arr[3]
            a = arr[4]

            if (p !~ /^[0-9]+$/) return 0
            if (canon_chr(c) != canon_chr(chr)) return 0
            if ((p + 0) != (pos + 0)) return 0
            if (toupper(r) != toupper(ref)) return 0
            if (toupper(a) != toupper(alt)) return 0
            return 1
        }
        NR==FNR {chr[$1]=$2; pos[$1]=$3; next}
        {
            idx = FNR
            if (idx in chr) {
                orig_chr = $1
                orig_pos = $2
                orig_id = $3
                if (should_rewrite_id(orig_id, orig_chr, orig_pos, $4, $5)) {
                    $3 = chr[idx] ":" pos[idx] ":" $4 ":" $5
                    rewritten++
                } else {
                    preserved++
                }
                $1 = chr[idx]
                $2 = pos[idx]
                print
            }
        }
        END {
            print (rewritten + 0) "\t" (preserved + 0) > stats_file
        }
    ' "$COORD_MAP" "$VCF_BODY" > "$VCF_BODY_LIFTED"

    if [[ -s "$ID_REWRITE_STATS" ]]; then
        ID_REWRITTEN_COUNT="$(awk 'NR==1 {print $1 + 0}' "$ID_REWRITE_STATS")"
        ID_PRESERVED_COUNT="$(awk 'NR==1 {print $2 + 0}' "$ID_REWRITE_STATS")"
    fi

    awk -v chain="$CHAIN_FILE" '
        /^#CHROM/ {
            print "##liftoverTool=UCSC_liftOver"
            print "##liftoverChain=" chain
            print "##liftoverNote=Coordinate update with positional-ID synchronization; verify REF/ALT consistency against target reference."
        }
        { print }
    ' "$VCF_HEADER" > "$VCF_UNSORTED"
    cat "$VCF_BODY_LIFTED" >> "$VCF_UNSORTED"

    # -----------------------------------------------------------------------------
    # STEP 6: Sort and write final VCF
    # Leverages bcftools to sort the final VCF output and generate index files (.csi/.tbi).
    # -----------------------------------------------------------------------------
    CURRENT_STEP="6/6 sort and write output"
    step "6/6 Sort and write lifted VCF"

    STEP6_LOG="${WORK_DIR}/vcf-step6-sort-write.log"
    LAST_LOG_FILE="$STEP6_LOG"

    case "$OUTPUT_VCF" in
        *.vcf.gz)
            bcftools sort -Oz -o "$OUTPUT_VCF" "$VCF_UNSORTED" > "$STEP6_LOG" 2>&1
            bcftools index -f "$OUTPUT_VCF" >> "$STEP6_LOG" 2>&1
            ;;
        *.vcf)
            bcftools sort -Ov -o "$OUTPUT_VCF" "$VCF_UNSORTED" > "$STEP6_LOG" 2>&1
            ;;
        *)
            die "Unsupported VCF output extension for -O. Use .vcf or .vcf.gz"
            ;;
    esac

    FINAL_VARIANTS="$(bcftools view -H "$OUTPUT_VCF" | wc -l | tr -d ' ')"

    line
    printf "%bRun Summary%b\n" "$BOLD" "$RESET"
    line
    print_kv "Mode" "VCF"
    print_kv "Direction" "${SOURCE_BUILD} -> ${TARGET_BUILD}"
    print_kv "Input build" "$DETECTED_INPUT_BUILD"
    print_kv "Primary chr set" "${INPUT_PRIMARY_CHR_LIST} (${INPUT_PRIMARY_CHR_COUNT})"
    print_kv "Chain source" "$SOURCE_BUILD"
    print_kv "Chain target" "$TARGET_BUILD"
    print_kv "Assembly check" "$CHAIN_COMPATIBILITY_STATUS"
    print_kv "Force mode" "$([[ "$FORCE_LIFTOVER" == true ]] && echo yes || echo no)"
    print_kv "Force all mode" "$([[ "$FORCE_ALL_LIFTOVER" == true ]] && echo yes || echo no)"
    print_kv "Input variants" "$TOTAL_VARIANTS"
    print_kv "Lifted variants" "$LIFTED_COUNT"
    print_kv "Unmapped" "${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"
    print_kv "Final variants" "$FINAL_VARIANTS"
    print_kv "IDs rewritten" "$ID_REWRITTEN_COUNT"
    print_kv "IDs preserved" "$ID_PRESERVED_COUNT"
    print_kv "Output VCF" "$OUTPUT_VCF"
    if [[ -f "${OUTPUT_VCF}.csi" ]]; then
        print_kv "Index" "${OUTPUT_VCF}.csi"
    elif [[ -f "${OUTPUT_VCF}.tbi" ]]; then
        print_kv "Index" "${OUTPUT_VCF}.tbi"
    fi
    if [[ -f "$UNMAPPED_VCF_BED" ]]; then
        print_kv "Unmapped BED" "$UNMAPPED_VCF_BED"
    fi
    print_kv "Run log" "$RUN_SUMMARY_FILE"
    line

    write_success_log_header
    log_kv "Mode" "VCF"
    log_kv "Direction" "${SOURCE_BUILD} -> ${TARGET_BUILD}"
    log_kv "Input build" "$DETECTED_INPUT_BUILD"
    log_kv "Build evidence" "$DETECTED_INPUT_BUILD_REASON"
    log_kv "Primary chr set" "${INPUT_PRIMARY_CHR_LIST} (${INPUT_PRIMARY_CHR_COUNT})"
    log_kv "Chain source" "$SOURCE_BUILD"
    log_kv "Chain target" "$TARGET_BUILD"
    log_kv "Assembly check" "$CHAIN_COMPATIBILITY_STATUS"
    log_kv "Force mode" "$([[ "$FORCE_LIFTOVER" == true ]] && echo yes || echo no)"
    log_kv "Force all mode" "$([[ "$FORCE_ALL_LIFTOVER" == true ]] && echo yes || echo no)"
    log_kv "Input VCF" "$INPUT_VCF"
    log_kv "Output VCF" "$OUTPUT_VCF"
    log_kv "Chain file" "$CHAIN_FILE"
    log_kv "Input variants" "$TOTAL_VARIANTS"
    log_kv "Lifted variants" "$LIFTED_COUNT"
    log_kv "Unmapped" "${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"
    log_kv "Final variants" "$FINAL_VARIANTS"
    log_kv "IDs rewritten" "$ID_REWRITTEN_COUNT"
    log_kv "IDs preserved" "$ID_PRESERVED_COUNT"
    if [[ -f "${OUTPUT_VCF}.csi" ]]; then
        log_kv "Index" "${OUTPUT_VCF}.csi"
    elif [[ -f "${OUTPUT_VCF}.tbi" ]]; then
        log_kv "Index" "${OUTPUT_VCF}.tbi"
    fi
    if [[ -f "$UNMAPPED_VCF_BED" ]]; then
        log_kv "Unmapped BED" "$UNMAPPED_VCF_BED"
    fi
    log_kv "Temp dir kept" "$([[ "$KEEP_INTERMEDIATE" == true ]] && echo yes || echo no)"
    log_line

    ok "VCF liftover completed successfully"
    exit 0
fi

# =============================================================================
# PLINK MODE PIPELINE
# 1) Assembly Check -> 2) Build UCSC BED (and optionally deduplicate) -> 
# 3) liftOver -> 4) Build update files -> 5) Apply updates via plink2
# =============================================================================

# -----------------------------------------------------------------------------
# STEP 1: Detect assembly build and validate chain direction
# Validate that the input variant coordinates match the chain source build.
# -----------------------------------------------------------------------------
CURRENT_STEP="1/5 assembly check"
step "1/5 Detect input assembly and validate chain direction"

STEP1_LOG="${WORK_DIR}/plink-step1-assembly-check.log"
LAST_LOG_FILE="$STEP1_LOG"
: > "$STEP1_LOG"

IFS=$'\t' read -r DETECTED_INPUT_BUILD DETECTED_INPUT_BUILD_REASON INPUT_PRIMARY_CHR_COUNT INPUT_PRIMARY_CHR_LIST \
    < <(detect_plink_build "$INPUT_FORMAT" "$VARIANT_FILE" "$SOURCE_BUILD")

[[ -n "$DETECTED_INPUT_BUILD" ]] || DETECTED_INPUT_BUILD="unknown"
[[ -n "$DETECTED_INPUT_BUILD_REASON" ]] || DETECTED_INPUT_BUILD_REASON="Unable to parse build markers from variant coordinates"
[[ "$INPUT_PRIMARY_CHR_COUNT" =~ ^[0-9]+$ ]] || INPUT_PRIMARY_CHR_COUNT=0
[[ -n "$INPUT_PRIMARY_CHR_LIST" ]] || INPUT_PRIMARY_CHR_LIST="none"

{
    echo "Detected build: ${DETECTED_INPUT_BUILD}"
    echo "Build evidence: ${DETECTED_INPUT_BUILD_REASON}"
    echo "Observed primary chromosome count: ${INPUT_PRIMARY_CHR_COUNT}"
    echo "Observed primary chromosomes: ${INPUT_PRIMARY_CHR_LIST}"
    echo "Chain source build: ${SOURCE_BUILD}"
    echo "Chain target build: ${TARGET_BUILD}"
} >> "$STEP1_LOG"

if [[ "$CHAIN_BUILD_GUARD_ENABLED" == true ]]; then
    if [[ "$DETECTED_INPUT_BUILD" == "unknown" ]]; then
        CHAIN_COMPATIBILITY_STATUS="FAILED"
        die "Could not confidently infer input build from input file (${DETECTED_INPUT_BUILD_REASON}). Verify input build metadata and requested target build. If you are sure of the builds, you can force liftover with --force (at your own risk)."
    else
        if [[ "$DETECTED_INPUT_BUILD" == "$TARGET_BUILD" && "$SOURCE_BUILD" != "$TARGET_BUILD" ]]; then
            if [[ "$FORCE_ALL_LIFTOVER" == true ]]; then
                warn "Input PLINK appears to be ${DETECTED_INPUT_BUILD}, while chain maps ${SOURCE_BUILD} -> ${TARGET_BUILD}. Continuing because --force-all is set."
            else
                CHAIN_COMPATIBILITY_STATUS="FAILED"
                die "Input PLINK appears to be ${DETECTED_INPUT_BUILD}, while chain maps ${SOURCE_BUILD} -> ${TARGET_BUILD}. Verify the input build and the build you want to lift over to. If you are absolutely sure of your actions, use --force-all to bypass."
            fi
        fi

        if [[ "$DETECTED_INPUT_BUILD" != "$SOURCE_BUILD" ]]; then
            if [[ "$FORCE_ALL_LIFTOVER" == true ]]; then
                warn "Input PLINK appears to be ${DETECTED_INPUT_BUILD}, but chain source build is ${SOURCE_BUILD} (${SOURCE_BUILD} -> ${TARGET_BUILD}). Continuing because --force-all is set."
            else
                CHAIN_COMPATIBILITY_STATUS="FAILED"
                die "Input PLINK appears to be ${DETECTED_INPUT_BUILD}, but chain source build is ${SOURCE_BUILD} (${SOURCE_BUILD} -> ${TARGET_BUILD}). Verify input build and desired target build. If you are absolutely sure of your actions, use --force-all to bypass."
            fi
        fi

        if [[ "$INPUT_PRIMARY_CHR_COUNT" == "1" && "$INPUT_PRIMARY_CHR_LIST" == "16" ]]; then
            warn "Only chr16 is present. The ~16kb size difference between GRCh37 and GRCh38 for this chromosome is too small for reliable build inference."
            warn "Assuming build ${DETECTED_INPUT_BUILD} entirely based on chain file ${CHAIN_FILE}."
            if [[ "$FORCE_LIFTOVER" != true ]]; then
                read -p "Do you want to continue anyways? (y/n) " -r
                if [[ ! "$REPLY" =~ ^[Yy]$ ]]; then
                    die "Liftover aborted by user."
                fi
            else
                warn "Continuing automatically because --force is set."
            fi
        fi

        CHAIN_COMPATIBILITY_STATUS="PASS (${DETECTED_INPUT_BUILD} matches chain source ${SOURCE_BUILD})"
        ok "Assembly check passed: ${CHAIN_COMPATIBILITY_STATUS}"
    fi
else
    CHAIN_COMPATIBILITY_STATUS="SKIPPED (chain not recognized as GRCh37/GRCh38)"
    warn "Assembly check skipped: chain direction could not be inferred as GRCh37/GRCh38."
    if [[ "$DETECTED_INPUT_BUILD" == "unknown" ]]; then
        warn "Input build could not be inferred: ${DETECTED_INPUT_BUILD_REASON}"
    else
        info "Detected input build: ${DETECTED_INPUT_BUILD} (${DETECTED_INPUT_BUILD_REASON})"
    fi
fi

# -----------------------------------------------------------------------------
# STEP 2: variant table -> UCSC BED
# Converts the variant lists into UCSC BED format.
# Handles deduplication (via plink2 --rm-dup) if identical variant IDs exist,
# to avoid mapping confusion.
# -----------------------------------------------------------------------------
CURRENT_STEP="2/5 extract variant coordinates"
step "2/5 Extract coordinates in UCSC BED format"

UCSC_BED="${WORK_DIR}/source.bed"
build_ucsc_bed "$INPUT_FORMAT" "$VARIANT_FILE" "$UCSC_BED"

TOTAL_VARIANTS="$(wc -l < "$UCSC_BED" | tr -d ' ')"
(( TOTAL_VARIANTS > 0 )) || die "No valid variants left after chromosome filtering. Use --keep-all-chr to keep XY/MT/non-primary contigs."

DUP_ID_COUNT="$(awk '{print $4}' "$UCSC_BED" | LC_ALL=C sort | uniq -d | wc -l | tr -d ' ')"
if (( DUP_ID_COUNT > 0 )); then
    ORIGINAL_TOTAL_VARIANTS="$TOTAL_VARIANTS"
    warn "Detected ${DUP_ID_COUNT} duplicated input variant IDs."
    warn "Duplicate IDs will be removed with plink2 --rm-dup exclude-mismatch."
    CURRENT_STEP="2/5 deduplicate duplicate IDs"
    step "2/5 Remove duplicate IDs with plink2 --rm-dup exclude-mismatch"

    DEDUP_PREFIX="${WORK_DIR}/tmp-rmdup"
    STEP1_RMDUP_LOG="${WORK_DIR}/plink2-step1-rmdup.log"
    LAST_LOG_FILE="$STEP1_RMDUP_LOG"
    PLINK2_ACTIVE_OUT_PREFIX="$DEDUP_PREFIX"
    PLINK2_ACTIVE_STEP_LOG="$STEP1_RMDUP_LOG"

    plink2 "${PLINK_INPUT_ARGS[@]}" \
        --rm-dup exclude-mismatch \
        --make-pgen \
        --out "$DEDUP_PREFIX" \
        > "$STEP1_RMDUP_LOG" 2>&1

    PLINK2_ACTIVE_OUT_PREFIX=""
    PLINK2_ACTIVE_STEP_LOG=""

    if [[ -f "${DEDUP_PREFIX}.log" ]]; then
        cp "${DEDUP_PREFIX}.log" "${LOG_PREFIX}.rmdup.log" 2>/dev/null || true
    else
        cp "$STEP1_RMDUP_LOG" "${LOG_PREFIX}.rmdup.log" 2>/dev/null || true
    fi

    PLINK_INPUT_ARGS=(--pfile "$DEDUP_PREFIX")
    INPUT_FORMAT="plink2"
    VARIANT_FILE="${DEDUP_PREFIX}.pvar"
    build_ucsc_bed "$INPUT_FORMAT" "$VARIANT_FILE" "$UCSC_BED"

    TOTAL_VARIANTS="$(wc -l < "$UCSC_BED" | tr -d ' ')"
    (( TOTAL_VARIANTS > 0 )) || die "No variants left after duplicate-ID removal."

    DUP_ID_COUNT="$(awk '{print $4}' "$UCSC_BED" | LC_ALL=C sort | uniq -d | wc -l | tr -d ' ')"
    (( DUP_ID_COUNT == 0 )) || die "Duplicate variant IDs remain after --rm-dup exclude-mismatch. See ${LOG_PREFIX}.rmdup.log"

    REMOVED_DUPS="$(( ORIGINAL_TOTAL_VARIANTS - TOTAL_VARIANTS ))"
    ok "Removed ${REMOVED_DUPS} variants during duplicate-ID cleanup. Log: ${LOG_PREFIX}.rmdup.log"
fi

ok "Prepared ${TOTAL_VARIANTS} coordinates"

# -----------------------------------------------------------------------------
# STEP 3: liftOver
# Invokes UCSC liftOver binary on the BED extracted from PLINK files.
# -----------------------------------------------------------------------------
CURRENT_STEP="3/5 liftover"
step "3/5 Run UCSC liftOver"

LIFTED_BED="${WORK_DIR}/lifted.bed"
UNMAPPED_RAW="${WORK_DIR}/unmapped.raw"
STEP2_LOG="${WORK_DIR}/liftover-step2.log"

LAST_LOG_FILE="$STEP2_LOG"
"$LIFTOVER_BIN" "$UCSC_BED" "$CHAIN_FILE" "$LIFTED_BED" "$UNMAPPED_RAW" \
    > "$STEP2_LOG" 2>&1

if [[ "$KEEP_ALL_CHR" == false ]]; then
    FILTERED_LIFTED_BED="${WORK_DIR}/lifted.primary-target.bed"
    LIFTED_COUNT_RAW="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
    filter_primary_target_contigs "$LIFTED_BED" "$FILTERED_LIFTED_BED"
    mv "$FILTERED_LIFTED_BED" "$LIFTED_BED"
    LIFTED_COUNT_FILTERED="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
    TARGET_FILTER_DROPPED="$(( LIFTED_COUNT_RAW - LIFTED_COUNT_FILTERED ))"
    if (( TARGET_FILTER_DROPPED > 0 )); then
        warn "Dropped ${TARGET_FILTER_DROPPED} lifted variants on non-primary target contigs before PLINK updates."
    fi
fi

LIFTED_COUNT="$(wc -l < "$LIFTED_BED" | tr -d ' ')"
UNMAPPED_COUNT="$(awk 'NF && $1 !~ /^#/' "$UNMAPPED_RAW" | wc -l | tr -d ' ')"
UNMAPPED_PCT="$(awk -v u="$UNMAPPED_COUNT" -v t="$TOTAL_VARIANTS" 'BEGIN {if (t > 0) printf "%.3f", (u/t)*100; else print "0.000"}')"

(( LIFTED_COUNT > 0 )) || die "liftOver mapped 0 variants after target chromosome filtering. Check chain/build compatibility or use --keep-all-chr."
ok "Lifted ${LIFTED_COUNT}; unmapped ${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"

UNMAPPED_BED_OUT="${OUTPUT_PREFIX}.unmapped.bed"
if (( UNMAPPED_COUNT > 0 )); then
    awk 'NF && $1 !~ /^#/' "$UNMAPPED_RAW" > "$UNMAPPED_BED_OUT"
    warn "Unmapped variants written: ${UNMAPPED_BED_OUT}"
fi

# -----------------------------------------------------------------------------
# STEP 4: Build PLINK update/extract files
# Produces formatted tables for --update-map and --update-chr, as well as an
# ID extract list to keep only variants that successfully mapped.
# -----------------------------------------------------------------------------
CURRENT_STEP="4/5 build update files"
step "4/5 Build chromosome/position update files"

POS_UPDATE="${WORK_DIR}/update-pos.tsv"   # variant_id  new_bp
CHR_UPDATE="${WORK_DIR}/update-chr.tsv"   # variant_id  new_chr
ID_EXTRACT="${WORK_DIR}/lifted-ids.txt"   # variant_id

awk -v OFS='\t' '{print $4, $3}' "$LIFTED_BED" > "$POS_UPDATE"
awk -v OFS='\t' '{
    chr = $1
    sub(/^chr/, "", chr)
    if (chr == "M") chr = "MT"
    print $4, chr
}' "$LIFTED_BED" > "$CHR_UPDATE"
awk '{print $4}' "$LIFTED_BED" > "$ID_EXTRACT"

ok "Update tables generated"

# -----------------------------------------------------------------------------
# STEP 5: Apply updates in PLINK and produce final output
# Uses plink2 to apply the coordinate/chromosome updates, rename variant IDs,
# and write the final sorted dataset in the requested format (pgen/bed).
# -----------------------------------------------------------------------------
CURRENT_STEP="5/5 apply updates"
step "5/5 Apply updates, sort in regular .pgen, then write final dataset"

TEMP_PREFIX="${WORK_DIR}/tmp-updated"
SORTED_PREFIX="${WORK_DIR}/tmp-sorted"
STEP4A_LOG="${WORK_DIR}/plink2-step4a.log"
STEP4B_LOG="${WORK_DIR}/plink2-step4b.log"
STEP4C_LOG="${WORK_DIR}/plink2-step4c.log"

LAST_LOG_FILE="$STEP4A_LOG"
PLINK2_ACTIVE_OUT_PREFIX="$TEMP_PREFIX"
PLINK2_ACTIVE_STEP_LOG="$STEP4A_LOG"
plink2 "${PLINK_INPUT_ARGS[@]}" \
    --update-map "$POS_UPDATE" \
    --update-chr "$CHR_UPDATE" \
    --extract "$ID_EXTRACT" \
    --sort-vars \
    --make-pgen \
    --out "$TEMP_PREFIX" \
    > "$STEP4A_LOG" 2>&1
PLINK2_ACTIVE_OUT_PREFIX=""
PLINK2_ACTIVE_STEP_LOG=""

LAST_LOG_FILE="$STEP4B_LOG"
PLINK2_ACTIVE_OUT_PREFIX="$SORTED_PREFIX"
PLINK2_ACTIVE_STEP_LOG="$STEP4B_LOG"
plink2 --pfile "$TEMP_PREFIX" \
    --set-all-var-ids "$VARID_TEMPLATE" \
    --sort-vars \
    --make-pgen \
    --out "$SORTED_PREFIX" \
    > "$STEP4B_LOG" 2>&1
PLINK2_ACTIVE_OUT_PREFIX=""
PLINK2_ACTIVE_STEP_LOG=""

LAST_LOG_FILE="$STEP4C_LOG"
PLINK2_ACTIVE_OUT_PREFIX="$OUTPUT_PREFIX"
PLINK2_ACTIVE_STEP_LOG="$STEP4C_LOG"
plink2 --pfile "$SORTED_PREFIX" \
    "${PLINK_OUTPUT_ARGS[@]}" \
    --out "$OUTPUT_PREFIX" \
    > "$STEP4C_LOG" 2>&1
PLINK2_ACTIVE_OUT_PREFIX=""
PLINK2_ACTIVE_STEP_LOG=""

# -----------------------------------------------------------------------------
# Final counts + summary
# -----------------------------------------------------------------------------
if [[ "$OUTPUT_FORMAT" == "plink2" ]]; then
    FINAL_VARIANTS="$(grep -vc '^#' "${OUTPUT_PREFIX}.pvar")"
else
    FINAL_VARIANTS="$(wc -l < "${OUTPUT_PREFIX}.bim" | tr -d ' ')"
fi

line
printf "%bRun Summary%b\n" "$BOLD" "$RESET"
line
print_kv "Direction" "${SOURCE_BUILD} -> ${TARGET_BUILD}"
print_kv "Input build" "$DETECTED_INPUT_BUILD"
print_kv "Primary chr set" "${INPUT_PRIMARY_CHR_LIST} (${INPUT_PRIMARY_CHR_COUNT})"
print_kv "Chain source" "$SOURCE_BUILD"
print_kv "Chain target" "$TARGET_BUILD"
print_kv "Assembly check" "$CHAIN_COMPATIBILITY_STATUS"
print_kv "Force mode" "$([[ "$FORCE_LIFTOVER" == true ]] && echo yes || echo no)"
print_kv "Force all mode" "$([[ "$FORCE_ALL_LIFTOVER" == true ]] && echo yes || echo no)"
print_kv "Input variants" "$TOTAL_VARIANTS"
print_kv "Lifted variants" "$LIFTED_COUNT"
print_kv "Unmapped" "${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"
print_kv "Final variants" "$FINAL_VARIANTS"
print_kv "Variant IDs" "$VARID_TEMPLATE"

if [[ "$OUTPUT_FORMAT" == "plink2" ]]; then
    print_kv "Output" "${OUTPUT_PREFIX}.pgen"
    print_kv "Output" "${OUTPUT_PREFIX}.pvar"
    print_kv "Output" "${OUTPUT_PREFIX}.psam"
else
    print_kv "Output" "${OUTPUT_PREFIX}.bed"
    print_kv "Output" "${OUTPUT_PREFIX}.bim"
    print_kv "Output" "${OUTPUT_PREFIX}.fam"
fi

if [[ -f "$UNMAPPED_BED_OUT" ]]; then
    print_kv "Unmapped BED" "$UNMAPPED_BED_OUT"
fi
print_kv "Run log" "$RUN_SUMMARY_FILE"

line

write_success_log_header
log_kv "Mode" "PLINK"
log_kv "Direction" "${SOURCE_BUILD} -> ${TARGET_BUILD}"
log_kv "Input build" "$DETECTED_INPUT_BUILD"
log_kv "Build evidence" "$DETECTED_INPUT_BUILD_REASON"
log_kv "Primary chr set" "${INPUT_PRIMARY_CHR_LIST} (${INPUT_PRIMARY_CHR_COUNT})"
log_kv "Chain source" "$SOURCE_BUILD"
log_kv "Chain target" "$TARGET_BUILD"
log_kv "Assembly check" "$CHAIN_COMPATIBILITY_STATUS"
log_kv "Force mode" "$([[ "$FORCE_LIFTOVER" == true ]] && echo yes || echo no)"
log_kv "Force all mode" "$([[ "$FORCE_ALL_LIFTOVER" == true ]] && echo yes || echo no)"
log_kv "Input prefix" "$INPUT_PREFIX"
log_kv "Input format" "$INPUT_FORMAT"
log_kv "Output prefix" "$OUTPUT_PREFIX"
log_kv "Output format" "$OUTPUT_FORMAT"
log_kv "Chain file" "$CHAIN_FILE"
log_kv "Var-ID template" "$VARID_TEMPLATE"
log_kv "Input variants" "$TOTAL_VARIANTS"
log_kv "Lifted variants" "$LIFTED_COUNT"
log_kv "Unmapped" "${UNMAPPED_COUNT} (${UNMAPPED_PCT}%)"
log_kv "Final variants" "$FINAL_VARIANTS"
if [[ -f "$UNMAPPED_BED_OUT" ]]; then
    log_kv "Unmapped BED" "$UNMAPPED_BED_OUT"
fi
if [[ -f "${LOG_PREFIX}.rmdup.log" ]]; then
    log_kv "RM-dup log" "${LOG_PREFIX}.rmdup.log"
fi
log_kv "Temp dir kept" "$([[ "$KEEP_INTERMEDIATE" == true ]] && echo yes || echo no)"
log_line

ok "Liftover completed successfully"

# Remove the final plink.log on success
rm -f "${OUTPUT_PREFIX}.log"
