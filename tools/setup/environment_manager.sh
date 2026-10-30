#!/usr/bin/env bash

# Resolve the inline CEP-24 OS selectors used by the canonical manifest before
# passing it to clients that otherwise treat selector keys as pip-like installers.
# Do not copy package names into shell code or silently accept unknown selectors.
postgwas_render_environment() {
    local source_file="$1"
    local host_os="$2"
    local target
    case "$host_os" in
        Linux) target=linux ;;
        Darwin) target=osx ;;
        *) printf 'ERROR: unsupported environment target: %s\n' "$host_os" >&2; return 2 ;;
    esac
    awk -v target="$target" '
        /^[[:space:]]*-[[:space:]]*sel\(/ {
            if ($0 !~ /^[[:space:]]*-[[:space:]]*sel\((linux|osx|unix|win)\):[[:space:]]*[^[:space:]#]/) {
                print "ERROR: unsupported or empty environment selector at line " NR > "/dev/stderr"
                exit 2
            }
            selector = $0
            sub(/^[^(]*\(/, "", selector)
            sub(/\).*/, "", selector)
            if (selector == target || selector == "unix") {
                sub(/sel\([^)]*\):[[:space:]]*/, "")
                print
            }
            next
        }
        { print }
    ' "$source_file"
}

# Shared by the bootstrap and isolated-runtime installers. Keep the native
# command first; only argparse's rejection of --yes permits a compatibility
# retry, before any environment transaction has started.
postgwas_environment() (
    local manager="$1"
    local operation="$2"
    shift 2
    local diagnostics
    local statuses
    diagnostics="$(mktemp "${TMPDIR:-/tmp}/postgwas-env.XXXXXX")" || return
    trap 'rm -f -- "$diagnostics"' EXIT

    # Stream diagnostics while retaining the exact parser error for inspection.
    if "$manager" env "$operation" --yes "$@" 2>&1 | tee "$diagnostics"; then
        return 0
    else
        statuses=("${PIPESTATUS[@]}")
    fi
    if [[ "${statuses[1]}" -ne 0 ]]; then
        return "${statuses[1]}"
    fi
    if [[ "${statuses[0]}" -ne 2 ]] || ! grep -Eq \
        '^[^[:space:]]+: error: unrecognized arguments: --yes[[:space:]]*$' \
        "$diagnostics"; then
        return "${statuses[0]}"
    fi

    printf 'Retrying without --yes: %s env %s rejected that option; using noninteractive environment settings.\n' \
        "$manager" "$operation"
    # Equivalent always_yes settings for Conda/legacy Mamba and current Mamba:
    # https://www.anaconda.com/blog/conda-configuration-engine-power-users
    # https://mamba.readthedocs.io/en/latest/user_guide/configuration.html
    # Scope them to the retry; do not change the user's configuration files.
    CONDA_ALWAYS_YES=true MAMBA_ALWAYS_YES=true \
        "$manager" env "$operation" "$@"
)
