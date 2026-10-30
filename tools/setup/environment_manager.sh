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

# Reuse native/source components only with matching content records. Staging is
# invocation-local; the top-level installer publishes records after verification.
postgwas_install_state() {
    [[ -n "${POSTGWAS_INSTALL_STATE_STAGE:-}" ]] || return 0
    "${CONDA_PREFIX:?}/bin/python" -I \
        "${POSTGWAS_SETUP_DIRECTORY:?}/installation_state.py" "$1" \
        --prefix "$CONDA_PREFIX" --staging "$POSTGWAS_INSTALL_STATE_STAGE" "${@:2}"
}

postgwas_component() {
    local component="$1"
    shift
    local arguments=(--component "$component"
        --input "$POSTGWAS_SETUP_DIRECTORY/environment_manager.sh"
        --input "$POSTGWAS_SETUP_DIRECTORY/installation_state.py")
    while [[ "${1:-}" != -- ]]; do
        [[ $# -gt 0 ]] || { printf 'ERROR: missing component command.\n' >&2; return 2; }
        arguments+=("$1")
        shift
    done
    shift
    if [[ "${POSTGWAS_INSTALL_UPDATE:-false}" == true && -n "${POSTGWAS_INSTALL_STATE_STAGE:-}" ]]; then
        if postgwas_install_state check "${arguments[@]}"; then
            return 0
        fi
    fi
    printf 'Installing component: %s\n' "$component"
    "$@"
    postgwas_install_state record "${arguments[@]}"
}

postgwas_runtime_environment() {
    local manager="$1" prefix="$2" manifest="$3"
    shift 3
    local operation=create
    if [[ -d "$prefix/conda-meta" ]]; then
        operation=update
    elif [[ -e "$prefix" ]]; then
        printf 'ERROR: existing runtime is not a Conda environment: %s\n' "$prefix" >&2
        return 2
    fi
    postgwas_environment "$manager" "$operation" "$@" --prefix "$prefix" --file "$manifest"
}

# Replace a complete installed source/library tree so removed upstream files
# cannot survive an update. Retain the old tree, including any local additions.
postgwas_replace_directory() (
    local source="$1" destination="$2" staging backup=""
    [[ -d "$source" ]] || { printf 'ERROR: source directory is missing: %s\n' "$source" >&2; return 2; }
    if [[ -e "$destination" && ! -d "$destination" ]]; then
        printf 'ERROR: installed path is not a directory: %s\n' "$destination" >&2
        return 2
    fi
    mkdir -p "$(dirname "$destination")"
    staging="$(mktemp -d "${destination}.staging.XXXXXX")"
    trap 'rm -rf -- "$staging"' EXIT
    cp -R "$source/." "$staging"
    if [[ -e "$destination" || -L "$destination" ]]; then
        backup="$(mktemp -d "${destination}.previous.XXXXXX")"
        rmdir "$backup"
        mv "$destination" "$backup"
    fi
    if ! mv "$staging" "$destination"; then
        [[ -z "$backup" ]] || mv "$backup" "$destination"
        return 1
    fi
    [[ -z "$backup" ]] || printf 'Previous software directory retained: %s\n' "$backup"
)
