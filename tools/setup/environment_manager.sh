#!/usr/bin/env bash

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
