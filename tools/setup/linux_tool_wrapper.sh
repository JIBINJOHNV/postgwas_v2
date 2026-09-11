#!/usr/bin/env bash
set -euo pipefail

install_prefix="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
command_name="$(basename "$0")"
binary_path="$install_prefix/share/postgwas/software/linux-tools/bin/$command_name"

if [[ ! -x "$binary_path" ]]; then
    printf 'ERROR: installed Linux tool binary is unavailable: %s\n' \
        "$binary_path" >&2
    exit 1
fi

export LD_LIBRARY_PATH="$install_prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$binary_path" "$@"
