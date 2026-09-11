#!/usr/bin/env bash
set -euo pipefail

command_name="$(basename "$0")"
case "$command_name" in
    ldsc.py|munge_sumstats.py)
        ;;
    *)
        printf 'ERROR: unsupported LDSC wrapper command: %s\n' "$command_name" >&2
        exit 2
        ;;
esac

exec micromamba run -n ldsc "$command_name" "$@"
