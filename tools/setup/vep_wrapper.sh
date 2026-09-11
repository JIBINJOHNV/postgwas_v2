#!/usr/bin/env bash
set -euo pipefail

wrapper_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
install_prefix="$(cd "$wrapper_directory/.." && pwd -P)"
vep_environment="$install_prefix/share/postgwas/environments/vep"
vep_executable="$vep_environment/bin/vep"
perl_executable="$vep_environment/bin/perl"

for executable in "$vep_executable" "$perl_executable"; do
    if [[ ! -x "$executable" ]]; then
        printf 'ERROR: the installed VEP runtime is incomplete: %s\n' \
            "$executable" >&2
        exit 1
    fi
done

# Keep the isolated VEP runtime independent of user and parent Conda Perl state.
unset PERL5LIB PERLLIB PERL5OPT PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT
export PATH="$vep_environment/bin:$PATH"

exec "$perl_executable" "$vep_executable" "$@"
