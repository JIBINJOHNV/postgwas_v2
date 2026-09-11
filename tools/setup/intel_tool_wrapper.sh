#!/usr/bin/env bash
set -euo pipefail

command_name="$(basename "$0")"
case "$command_name" in
    finemap|ldstore)
        ;;
    *)
        printf 'ERROR: unsupported Intel-mac tool wrapper name: %s\n' \
            "$command_name" >&2
        exit 1
        ;;
esac

install_prefix="$(cd "$(dirname "$0")/.." && pwd -P)"
runtime_root="$install_prefix/share/postgwas/environments/intel-tools"
software_root="$install_prefix/share/postgwas/software/intel-tools"
runtime_library="$software_root/gcc/runtime"
tool_binary="$software_root/bin/$command_name"

for required_path in \
    "$runtime_root/lib" \
    "$runtime_root/lib/libzstd.1.dylib" \
    "$runtime_root/lib/libopenblas.0.dylib" \
    "$runtime_root/lib/libomp.dylib" \
    "$runtime_library/libstdc++.6.dylib" \
    "$runtime_library/libgfortran.5.dylib" \
    "$runtime_library/libquadmath.0.dylib" \
    "$runtime_library/libgcc_s.1.dylib" \
    "$runtime_library/libgomp.1.dylib" \
    "$tool_binary"; do
    if [[ ! -e "$required_path" ]]; then
        printf 'ERROR: required Intel-mac runtime path is missing: %s\n' \
            "$required_path" >&2
        exit 1
    fi
done

exec /usr/bin/arch -x86_64 "$tool_binary" "$@"
