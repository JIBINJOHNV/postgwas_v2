#!/usr/bin/env bash
set -Eeuo pipefail

# Isolated runtimes must resolve only packages installed for this checkout.
export PYTHONNOUSERSITE=1
unset PYTHONPATH

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd)"
# shellcheck source=tools/setup/software_versions.env
source "$script_directory/software_versions.env"

mamba_executable=""
build_jobs="${POSTGWAS_INSTALL_JOBS:-2}"

usage() {
    cat <<'EOF'
Usage: install_external_tools.sh --mamba PATH [--jobs N]

Install checksum- or commit-verified scientific programs into the active
PostGWAS Conda prefix on Apple-silicon macOS or glibc Linux x86-64. This
script is called by install_postgwas.sh --all-tools; it is not a standalone
environment creator.
EOF
}

while (($#)); do
    case "$1" in
        --mamba)
            mamba_executable="${2:-}"
            shift 2
            ;;
        --jobs)
            build_jobs="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'ERROR: unknown option: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$mamba_executable" || ! -x "$mamba_executable" ]]; then
    printf 'ERROR: --mamba must name an executable Mamba binary.\n' >&2
    exit 2
fi
if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: --jobs must be a positive integer; found %s.\n' \
        "$build_jobs" >&2
    exit 2
fi

install_prefix="${CONDA_PREFIX:?Run this script inside the target Conda environment.}"
host_os="$(uname -s)"
host_architecture="$(uname -m)"
case "$host_os:$host_architecture" in
    Linux:x86_64|Linux:amd64)
        ;;
    Darwin:arm64)
        /usr/bin/arch -x86_64 /usr/bin/true || {
            printf 'ERROR: Rosetta 2 is required for pinned Intel-only macOS tools.\n' >&2
            exit 1
        }
        xcrun --find clang >/dev/null || {
            printf 'ERROR: Apple Command Line Tools are required for source builds.\n' >&2
            exit 1
        }
        xcrun --find install_name_tool >/dev/null || {
            printf 'ERROR: Apple install_name_tool is required to relocate pinned Intel tool libraries.\n' >&2
            exit 1
        }
        xcrun --find codesign >/dev/null || {
            printf 'ERROR: Apple codesign is required to sign relocated Intel tool binaries.\n' >&2
            exit 1
        }
        ;;
    *)
        printf 'ERROR: the complete tool stack supports Apple-silicon macOS or glibc Linux x86-64.\n' >&2
        exit 1
        ;;
esac
if [[ "$host_os" == "Linux" ]] && ldd --version 2>&1 | grep -qi musl; then
    printf 'ERROR: Alpine/musl Linux is not supported.\n' >&2
    exit 1
fi

for executable in curl git install make cmake tar unzip python; do
    command -v "$executable" >/dev/null 2>&1 || {
        printf 'ERROR: required build command is unavailable: %s\n' \
            "$executable" >&2
        exit 1
    }
done

temporary_root="$(mktemp -d "${TMPDIR:-/tmp}/postgwas-tools.XXXXXX")"
cleanup() {
    rm -rf -- "$temporary_root"
}
trap cleanup EXIT

verify_sha256() {
    local expected="$1"
    local path="$2"
    local observed
    if command -v sha256sum >/dev/null 2>&1; then
        observed="$(sha256sum "$path" | awk '{print $1}')"
    elif command -v shasum >/dev/null 2>&1; then
        observed="$(shasum -a 256 "$path" | awk '{print $1}')"
    else
        printf 'ERROR: sha256sum or shasum is required to verify downloads.\n' >&2
        exit 1
    fi
    if [[ "$observed" != "$expected" ]]; then
        printf 'ERROR: SHA-256 mismatch for %s; expected %s, observed %s.\n' \
            "$path" "$expected" "$observed" >&2
        exit 1
    fi
}

download() {
    local url="$1"
    local expected="$2"
    local destination="$3"
    curl --fail --location --silent --show-error "$url" -o "$destination"
    verify_sha256 "$expected" "$destination"
}

download_ghcr_blob() {
    local repository="$1"
    local sha256="$2"
    local destination="$3"
    local token
    token="$(
        curl --fail --location --silent --show-error \
            "https://ghcr.io/token?scope=repository:${repository}:pull" | \
            python -c 'import json,sys; print(json.load(sys.stdin)["token"])'
    )"
    curl --fail --location --silent --show-error \
        --header "Authorization: Bearer $token" \
        "https://ghcr.io/v2/${repository}/blobs/sha256:${sha256}" \
        -o "$destination"
    verify_sha256 "$sha256" "$destination"
}

install_binary_archives() {
    local install_root="$temporary_root/binaries"
    local linux_tools_root="$install_prefix/share/postgwas/software/linux-tools"
    mkdir -p "$install_root/magma" "$install_root/finemap" "$install_root/ldstore"

    if [[ "$host_os" == "Linux" ]]; then
        download \
            "$POSTGWAS_MAGMA_URL" "$POSTGWAS_MAGMA_SHA256" \
            "$install_root/magma/magma.zip"
        unzip -q "$install_root/magma/magma.zip" -d "$install_root/magma"
        install -m 0755 "$install_root/magma/magma" "$install_prefix/bin/magma"

        download \
            "$POSTGWAS_FINEMAP_URL" "$POSTGWAS_FINEMAP_SHA256" \
            "$install_root/finemap/finemap.tgz"
        tar -xzf "$install_root/finemap/finemap.tgz" -C "$install_root/finemap"
        mkdir -p "$linux_tools_root/bin"
        install -m 0755 \
            "$install_root/finemap/finemap_v${POSTGWAS_FINEMAP_VERSION}_x86_64/finemap_v${POSTGWAS_FINEMAP_VERSION}_x86_64" \
            "$linux_tools_root/bin/finemap"

        download \
            "$POSTGWAS_LDSTORE_URL" "$POSTGWAS_LDSTORE_SHA256" \
            "$install_root/ldstore/ldstore.tgz"
        tar -xzf "$install_root/ldstore/ldstore.tgz" -C "$install_root/ldstore"
        install -m 0755 \
            "$install_root/ldstore/ldstore_v${POSTGWAS_LDSTORE_VERSION}_x86_64/ldstore_v${POSTGWAS_LDSTORE_VERSION}_x86_64" \
            "$linux_tools_root/bin/ldstore"
        for command_name in finemap ldstore; do
            install -m 0755 \
                "$script_directory/linux_tool_wrapper.sh" \
                "$install_prefix/bin/$command_name"
        done
        return
    fi

    install_macos_binary_archives "$install_root"
}

install_macos_binary_archives() {
    local install_root="$1"
    local software_root="$install_prefix/share/postgwas/software/intel-tools"
    local gcta_root="$install_prefix/share/postgwas/software/gcta"
    local gcc_archive="$install_root/gcc10.tar.gz"
    local gcc_extraction="$install_root/gcc10"

    download \
        "$POSTGWAS_MAGMA_MACOS_URL" "$POSTGWAS_MAGMA_MACOS_SHA256" \
        "$install_root/magma/magma.zip"
    unzip -q "$install_root/magma/magma.zip" -d "$install_root/magma"
    install -m 0755 "$install_root/magma/magma" "$install_prefix/bin/magma"

    download \
        "$POSTGWAS_FINEMAP_MACOS_URL" "$POSTGWAS_FINEMAP_MACOS_SHA256" \
        "$install_root/finemap/finemap.tgz"
    tar -xzf "$install_root/finemap/finemap.tgz" -C "$install_root/finemap"
    download \
        "$POSTGWAS_LDSTORE_MACOS_URL" "$POSTGWAS_LDSTORE_MACOS_SHA256" \
        "$install_root/ldstore/ldstore.tgz"
    tar -xzf "$install_root/ldstore/ldstore.tgz" -C "$install_root/ldstore"

    mkdir -p "$software_root/bin" "$gcc_extraction"
    install -m 0755 \
        "$install_root/finemap/finemap_v${POSTGWAS_FINEMAP_VERSION}_MacOSX/finemap_v${POSTGWAS_FINEMAP_VERSION}_MacOSX" \
        "$software_root/bin/finemap"
    install -m 0755 \
        "$install_root/ldstore/ldstore_v${POSTGWAS_LDSTORE_VERSION}_MacOSX/ldstore_v${POSTGWAS_LDSTORE_VERSION}_MacOSX" \
        "$software_root/bin/ldstore"
    download_ghcr_blob \
        "$POSTGWAS_GCC10_MACOS_REPOSITORY" \
        "$POSTGWAS_GCC10_MACOS_SHA256" \
        "$gcc_archive"
    tar -xzf "$gcc_archive" -C "$gcc_extraction"
    mkdir -p "$software_root/gcc"
    cp -R \
        "$gcc_extraction/$POSTGWAS_GCC10_MACOS_PREFIX/$POSTGWAS_GCC10_MACOS_LIBRARY_SUBDIRECTORY" \
        "$software_root/gcc/runtime"
    for command_name in finemap ldstore; do
        install -m 0755 \
            "$script_directory/intel_tool_wrapper.sh" \
            "$install_prefix/bin/$command_name"
    done

    download \
        "$POSTGWAS_PLINK2_MACOS_URL" "$POSTGWAS_PLINK2_MACOS_SHA256" \
        "$install_root/plink2.zip"
    unzip -q "$install_root/plink2.zip" -d "$install_root/plink2"
    install -m 0755 "$install_root/plink2/plink2" "$install_prefix/bin/plink2"

    download \
        "$POSTGWAS_GCTA_MACOS_URL" "$POSTGWAS_GCTA_MACOS_SHA256" \
        "$install_root/gcta.zip"
    unzip -q "$install_root/gcta.zip" -d "$install_root/gcta"
    mkdir -p "$gcta_root"
    cp -R \
        "$install_root/gcta/gcta-${POSTGWAS_GCTA_MACOS_VERSION}-macOS-arm64/." \
        "$gcta_root/"
    ln -s "$gcta_root/bin/gcta64" "$install_prefix/bin/gcta64"
}

checkout_commit() {
    local repository="$1"
    local commit="$2"
    local destination="$3"
    mkdir -p "$destination"
    git -C "$destination" init --quiet
    git -C "$destination" remote add origin "$repository"
    git -C "$destination" fetch --quiet --depth 1 origin "$commit"
    git -C "$destination" checkout --quiet --detach FETCH_HEAD
    if [[ "$(git -C "$destination" rev-parse HEAD)" != "$commit" ]]; then
        printf 'ERROR: downloaded commit does not match %s.\n' "$commit" >&2
        exit 1
    fi
    rm -rf -- "$destination/.git"
}

create_isolated_environments() {
    local environments_root="$install_prefix/share/postgwas/environments"
    local software_root="$install_prefix/share/postgwas/software"
    local ldsc_environment="$environments_root/ldsc"
    local enrichment_environment="$environments_root/enrichment"
    local vep_environment="$environments_root/vep"
    local intel_tools_environment="$environments_root/intel-tools"
    local mixer_environment="$environments_root/mixer"
    local ldsc_source="$software_root/ldsc/$POSTGWAS_LDSC_COMMIT"
    local mixer_matplotlib_venn_archive="${temporary_root}/matplotlib-venn-${POSTGWAS_MIXER_MATPLOTLIB_VENN_VERSION}.zip"
    local environment_platform=()

    if [[ "$host_os" == "Darwin" ]]; then
        environment_platform=(--platform osx-64)
    fi

    mkdir -p "$environments_root" "$(dirname "$ldsc_source")"
    "$mamba_executable" env create --yes "${environment_platform[@]}" \
        --prefix "$ldsc_environment" \
        --file "$script_directory/environments/ldsc.yml"
    checkout_commit \
        "$POSTGWAS_LDSC_REPOSITORY" "$POSTGWAS_LDSC_COMMIT" "$ldsc_source"
    "$mamba_executable" run --prefix "$ldsc_environment" \
        python -m pip install --no-deps --no-build-isolation "$ldsc_source"
    "$mamba_executable" run --prefix "$ldsc_environment" python -m pip check
    ln -s "$ldsc_environment/bin/ldsc.py" "$install_prefix/bin/ldsc.py"
    ln -s \
        "$ldsc_environment/bin/munge_sumstats.py" \
        "$install_prefix/bin/munge_sumstats.py"

    "$mamba_executable" env create --yes "${environment_platform[@]}" \
        --prefix "$enrichment_environment" \
        --file "$script_directory/environments/enrichment.yml"
    "$mamba_executable" run --prefix "$enrichment_environment" \
        python -m pip install --no-deps --no-build-isolation "$repository_root"
    "$mamba_executable" run --prefix "$enrichment_environment" python -m pip check

    "$mamba_executable" create --yes "${environment_platform[@]}" \
        --prefix "$vep_environment" \
        --channel conda-forge --channel bioconda --override-channels \
        "ensembl-vep=$POSTGWAS_VEP_VERSION"
    install -m 0755 \
        "$script_directory/vep_wrapper.sh" \
        "$install_prefix/bin/vep"

    if [[ "$host_os" == "Darwin" ]]; then
        "$mamba_executable" env create --yes --platform osx-64 \
            --prefix "$intel_tools_environment" \
            --file "$script_directory/environments/intel-tools-macos.yml"
        "$mamba_executable" env create --yes --platform osx-64 \
            --prefix "$mixer_environment" \
            --file "$script_directory/environments/mixer-macos.yml"
        download \
            "$POSTGWAS_MIXER_MATPLOTLIB_VENN_URL" \
            "$POSTGWAS_MIXER_MATPLOTLIB_VENN_SHA256" \
            "$mixer_matplotlib_venn_archive"
        "$mamba_executable" run --prefix "$mixer_environment" \
            python -m pip install --no-deps --no-build-isolation \
            "$mixer_matplotlib_venn_archive"
        "$mamba_executable" run --prefix "$mixer_environment" \
            python -m pip check
    fi

    mkdir -p "$install_prefix/etc/conda/activate.d"
    printf "export POSTGWAS_ENRICHMENT_PYTHON='%s/bin/python'\n" \
        "$enrichment_environment" \
        >"$install_prefix/etc/conda/activate.d/postgwas-enrichment.sh"
}

relocate_macos_dependency() {
    local install_name_tool_executable="$1"
    local binary_path="$2"
    local original_path="$3"
    local relocated_path="$4"

    if ! otool -L "$binary_path" | awk '{print $1}' | \
        grep -F -x -q "$original_path"; then
        printf 'ERROR: pinned macOS dependency was not found in %s: %s\n' \
            "$binary_path" "$original_path" >&2
        exit 1
    fi
    "$install_name_tool_executable" \
        -change "$original_path" "$relocated_path" "$binary_path"
}

relocate_macos_library_id() {
    local install_name_tool_executable="$1"
    local library_path="$2"
    local expected_id="$3"
    local current_id

    current_id="$(otool -D "$library_path" | tail -n 1)"
    if [[ "$current_id" != "$expected_id" ]]; then
        printf 'ERROR: unexpected pinned macOS library ID in %s: %s\n' \
            "$library_path" "$current_id" >&2
        exit 1
    fi
    "$install_name_tool_executable" \
        -id "@rpath/$(basename "$library_path")" "$library_path"
}

relocate_macos_intel_tools() {
    [[ "$host_os" == "Darwin" ]] || return 0

    local software_root="$install_prefix/share/postgwas/software/intel-tools"
    local gcc_runtime="$software_root/gcc/runtime"
    local finemap_binary="$software_root/bin/finemap"
    local ldstore_binary="$software_root/bin/ldstore"
    local install_name_tool_executable
    local codesign_executable
    local binary_path
    local library_name
    local library_path
    install_name_tool_executable="$(xcrun --find install_name_tool)"
    codesign_executable="$(xcrun --find codesign)"

    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/lib/libzstd.1.dylib \
        @rpath/libzstd.1.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/opt/gcc/lib/gcc/10/libgfortran.5.dylib \
        @rpath/libgfortran.5.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/opt/gcc/lib/gcc/10/libquadmath.0.dylib \
        @rpath/libquadmath.0.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/opt/gcc/lib/gcc/10/libstdc++.6.dylib \
        @rpath/libstdc++.6.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/opt/libomp/lib/libomp.dylib \
        @rpath/libomp.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$finemap_binary" /usr/local/lib/gcc/10/libgcc_s.1.dylib \
        @rpath/libgcc_s.1.dylib

    relocate_macos_dependency "$install_name_tool_executable" \
        "$ldstore_binary" /usr/local/lib/libzstd.1.dylib \
        @rpath/libzstd.1.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$ldstore_binary" /usr/local/opt/openblas/lib/libopenblas.0.dylib \
        @rpath/libopenblas.0.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$ldstore_binary" /usr/local/opt/gcc/lib/gcc/9/libstdc++.6.dylib \
        @rpath/libstdc++.6.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$ldstore_binary" /usr/local/opt/gcc/lib/gcc/9/libgomp.1.dylib \
        @rpath/libgomp.1.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$ldstore_binary" /usr/local/lib/gcc/9/libgcc_s.1.dylib \
        @rpath/libgcc_s.1.dylib

    relocate_macos_dependency "$install_name_tool_executable" \
        "$gcc_runtime/libgfortran.5.dylib" \
        '@@HOMEBREW_CELLAR@@/gcc@10/10.5.0/lib/gcc/10/libquadmath.0.dylib' \
        @rpath/libquadmath.0.dylib
    relocate_macos_dependency "$install_name_tool_executable" \
        "$gcc_runtime/libgcc_s.1.dylib" \
        '@@HOMEBREW_CELLAR@@/gcc@10/10.5.0/lib/gcc/10/libgcc_s.1.1.dylib' \
        @rpath/libgcc_s.1.1.dylib

    for binary_path in "$finemap_binary" "$ldstore_binary"; do
        "$install_name_tool_executable" \
            -add_rpath '@loader_path/../gcc/runtime' "$binary_path"
        "$install_name_tool_executable" \
            -add_rpath '@loader_path/../../../environments/intel-tools/lib' \
            "$binary_path"
    done

    for library_name in \
        libgfortran.5.dylib \
        libquadmath.0.dylib \
        libstdc++.6.dylib \
        libgcc_s.1.dylib \
        libgcc_s.1.1.dylib \
        libgomp.1.dylib; do
        library_path="$gcc_runtime/$library_name"
        relocate_macos_library_id "$install_name_tool_executable" \
            "$library_path" \
            "@@HOMEBREW_PREFIX@@/opt/gcc@10/lib/gcc/10/$library_name"
        "$codesign_executable" --force --sign - "$library_path"
    done
    "$codesign_executable" --force --sign - "$finemap_binary"
    "$codesign_executable" --force --sign - "$ldstore_binary"
}

install_bgenix_macos() {
    local archive="$temporary_root/bgenix.tar.gz"
    local source_root="$temporary_root/bgenix-${POSTGWAS_BGENIX_VERSION}"
    local wscript="$source_root/wscript"

    [[ "$host_os" == "Darwin" ]] || return 0
    download \
        "$POSTGWAS_BGENIX_SOURCE_URL" "$POSTGWAS_BGENIX_SOURCE_SHA256" \
        "$archive"
    tar -xzf "$archive" -C "$temporary_root"
    python - "$wscript" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
needle = (
    "cfg.env.CXXFLAGS = [ '-Wall', '-pedantic', "
    "'-Wno-unused-local-typedefs', '-Wno-c++11-long-long', "
    "'-Wno-deprecated-declarations', '-Wno-long-long', '-fPIC' ]"
)
replacement = needle[:-2] + (
    ", '-D_LIBCPP_ENABLE_CXX17_REMOVED_AUTO_PTR', "
    "'-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION' ]"
)
if text.count(needle) != 1:
    raise SystemExit("Unexpected pinned bgenix compiler configuration")
path.write_text(text.replace(needle, replacement), encoding="utf-8")
PY
    (
        cd "$source_root"
        ./waf configure --prefix="$install_prefix"
        ./waf -j"$build_jobs"
        ./build/test/unit/test_bgen
        ./waf install
    )
    bgenix -help >/dev/null
}

install_kpops_and_caldera() {
    local resource_root="$install_prefix/share/postgwas/upstream"
    mkdir -p "$resource_root"
    bash "$repository_root/tools/resource_preparation/prepare_kpops_caldera_resources.sh" \
        --resource-root "$resource_root" \
        --python "$install_prefix/bin/python" \
        --software-only
}

install_mixer() {
    if [[ "$host_os" == "Darwin" ]]; then
        install_mixer_macos
        return
    fi

    local build_root="$temporary_root/mixer"
    local boost_archive="$build_root/boost.tar.gz"
    local boost_root="$build_root/boost_1_69_0"
    local source_root="$build_root/gsa-mixer"
    local mixer_root="$install_prefix/share/postgwas/mixer"
    mkdir -p "$build_root"

    download \
        "$POSTGWAS_MIXER_BOOST_URL" "$POSTGWAS_MIXER_BOOST_SHA256" \
        "$boost_archive"
    tar -xzf "$boost_archive" -C "$build_root"
    (
        cd "$boost_root"
        ./bootstrap.sh \
            --with-libraries=program_options,filesystem,system,date_time
        ./b2 \
            -j"$build_jobs" \
            cxxflags=-fPIC \
            variant=release \
            link=static \
            threading=multi \
            --with-program_options \
            --with-filesystem \
            --with-system \
            --with-date_time
    )

    checkout_commit \
        "$POSTGWAS_MIXER_REPOSITORY" "$POSTGWAS_MIXER_COMMIT" "$source_root"
    sed -i 's/-march=native/-mssse3 -msse4.1/g' "$source_root/src/CMakeLists.txt"
    local bitutil="$source_root/src/TurboPFor/bitutil.c"
    local nlopt_stop="$source_root/src/nlopt/stop.c"
    [[ "$(grep -c 'const _md = _md_;' "$bitutil")" -eq 2 ]]
    sed -i 's/const _md = _md_;/const _t_ _md = _md_;/g' "$bitutil"
    [[ "$(grep -c 'const _t_ _md = _md_;' "$bitutil")" -eq 2 ]]
    [[ "$(grep -c '^#include <stdarg.h>$' "$nlopt_stop")" -eq 1 ]]
    [[ "$(grep -c '^#include <stdlib.h>$' "$nlopt_stop")" -eq 0 ]]
    sed -i '/^#include <stdarg.h>$/a #include <stdlib.h>' "$nlopt_stop"
    [[ "$(grep -c '^#include <stdlib.h>$' "$nlopt_stop")" -eq 1 ]]

    cmake \
        -S "$source_root/src" \
        -B "$source_root/src/build" \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DBoost_NO_BOOST_CMAKE=ON \
        -DBOOST_ROOT="$boost_root" \
        -DBoost_USE_STATIC_LIBS=ON
    cmake \
        --build "$source_root/src/build" \
        --target bgmg \
        --parallel "$build_jobs"

    mkdir -p "$mixer_root/precimed" "$mixer_root/lib"
    cp -R "$source_root/precimed/." "$mixer_root/precimed/"
    install -m 0755 "$source_root/src/build/lib/libbgmg.so" "$mixer_root/lib/libbgmg.so"
    install -m 0644 "$source_root/LICENSE" "$mixer_root/LICENSE"
    rm -rf -- "$mixer_root/precimed/mixer-test" "$mixer_root/precimed/__pycache__"
    find "$mixer_root/precimed" -type f \
        \( -name '*.pyc' -o -name '*.pyo' \) -delete
    for command_name in mixer.py mixer_dev.py mixer_figures.py; do
        install -m 0755 \
            "$script_directory/mixer_wrapper.py" \
            "$install_prefix/bin/$command_name"
    done
}

install_mixer_macos() {
    local build_root="$temporary_root/mixer"
    local boost_archive="$build_root/boost.tar.gz"
    local boost_root="$build_root/boost_1_69_0"
    local source_root="$build_root/gsa-mixer"
    local source_build="$source_root/src/build"
    local mixer_root="$install_prefix/share/postgwas/mixer"
    local mixer_environment="$install_prefix/share/postgwas/environments/mixer"
    local compile_flags="-arch x86_64 -I${mixer_environment}/include"
    local linker_flags="-arch x86_64 -L${mixer_environment}/lib -lomp -Wl,-rpath,${mixer_environment}/lib"
    mkdir -p "$build_root"

    download \
        "$POSTGWAS_MIXER_BOOST_URL" "$POSTGWAS_MIXER_BOOST_SHA256" \
        "$boost_archive"
    tar -xzf "$boost_archive" -C "$build_root"
    checkout_commit \
        "$POSTGWAS_MIXER_REPOSITORY" "$POSTGWAS_MIXER_COMMIT" "$source_root"

    python - "$boost_root" "$source_root" <<'PY'
from pathlib import Path
import sys

boost_root = Path(sys.argv[1])
source_root = Path(sys.argv[2])

def replace(path: Path, old: str, new: str, count: int = 1) -> None:
    text = path.read_text(encoding="utf-8")
    if text.count(old) != count:
        raise SystemExit(
            "Unexpected pinned source while patching %s: %r" % (path, old)
        )
    path.write_text(text.replace(old, new), encoding="utf-8")

replace(
    boost_root / "tools/build/src/tools/darwin.jam",
    "        flags darwin.compile.c++ OPTIONS $(condition) : -fcoalesce-templates ;",
    "        # Modern Apple Clang removed -fcoalesce-templates.",
)
replace(
    boost_root / "tools/build/src/tools/darwin.jam",
    "        flags darwin.compile OPTIONS $(condition) : -Wno-long-double ;",
    "        # Modern Apple Clang removed -Wno-long-double.",
)
replace(
    boost_root / "boost/mpl/aux_/integral_wrapper.hpp",
    "#if BOOST_WORKAROUND(__EDG_VERSION__, <= 243)",
    "#if BOOST_WORKAROUND(__EDG_VERSION__, <= 243) \\\n    || __cplusplus >= 201103L",
)
replace(
    source_root / "src/CMakeLists.txt",
    "-march=native -mfpmath=sse",
    "-mssse3 -msse4.1",
    count=2,
)
replace(
    source_root / "src/TurboPFor/bitutil.c",
    "const _md = _md_;",
    "const _t_ _md = _md_;",
    count=2,
)
replace(
    source_root / "src/nlopt/stop.c",
    "#include <stdarg.h>\n",
    "#include <stdarg.h>\n#include <stdlib.h>\n",
)
replace(
    source_root / "src/zlib/zutil.h",
    "#if defined(MACOS) || defined(TARGET_OS_MAC)",
    "#if defined(MACOS)",
)
replace(
    source_root / "src/strict_fstream.hpp",
    "#elif (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE",
    "#elif defined(__APPLE__) || \\\n    ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && ! _GNU_SOURCE)",
)
PY

    "$mamba_executable" run --prefix "$mixer_environment" \
        bash -c \
        'cd "$1" && ./bootstrap.sh --with-libraries=program_options,filesystem,system,date_time' \
        bash "$boost_root"
    python - "$boost_root/project-config.jam" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
old = "    using darwin ; "
new = "    using darwin : : /usr/bin/clang++ ;"
if text.count(old) != 1:
    raise SystemExit("Unexpected Boost project-config.jam compiler declaration")
path.write_text(text.replace(old, new), encoding="utf-8")
PY
    "$mamba_executable" run --prefix "$mixer_environment" \
        bash -c \
        'cd "$1" && ./b2 -j"$2" cxxflags="-fPIC -arch x86_64" linkflags="-arch x86_64" architecture=x86 address-model=64 variant=release link=static threading=multi --with-program_options --with-filesystem --with-system --with-date_time' \
        bash "$boost_root" "$build_jobs"

    "$mamba_executable" run --prefix "$mixer_environment" \
        env \
        CC=/usr/bin/clang \
        CXX=/usr/bin/clang++ \
        CFLAGS="$compile_flags" \
        CXXFLAGS="$compile_flags" \
        cmake \
        -S "$source_root/src" \
        -B "$source_build" \
        -DCMAKE_OSX_ARCHITECTURES=x86_64 \
        -DCMAKE_SHARED_LINKER_FLAGS="$linker_flags" \
        -DCMAKE_EXE_LINKER_FLAGS="$linker_flags" \
        -DBoost_NO_SYSTEM_PATHS=ON \
        -DBOOST_ROOT="$boost_root" \
        -DBoost_LIBRARY_DIR="$boost_root/stage/lib"
    "$mamba_executable" run --prefix "$mixer_environment" \
        cmake --build "$source_build" \
        --target bgmg bgmg-test --parallel "$build_jobs"

    DYLD_LIBRARY_PATH="$mixer_environment/lib" \
        /usr/bin/arch -x86_64 "$source_build/bin/bgmg-test" \
        --gtest_filter='BgmgTest.CalcLikelihood:TestLd.SingleMarker:UgmgTest.CalcUnifiedGaussianLikelihood'

    mkdir -p "$mixer_root/precimed" "$mixer_root/lib"
    cp -R "$source_root/precimed/." "$mixer_root/precimed/"
    install -m 0755 "$source_build/lib/libbgmg.dylib" \
        "$mixer_root/lib/libbgmg.dylib"
    install -m 0644 "$source_root/LICENSE" "$mixer_root/LICENSE"
    rm -rf -- "$mixer_root/precimed/mixer-test" "$mixer_root/precimed/__pycache__"
    find "$mixer_root/precimed" -type f \
        \( -name '*.pyc' -o -name '*.pyo' \) -delete
    for command_name in mixer.py mixer_dev.py mixer_figures.py; do
        install -m 0755 \
            "$script_directory/mixer_wrapper.py" \
            "$install_prefix/bin/$command_name"
    done
}

printf 'Installing pinned external scientific programs into %s\n' "$install_prefix"
install_binary_archives
create_isolated_environments
relocate_macos_intel_tools
install_bgenix_macos
install_kpops_and_caldera
install_mixer
printf 'Installed external scientific programs successfully.\n'
