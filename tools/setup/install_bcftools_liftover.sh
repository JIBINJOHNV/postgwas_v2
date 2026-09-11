#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools/setup/software_versions.env
source "$script_directory/software_versions.env"

default_bcftools_version="$POSTGWAS_BCFTOOLS_VERSION"
default_score_revision="$POSTGWAS_SCORE_COMMIT"
default_bcftools_sha256="$POSTGWAS_BCFTOOLS_SHA256"
default_liftover_sha256="$POSTGWAS_LIFTOVER_SHA256"

bcftools_version="${1:-$default_bcftools_version}"
score_revision="${2:-$default_score_revision}"
bcftools_sha256="${3:-}"
liftover_sha256="${4:-}"
build_jobs="${POSTGWAS_INSTALL_JOBS:-2}"

if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: POSTGWAS_INSTALL_JOBS must be a positive integer; found %s.\n' \
        "$build_jobs" >&2
    exit 2
fi

if [[ -z "$bcftools_sha256" ]]; then
    if [[ "$bcftools_version" != "$default_bcftools_version" ]]; then
        printf 'ERROR: provide the expected bcftools archive SHA-256 as argument 3 when overriding version %s.\n' \
            "$bcftools_version" >&2
        exit 2
    fi
    bcftools_sha256="$default_bcftools_sha256"
fi
if [[ -z "$liftover_sha256" ]]; then
    if [[ "$score_revision" != "$default_score_revision" ]]; then
        printf 'ERROR: provide the expected liftover.c SHA-256 as argument 4 when overriding score revision %s.\n' \
            "$score_revision" >&2
        exit 2
    fi
    liftover_sha256="$default_liftover_sha256"
fi

install_prefix="${CONDA_PREFIX:?Activate the target conda environment first.}"
build_directory="$(mktemp -d /tmp/postgwas-liftover.XXXXXX)"

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

cleanup() {
    rm -rf -- "$build_directory"
}
trap cleanup EXIT

archive="$build_directory/bcftools.tar.bz2"
source_directory="$build_directory/bcftools-$bcftools_version"
curl -fsSL \
    "https://github.com/samtools/bcftools/releases/download/$bcftools_version/bcftools-$bcftools_version.tar.bz2" \
    -o "$archive"
verify_sha256 "$bcftools_sha256" "$archive"
tar -xjf "$archive" -C "$build_directory"
curl -fsSL \
    "https://raw.githubusercontent.com/freeseek/score/$score_revision/liftover.c" \
    -o "$source_directory/plugins/liftover.c"
verify_sha256 "$liftover_sha256" "$source_directory/plugins/liftover.c"

cd "$source_directory"
./configure --prefix="$install_prefix" --with-htslib="$install_prefix"
make -j "$build_jobs" bcftools
make plugins/liftover.so
install -d "$install_prefix/bin" "$install_prefix/share/man/man1"
install -m 755 bcftools "$install_prefix/bin/bcftools"
install -m 755 misc/* "$install_prefix/bin/"
install -m 644 doc/bcftools.1 "$install_prefix/share/man/man1/bcftools.1"
install -d "$install_prefix/libexec/bcftools"
install -m 755 plugins/liftover.so "$install_prefix/libexec/bcftools/liftover.so"

installed_bcftools_version="$(bcftools --version-only)"
installed_bcftools_version="${installed_bcftools_version%%+*}"
if [[ "$installed_bcftools_version" != "$bcftools_version" ]]; then
    printf 'ERROR: installed bcftools version does not match %s.\n' \
        "$bcftools_version" >&2
    exit 1
fi
bcftools plugin -l | grep -qx liftover
printf 'Installed bcftools %s and its liftover plugin in %s\n' \
    "$bcftools_version" "$install_prefix"
