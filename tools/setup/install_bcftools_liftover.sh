#!/usr/bin/env bash
set -euo pipefail

bcftools_version="${1:-1.23.1}"
score_revision="${2:-909d23019e19aeadf3bf6fe1407fd6afc094592a}"
install_prefix="${CONDA_PREFIX:?Activate the target conda environment first.}"
build_directory="$(mktemp -d /tmp/postgwas-liftover.XXXXXX)"

cleanup() {
    rm -rf -- "$build_directory"
}
trap cleanup EXIT

archive="$build_directory/bcftools.tar.bz2"
source_directory="$build_directory/bcftools-$bcftools_version"
curl -fsSL \
    "https://github.com/samtools/bcftools/releases/download/$bcftools_version/bcftools-$bcftools_version.tar.bz2" \
    -o "$archive"
tar -xjf "$archive" -C "$build_directory"
curl -fsSL \
    "https://raw.githubusercontent.com/freeseek/score/$score_revision/liftover.c" \
    -o "$source_directory/plugins/liftover.c"

cd "$source_directory"
./configure --prefix="$install_prefix" --with-htslib="$install_prefix"
make plugins/liftover.so
install -d "$install_prefix/libexec/bcftools"
install -m 755 plugins/liftover.so "$install_prefix/libexec/bcftools/liftover.so"

bcftools plugin -l | grep -qx liftover
printf 'Installed bcftools liftover plugin in %s\n' "$install_prefix"
