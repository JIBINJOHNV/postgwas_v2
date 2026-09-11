#!/usr/bin/env bash
set -Eeuo pipefail

fail() {
    printf 'ERROR: %s\n' "$1" >&2
    exit 1
}

report_failed_check() {
    local status="$1"
    local line_number="$2"
    local command_text="$3"
    trap - ERR
    printf 'ERROR: software verification failed at line %s (exit %s): %s\n' \
        "$line_number" "$status" "$command_text" >&2
    exit "$status"
}

trap 'report_failed_check "$?" "$LINENO" "$BASH_COMMAND"' ERR

require_command() {
    local command_name="$1"
    command -v "$command_name" >/dev/null 2>&1 || \
        fail "required command is unavailable: $command_name"
}

verify_gcta() {
    local temporary_directory
    local version_output
    temporary_directory="$(mktemp -d)"
    version_output="$(
        gcta64 --out "${temporary_directory}/verification" 2>&1 || true
    )"
    rm -rf -- "$temporary_directory"
    printf '%s\n' "$version_output" | \
        grep -Eq '^\* version v[0-9]+([.][0-9]+)+'
}

printf 'Verifying PostGWAS container software stack\n'

required_commands=(
    python postgwas Rscript
    bcftools tabix bedtools pigz
    plink plink2 gcta64
    magma finemap ldstore bgenix
    ldsc.py munge_sumstats.py scdrs k-pops.py
    micromamba
)
for executable in "${required_commands[@]}"; do
    require_command "$executable"
done

python -m pip check
postgwas --help >/dev/null
postgwas config --help >/dev/null
python -c \
    "import postgwas, numpy, pandas, scipy, pysam, sklearn, torch, xgboost"

Rscript -e \
    "for (p in c('data.table','dplyr','ggplot2','susieR')) stopifnot(requireNamespace(p, quietly=TRUE))" \
    >/dev/null

bcftools --version >/dev/null
bcftools plugin -l | grep -qx liftover
tabix --version >/dev/null 2>&1
plink --version >/dev/null
plink2 --version >/dev/null
verify_gcta
magma --version >/dev/null 2>&1
python "$(command -v k-pops.py)" --help >/dev/null

micromamba run -n ldsc python -m pip check
ldsc.py --help >/dev/null
munge_sumstats.py --help >/dev/null

micromamba run -n enricher python -m pip check
micromamba run -n enricher python -c \
    "import gprofiler, gseapy, networkx, omnipath, postgwas.modules.enrichment.service, rpy2, zeep"
micromamba run -n enricher Rscript -e \
    "stopifnot(requireNamespace('WebGestaltR', quietly=TRUE))" \
    >/dev/null

micromamba run -n vep vep --help >/dev/null 2>&1

caldera_root="/opt/conda/envs/postgwas/share/postgwas/caldera"
test -s "$caldera_root/z_caldera.R" || \
    fail "CALDERA runner is unavailable"
test -s "$caldera_root/trained_models/caldera_model_no_covs.rds" || \
    fail "CALDERA trained model is unavailable"
Rscript -e "parse(file='$caldera_root/z_caldera.R')" >/dev/null

test -s /tools/mixer/lib/libbgmg.so || fail "MiXeR native library is unavailable"
if ldd /tools/mixer/lib/libbgmg.so | grep -q 'not found'; then
    ldd /tools/mixer/lib/libbgmg.so >&2
    fail "MiXeR native library has unresolved dependencies"
fi
python /tools/mixer/precimed/mixer.py --version >/dev/null
python /tools/mixer/precimed/mixer.py --help >/dev/null
python /tools/mixer/precimed/mixer_figures.py --help >/dev/null

printf 'All PostGWAS container software checks passed.\n'
printf 'Scientific references, VEP cache data, and service credentials must be supplied per analysis.\n'
