#!/usr/bin/env bash
set -Eeuo pipefail

export PYTHONNOUSERSITE=1
unset PYTHONPATH

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd -P)"
# shellcheck source=tools/setup/software_versions.env
source "$script_directory/software_versions.env"

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
    local command_path
    command_path="$(command -v "$command_name")" || \
        fail "required command is unavailable: $command_name"
    case "$command_path" in
        "$install_prefix"/*)
            ;;
        *)
            fail "required command resolved outside the active environment: $command_name -> $command_path"
            ;;
    esac
}

verify_postgwas_import_origins() {
    local python_executable="$1"
    shift
    "$python_executable" - "$@" <<'PY'
from importlib import import_module
from pathlib import Path
import sys

allowed_roots = [Path(value).resolve() for value in sys.argv[1:]]
module_names = (
    "postgwas",
    "postgwas.modules.flames.service",
    "postgwas.modules.pops.service",
)
for module_name in module_names:
    module_path = Path(import_module(module_name).__file__).resolve()
    if not any(module_path.is_relative_to(root) for root in allowed_roots):
        raise SystemExit(
            "%s resolved outside the installed PostGWAS roots: %s"
            % (module_name, module_path)
        )
    print("Verified import origin: %s -> %s" % (module_name, module_path))
PY
}

verify_gcta() {
    local temporary_directory
    local version_output
    temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/postgwas-gcta.XXXXXX")"
    version_output="$(
        gcta64 --out "${temporary_directory}/verification" 2>&1 || true
    )"
    rm -rf -- "$temporary_directory"
    printf '%s\n' "$version_output" | \
        grep -Eq '^\* version v[0-9]+([.][0-9]+)+'
}

install_prefix="${CONDA_PREFIX:-$(python -c 'import sys; print(sys.prefix)')}"
host_os="$(uname -s)"
ldsc_environment="$install_prefix/share/postgwas/environments/ldsc"
enrichment_environment="$install_prefix/share/postgwas/environments/enrichment"
vep_environment="$install_prefix/share/postgwas/environments/vep"
mixer_environment="$install_prefix/share/postgwas/environments/mixer"
caldera_root="$install_prefix/share/postgwas/caldera"
mixer_root="$install_prefix/share/postgwas/mixer"
intel_tools_root="$install_prefix/share/postgwas/software/intel-tools"
linux_tools_root="$install_prefix/share/postgwas/software/linux-tools"

printf 'Verifying the complete PostGWAS Mamba installation\n'

required_commands=(
    python postgwas Rscript
    bcftools tabix bedtools pigz
    plink plink2 gcta64
    magma finemap ldstore bgenix
    ldsc.py munge_sumstats.py scdrs k-pops.py
    mixer.py mixer_figures.py vep
)
for executable in "${required_commands[@]}"; do
    require_command "$executable"
done

python -m pip check
postgwas --help >/dev/null
postgwas config --help >/dev/null
python -c \
    "import postgwas, numpy, pandas, scipy, pysam, sklearn, torch, xgboost; from postgwas.modules.flames import service as flames; from postgwas.modules.pops import service as pops"
python -c \
    "from importlib.resources import files; model = files('postgwas.modules.flames.resources.model'); assert model.joinpath('FLAMES_XGB_model.sav').is_file(); assert model.joinpath('features.txt').read_text(encoding='utf-8').strip()"
verify_postgwas_import_origins python "$install_prefix" "$repository_root"

Rscript -e \
    "for (p in c('data.table','dplyr','ggplot2','susieR','survey','xml2')) stopifnot(requireNamespace(p, quietly=TRUE))" \
    >/dev/null

bcftools_version="$(bcftools --version-only)"
bcftools_version="${bcftools_version%%+*}"
[[ "$bcftools_version" == "$POSTGWAS_BCFTOOLS_VERSION" ]]
bcftools plugin -l | grep -qx liftover
tabix --version >/dev/null 2>&1
plink --version >/dev/null
plink2_version="$(plink2 --version 2>&1)"
if [[ "$host_os" == "Darwin" ]]; then
    printf '%s\n' "$plink2_version" | \
        grep -Fq "PLINK v${POSTGWAS_PLINK2_MACOS_VERSION}"
fi
verify_gcta

magma_version="$(magma --version 2>&1)"
printf '%s\n' "$magma_version" | \
    grep -Eq "MAGMA version: v${POSTGWAS_MAGMA_VERSION}"
if [[ "$host_os" == "Darwin" ]]; then
    file "$(command -v magma)" | grep -q 'x86_64'
    file "$(command -v plink2)" | grep -q 'arm64'
    file "$(command -v gcta64)" | grep -q 'arm64'
    file "$(command -v bgenix)" | grep -q 'arm64'
    file "$vep_environment/bin/perl" | grep -q 'x86_64'
    for binary_path in \
        "$intel_tools_root/bin/finemap" \
        "$intel_tools_root/bin/ldstore"; do
        file "$binary_path" | grep -q 'x86_64'
        if otool -L "$binary_path" | grep -Eq '/usr/local|@@HOMEBREW'; then
            otool -L "$binary_path" >&2
            fail "Intel-mac tool contains an external absolute library path: $binary_path"
        fi
        otool -l "$binary_path" | \
            grep -Fq '@loader_path/../gcc/runtime'
        otool -l "$binary_path" | \
            grep -Fq '@loader_path/../../../environments/intel-tools/lib'
    done
else
    for binary_path in \
        "$linux_tools_root/bin/finemap" \
        "$linux_tools_root/bin/ldstore"; do
        test -x "$binary_path" || \
            fail "installed Linux tool binary is unavailable: $binary_path"
        if LD_LIBRARY_PATH="$install_prefix/lib" ldd "$binary_path" | \
            grep -q 'not found'; then
            LD_LIBRARY_PATH="$install_prefix/lib" ldd "$binary_path" >&2
            fail "Linux tool has unresolved environment dependencies: $binary_path"
        fi
    done
fi
finemap --help >/dev/null
ldstore --help >/dev/null
bgenix -help >/dev/null

"$ldsc_environment/bin/python" -m pip check
ldsc.py --help >/dev/null
munge_sumstats.py --help >/dev/null

(
    export PATH="$enrichment_environment/bin:$PATH"
    export R_HOME="$enrichment_environment/lib/R"
    export R_LIBS_USER=/dev/null
    "$enrichment_environment/bin/python" -m pip check
    "$enrichment_environment/bin/python" -c \
        "import gprofiler, gseapy, networkx, omnipath, postgwas.modules.enrichment.service, rpy2, zeep"
    verify_postgwas_import_origins \
        "$enrichment_environment/bin/python" "$enrichment_environment"
    "$enrichment_environment/bin/Rscript" -e \
        "stopifnot(requireNamespace('WebGestaltR', quietly=TRUE))" \
        >/dev/null
)

vep --help >/dev/null

python "$(command -v k-pops.py)" --help >/dev/null
test -s "$caldera_root/z_caldera.R" || fail "CALDERA runner is unavailable"
test -s "$caldera_root/trained_models/caldera_model_no_covs.rds" || \
    fail "CALDERA trained model is unavailable"
Rscript -e "parse(file='$caldera_root/z_caldera.R')" >/dev/null

if [[ "$host_os" == "Darwin" ]]; then
    mixer_library="$mixer_root/lib/libbgmg.dylib"
    test -s "$mixer_library" || fail "MiXeR native library is unavailable"
    file "$mixer_library" | grep -q 'x86_64'
    otool -L "$mixer_library" | grep -q '@rpath/libomp.dylib'
    "$mixer_environment/bin/python" -m pip check
    "$mixer_environment/bin/python" -c \
        "from importlib.metadata import version; from matplotlib_venn import venn2; assert version('matplotlib-venn') == '$POSTGWAS_MIXER_MATPLOTLIB_VENN_VERSION'; assert callable(venn2)"
else
    mixer_library="$mixer_root/lib/libbgmg.so"
    test -s "$mixer_library" || fail "MiXeR native library is unavailable"
    if ldd "$mixer_library" | grep -q 'not found'; then
        ldd "$mixer_library" >&2
        fail "MiXeR native library has unresolved dependencies"
    fi
fi
mixer.py --version >/dev/null
mixer.py --help >/dev/null
mixer_figures.py --help >/dev/null

test ! -e "$install_prefix/bin/gwas2vcf" || \
    fail "a separate gwas2vcf executable was unexpectedly installed"

printf 'All complete PostGWAS Mamba installation checks passed.\n'
printf 'Scientific references, VEP cache data, and service credentials must be supplied per analysis.\n'
