#!/usr/bin/env bash
set -euo pipefail

# Installation must not inherit modules from a sibling checkout or user site.
export PYTHONNOUSERSITE=1
unset PYTHONPATH

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd)"
environment_file="$repository_root/environment.yml"
editable=false
update=false
all_tools=false
build_jobs="${POSTGWAS_INSTALL_JOBS:-2}"

usage() {
    cat <<'EOF'
Usage: bash tools/setup/install_postgwas.sh [options]

Create a new PostGWAS environment, install the current checkout and its
dependencies, and verify the installation. --all-tools also downloads and
installs the complete scientific-software stack on Apple-silicon macOS or
glibc Linux x86-64.

Options:
  --all-tools       Install all supported software into one user-facing Mamba env.
  --name NAME       Override the environment.yml environment name.
  --jobs N          Use N parallel jobs for source builds (default: 2).
  --editable        Install the checkout in editable development mode.
  --update          Update an existing portable environment instead of refusing.
  -h, --help        Show this help message.

Supported hosts:
  --all-tools: Apple-silicon macOS or glibc Linux x86-64 with Mamba.
  Portable mode: macOS arm64/x86_64 or glibc-based Linux x86-64.

Scientific reference datasets, API credentials, and licensed cache data are
analysis inputs and are not installed. No separate gwas2vcf package or binary
is downloaded.
EOF
}

environment_name="$(awk '$1 == "name:" { print $2; exit }' "$environment_file")"

while (($#)); do
    case "$1" in
        --name)
            if (($# < 2)); then
                printf 'ERROR: --name requires a value.\n' >&2
                exit 2
            fi
            environment_name="$2"
            shift 2
            ;;
        --editable)
            editable=true
            shift
            ;;
        --update)
            update=true
            shift
            ;;
        --all-tools)
            all_tools=true
            shift
            ;;
        --jobs)
            if (($# < 2)); then
                printf 'ERROR: --jobs requires a value.\n' >&2
                exit 2
            fi
            build_jobs="$2"
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

if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'ERROR: --jobs must be a positive integer; found %s.\n' \
        "$build_jobs" >&2
    exit 2
fi
if $all_tools && $update; then
    printf 'ERROR: --all-tools always creates a new environment and cannot be combined with --update.\n' >&2
    exit 2
fi
if [[ ! "$environment_name" =~ ^[A-Za-z0-9._-]+$ ]] ||
   [[ "$environment_name" == "base" || "$environment_name" == "root" ]]; then
    printf 'ERROR: invalid or protected Conda environment name: %s\n' \
        "$environment_name" >&2
    exit 2
fi

host_os="$(uname -s)"
host_architecture="$(uname -m)"
case "$host_os:$host_architecture" in
    Darwin:arm64)
        if $all_tools && ! /usr/bin/arch -x86_64 /usr/bin/true; then
            printf 'ERROR: --all-tools requires Rosetta 2 for pinned Intel-only macOS tools. Install it with softwareupdate --install-rosetta and rerun.\n' >&2
            exit 1
        fi
        ;;
    Darwin:x86_64)
        if $all_tools; then
            printf 'ERROR: --all-tools supports Apple-silicon macOS and glibc Linux x86-64; found Intel macOS.\n' >&2
            exit 1
        fi
        ;;
    Linux:x86_64|Linux:amd64)
        if ldd --version 2>&1 | grep -qi musl; then
            printf 'ERROR: Alpine/musl Linux is not supported; use a glibc-based Linux host.\n' >&2
            exit 1
        fi
        ;;
    *)
        printf 'ERROR: unsupported native host %s/%s.\n' \
            "$host_os" "$host_architecture" >&2
        printf 'Supported hosts are macOS arm64/x86_64 and glibc Linux x86_64.\n' >&2
        exit 1
        ;;
esac

environment_manager=conda
if $all_tools; then
    if ! command -v mamba >/dev/null 2>&1; then
        printf 'ERROR: mamba was not found. Install and initialise Miniforge with Mamba, then rerun this command.\n' >&2
        exit 1
    fi
    if ! command -v conda >/dev/null 2>&1; then
        printf 'ERROR: conda was not found beside Mamba. Install Miniforge, then rerun this command.\n' >&2
        exit 1
    fi
    environment_manager=mamba
elif command -v mamba >/dev/null 2>&1; then
    environment_manager=mamba
elif ! command -v conda >/dev/null 2>&1; then
    printf 'ERROR: Conda was not found. Install and initialise Miniforge, then rerun this command.\n' >&2
    exit 1
fi
environment_manager_path="$(command -v "$environment_manager")"
environment_runner_path="$environment_manager_path"
if $all_tools; then
    environment_runner_path="$(command -v conda)"
fi

environment_exists=false
if "$environment_manager" run -n "$environment_name" true >/dev/null 2>&1; then
    environment_exists=true
fi

printf 'PostGWAS installation\n'
printf '  Host: %s/%s\n' "$host_os" "$host_architecture"
printf '  Environment: %s\n' "$environment_name"
printf '  Solver: %s\n' "$environment_manager"
printf '  Source-build jobs: %s\n' "$build_jobs"
if $all_tools; then
    printf '  Software coverage: complete native stack\n'
else
    printf '  Software coverage: portable stack\n'
fi

if $environment_exists; then
    if ! $update; then
        printf 'ERROR: Conda environment %s already exists.\n' "$environment_name" >&2
        printf 'Select a new environment with --name. Portable installs may instead use --update.\n' >&2
        exit 1
    fi
    "$environment_manager" env update \
        --yes \
        --name "$environment_name" \
        --file "$environment_file"
else
    "$environment_manager" env create \
        --yes \
        --name "$environment_name" \
        --file "$environment_file"
fi

if [[ "$host_os" == "Linux" ]]; then
    linux_packages=()
    while IFS= read -r package_name; do
        if [[ -n "$package_name" ]]; then
            linux_packages+=("$package_name")
        fi
    done < <(
        awk '$1 == "-" && $2 == "sel(linux):" { print $3 }' \
            "$environment_file"
    )
    if ((${#linux_packages[@]} == 0)); then
        printf 'ERROR: environment.yml contains no Linux-selected package definitions.\n' >&2
        exit 1
    fi
    "$environment_manager" install \
        --yes \
        --name "$environment_name" \
        --channel conda-forge \
        --channel bioconda \
        --override-channels \
        "${linux_packages[@]}"
fi

pip_arguments=(python -m pip install --no-deps --no-build-isolation)
if $editable; then
    pip_arguments+=(-e)
fi
pip_arguments+=("$repository_root")

"$environment_manager" run -n "$environment_name" "${pip_arguments[@]}"
POSTGWAS_INSTALL_JOBS="$build_jobs" \
    "$environment_manager" run -n "$environment_name" \
    bash "$script_directory/install_bcftools_liftover.sh"

if $all_tools; then
    # Use Conda only as the environment launcher here. Running this through
    # `mamba run` would keep a parent Mamba process alive while the child
    # installer creates its isolated Mamba runtimes, which can contend for the
    # shared package-cache lock on Linux.
    "$environment_runner_path" run --no-capture-output \
        -n "$environment_name" \
        bash "$script_directory/install_external_tools.sh" \
        --mamba "$environment_manager_path" \
        --jobs "$build_jobs"
fi

printf 'Verifying installation\n'
"$environment_manager" run -n "$environment_name" python -m pip check
"$environment_manager" run -n "$environment_name" postgwas --help >/dev/null
"$environment_manager" run -n "$environment_name" \
    bcftools plugin -l | grep -qx liftover

portable_commands=(python Rscript bcftools tabix bedtools pigz plink scdrs)
for executable in "${portable_commands[@]}"; do
    "$environment_manager" run -n "$environment_name" \
        bash -c 'command -v "$1" >/dev/null' bash "$executable"
done

if [[ "$host_os" == "Linux" ]]; then
    linux_commands=(plink2 gcta64 bgenix)
    for executable in "${linux_commands[@]}"; do
        "$environment_manager" run -n "$environment_name" \
            bash -c 'command -v "$1" >/dev/null' bash "$executable"
    done
fi

if $all_tools; then
    "$environment_manager" run -n "$environment_name" \
        bash "$script_directory/verify_mamba_installation.sh"
fi

printf '\nPostGWAS software installation completed successfully.\n'
printf 'Activate it with: conda activate %s\n' "$environment_name"
if $all_tools; then
    printf 'MAGMA, LDSC, FLaMEs, PoPS, FINEMAP, MiXeR, K-POPS, CALDERA, and supporting tools were verified.\n'
else
    printf 'This is the portable subset. Use --all-tools for the complete stack.\n'
fi
printf 'Scientific reference datasets, API credentials, and licensed cache data are not installed by this command.\n'
