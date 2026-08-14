#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository_root="$(cd "$script_directory/../.." && pwd)"
environment_file="$repository_root/environment.yml"
editable=false
update=false

usage() {
    cat <<'EOF'
Usage: bash tools/setup/install_postgwas.sh [options]

Create one PostGWAS Conda environment, install the current checkout, build the
bcftools liftover plugin inside that environment, and run installation checks.

Options:
  --name NAME   Override the environment.yml environment name.
  --editable    Install the checkout in editable development mode.
  --update      Update an existing environment instead of refusing to change it.
  -h, --help    Show this help message.

Supported native hosts:
  macOS arm64, macOS x86_64, and glibc-based Linux x86_64.
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

if [[ ! "$environment_name" =~ ^[A-Za-z0-9._-]+$ ]] ||
   [[ "$environment_name" == "base" || "$environment_name" == "root" ]]; then
    printf 'ERROR: invalid or protected Conda environment name: %s\n' \
        "$environment_name" >&2
    exit 2
fi

host_os="$(uname -s)"
host_architecture="$(uname -m)"
case "$host_os:$host_architecture" in
    Darwin:arm64|Darwin:x86_64)
        ;;
    Linux:x86_64|Linux:amd64)
        if ldd --version 2>&1 | grep -qi musl; then
            printf 'ERROR: Alpine/musl Linux is not supported; use a glibc-based Linux host or container.\n' >&2
            exit 1
        fi
        ;;
    *)
        printf 'ERROR: unsupported native host %s/%s.\n' \
            "$host_os" "$host_architecture" >&2
        printf 'Supported hosts are macOS arm64/x86_64 and glibc-based Linux x86_64.\n' >&2
        exit 1
        ;;
esac

if ! command -v conda >/dev/null 2>&1; then
    printf 'ERROR: conda was not found. Install and initialise Miniforge, then rerun this command.\n' >&2
    exit 1
fi

environment_manager=conda
if command -v mamba >/dev/null 2>&1; then
    environment_manager=mamba
fi

environment_exists=false
if conda run -n "$environment_name" true >/dev/null 2>&1; then
    environment_exists=true
fi

printf 'PostGWAS installation\n'
printf '  Host: %s/%s\n' "$host_os" "$host_architecture"
printf '  Environment: %s\n' "$environment_name"
printf '  Solver: %s\n' "$environment_manager"

if $environment_exists; then
    if ! $update; then
        printf 'ERROR: Conda environment %s already exists.\n' "$environment_name" >&2
        printf 'Rerun with --update, or select a new environment with --name.\n' >&2
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

pip_arguments=(python -m pip install --no-deps --no-build-isolation)
if $editable; then
    pip_arguments+=(-e)
fi
pip_arguments+=("$repository_root")

conda run -n "$environment_name" "${pip_arguments[@]}"
conda run -n "$environment_name" \
    bash "$script_directory/install_bcftools_liftover.sh"

printf 'Verifying installation\n'
conda run -n "$environment_name" python -m pip check
conda run -n "$environment_name" postgwas --help >/dev/null
conda run -n "$environment_name" bcftools plugin -l | grep -qx liftover

printf '\nPostGWAS software installation completed successfully.\n'
printf 'Activate it with: conda activate %s\n' "$environment_name"
printf 'Scientific reference datasets and separately licensed tools are not installed by this command.\n'
