#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=tools/setup/software_versions.env
source "$script_directory/../setup/software_versions.env"

KPOPS_COMMIT="$POSTGWAS_KPOPS_COMMIT"
KPOPS_ARCHIVE_SHA256="$POSTGWAS_KPOPS_SHA256"
CALDERA_COMMIT="$POSTGWAS_CALDERA_COMMIT"
CALDERA_ARCHIVE_SHA256="$POSTGWAS_CALDERA_SHA256"
GENE_ANNOTATION_SHA256="940086e0781e9b244da76d27f02cc25dd45c957ee4415915cd1ced37b63d1e72"

usage() {
    printf '%s\n' \
        "Usage: $0 --resource-root PATH [--python PATH] [--software-only] [--prepare-linear-kernel]" \
        "" \
        "Downloads pinned K-POPS and CALDERA sources into versioned resource" \
        "directories. Supplying --python installs both tools into that environment." \
        "--software-only skips validation of separately supplied PoPS feature data." \
        "With --prepare-linear-kernel, it also derives the standardized GRCh37" \
        "linear kernel from an existing PoPS munged feature matrix."
}

resource_root=""
python_bin=""
prepare_kernel=false
software_only=false
while (($#)); do
    case "$1" in
        --resource-root)
            resource_root="${2:-}"
            shift 2
            ;;
        --python)
            python_bin="${2:-}"
            shift 2
            ;;
        --prepare-linear-kernel)
            prepare_kernel=true
            shift
            ;;
        --software-only)
            software_only=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            printf 'Unknown argument: %s\n' "$1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$resource_root" ]]; then
    printf '%s\n' "--resource-root is required" >&2
    exit 2
fi
if [[ "$prepare_kernel" == true && -z "$python_bin" ]]; then
    printf '%s\n' "--python is required with --prepare-linear-kernel" >&2
    exit 2
fi
if [[ "$prepare_kernel" == true && "$software_only" == true ]]; then
    printf '%s\n' "--prepare-linear-kernel cannot be combined with --software-only" >&2
    exit 2
fi
for executable in curl tar shasum mktemp; do
    command -v "$executable" >/dev/null || {
        printf 'Required executable not found: %s\n' "$executable" >&2
        exit 1
    }
done
if [[ -n "$python_bin" && ! -x "$python_bin" ]]; then
    printf 'Configured Python is not executable: %s\n' "$python_bin" >&2
    exit 1
fi
if [[ -n "$python_bin" ]]; then
    for executable in cmp cp install; do
        command -v "$executable" >/dev/null || {
            printf 'Required executable not found: %s\n' "$executable" >&2
            exit 1
        }
    done
fi

resource_root="$(cd "$resource_root" && pwd -P)"
work_directory="$(mktemp -d "${resource_root}/.kpops-caldera-preparation.XXXXXX")"
caldera_staging=""
cleanup() {
    if [[ -n "${work_directory:-}" && -d "$work_directory" ]]; then
        rm -rf -- "$work_directory"
    fi
    if [[ -n "${caldera_staging:-}" && -d "$caldera_staging" ]]; then
        rm -rf -- "$caldera_staging"
    fi
}
trap cleanup EXIT

install_archive() {
    local name="$1"
    local commit="$2"
    local expected_sha256="$3"
    local archive_url="$4"
    local archive_root="$5"
    local destination="$6"
    shift 6
    local marker="${destination}/POSTGWAS_UPSTREAM_COMMIT"
    local checksum_marker="${destination}/POSTGWAS_ARCHIVE_SHA256"
    local required

    if [[ -d "$destination" ]]; then
        if [[ ! -f "$marker" || "$(<"$marker")" != "$commit" ]] \
            || [[ ! -f "$checksum_marker" ]] \
            || [[ "$(<"$checksum_marker")" != "$expected_sha256" ]]; then
            printf '%s\n' \
                "Refusing to replace an unverified existing directory: $destination" >&2
            exit 1
        fi
        for required in "$@"; do
            [[ -s "${destination}/${required}" ]] || {
                printf 'Existing %s resource is incomplete: %s\n' \
                    "$name" "${destination}/${required}" >&2
                exit 1
            }
        done
        printf 'Validated existing %s source: %s\n' "$name" "$destination"
        return
    fi

    local archive="${work_directory}/${name}.tar.gz"
    local extracted="${work_directory}/${archive_root}"
    printf 'Downloading pinned %s source...\n' "$name"
    curl --fail --location --silent --show-error "$archive_url" -o "$archive"
    printf '%s  %s\n' "$expected_sha256" "$archive" | shasum -a 256 --check --status
    tar -xzf "$archive" -C "$work_directory"
    for required in "$@"; do
        [[ -s "${extracted}/${required}" ]] || {
            printf 'Downloaded %s archive is missing: %s\n' "$name" "$required" >&2
            exit 1
        }
    done
    mkdir -p "$(dirname "$destination")"
    printf '%s\n' "$commit" >"${extracted}/POSTGWAS_UPSTREAM_COMMIT"
    printf '%s\n' "$expected_sha256" >"${extracted}/POSTGWAS_ARCHIVE_SHA256"
    mv "$extracted" "$destination"
    printf 'Installed %s source: %s\n' "$name" "$destination"
}

kpops_destination="${resource_root}/kpops/software/${KPOPS_COMMIT}"
caldera_destination="${resource_root}/caldera/${CALDERA_COMMIT}"
install_archive \
    "kpops" "$KPOPS_COMMIT" "$KPOPS_ARCHIVE_SHA256" \
    "https://github.com/JasonTan-code/k-pops/archive/${KPOPS_COMMIT}.tar.gz" \
    "k-pops-${KPOPS_COMMIT}" "$kpops_destination" \
    "k-pops.py" "prepare_kernel.py" "data/gene_annot_jun10.txt"
install_archive \
    "caldera" "$CALDERA_COMMIT" "$CALDERA_ARCHIVE_SHA256" \
    "https://github.com/kheilbron/caldera/archive/${CALDERA_COMMIT}.tar.gz" \
    "caldera-${CALDERA_COMMIT}" "$caldera_destination" \
    "z_caldera.R" "data/high_h2_coding_SNPs.tsv" \
    "data/gene_locations_gencode_v47_grch37_all_protein_coding.tsv" \
    "data/gene_locations_gencode_v47_grch38_all_protein_coding.tsv" \
    "trained_models/caldera_model_no_covs.rds"

if [[ -n "$python_bin" ]]; then
    kpops_command="$(dirname "$python_bin")/k-pops.py"
    if [[ -e "$kpops_command" ]] && ! cmp -s \
        "${kpops_destination}/k-pops.py" "$kpops_command"; then
        printf 'Refusing to replace a different installed K-POPS command: %s\n' \
            "$kpops_command" >&2
        exit 1
    fi
    install -m 0755 "${kpops_destination}/k-pops.py" "$kpops_command"
    printf 'Installed K-POPS command: %s\n' "$kpops_command"

    python_environment="$(cd "$(dirname "$python_bin")/.." && pwd -P)"
    caldera_install="${python_environment}/share/postgwas/caldera"
    if [[ -e "$caldera_install" && ! -d "$caldera_install" ]]; then
        printf 'CALDERA installation path is not a directory: %s\n' \
            "$caldera_install" >&2
        exit 1
    fi
    if [[ -d "$caldera_install" ]]; then
        if [[ ! -f "${caldera_install}/POSTGWAS_UPSTREAM_COMMIT" ]] \
            || [[ "$(<"${caldera_install}/POSTGWAS_UPSTREAM_COMMIT")" != "$CALDERA_COMMIT" ]] \
            || [[ ! -f "${caldera_install}/POSTGWAS_ARCHIVE_SHA256" ]] \
            || [[ "$(<"${caldera_install}/POSTGWAS_ARCHIVE_SHA256")" != "$CALDERA_ARCHIVE_SHA256" ]]; then
            printf 'Refusing to replace a different CALDERA installation: %s\n' \
                "$caldera_install" >&2
            exit 1
        fi
    else
        mkdir -p "$(dirname "$caldera_install")"
        caldera_staging="$(mktemp -d "${caldera_install}.XXXXXX")"
        cp -R "${caldera_destination}/." "$caldera_staging"
        mv "$caldera_staging" "$caldera_install"
        caldera_staging=""
    fi
    for required in z_caldera.R \
        data/high_h2_coding_SNPs.tsv \
        data/gene_locations_gencode_v47_grch37_all_protein_coding.tsv \
        data/gene_locations_gencode_v47_grch38_all_protein_coding.tsv \
        trained_models/caldera_model_no_covs.rds; do
        [[ -s "${caldera_install}/${required}" ]] || {
            printf 'Installed CALDERA resource is incomplete: %s\n' \
                "${caldera_install}/${required}" >&2
            exit 1
        }
    done
    printf 'Installed CALDERA repository: %s\n' "$caldera_install"
fi

gene_annotation="${resource_root}/pops/GRCh37_gene_annot_jun10.txt"
feature_prefix="${resource_root}/pops/features_munged/pops_features"
if [[ "$software_only" == false ]]; then
    printf '%s  %s\n' "$GENE_ANNOTATION_SHA256" "$gene_annotation" \
        | shasum -a 256 --check --status || {
            printf 'The existing GRCh37 PoPS annotation is missing or differs from upstream: %s\n' \
                "$gene_annotation" >&2
            exit 1
        }
fi

if [[ "$prepare_kernel" == true ]]; then
    "$python_bin" - "$gene_annotation" "$feature_prefix" <<'PY'
from pathlib import Path
import glob
import sys

import numpy as np
import pandas as pd

annotation_path, prefix = sys.argv[1:]
annotation = pd.read_csv(annotation_path, sep=r"\s+", usecols=["ENSGID"])
rows_path = Path(prefix + ".rows.txt")
rows = pd.read_csv(rows_path, header=None, names=["ENSGID"], dtype=str)
matrices = sorted(glob.glob(prefix + ".mat.*.npy"))
columns = sorted(glob.glob(prefix + ".cols.*.txt"))
if not matrices or len(matrices) != len(columns):
    raise SystemExit("PoPS munged feature chunks are missing or unpaired")
if rows.empty or rows["ENSGID"].duplicated().any():
    raise SystemExit("PoPS feature row identifiers must be non-empty and unique")
if not set(rows["ENSGID"]).issubset(set(annotation["ENSGID"].astype(str))):
    raise SystemExit("PoPS feature genes are not a subset of the GRCh37 annotation")
row_count = len(rows)
for path in matrices:
    matrix = np.load(path, mmap_mode="r")
    if matrix.ndim != 2 or matrix.shape[0] != row_count:
        raise SystemExit(f"Invalid PoPS feature chunk shape: {path} {matrix.shape}")
print(f"Validated {row_count} genes across {len(matrices)} PoPS feature chunks")
PY

    kernel_directory="${resource_root}/kpops/kernels/GRCh37"
    kernel_prefix="${kernel_directory}/pops_features_standardized_linear"
    kernel_file="${kernel_prefix}.bin"
    genes_file="${kernel_prefix}.genes"
    if [[ -e "$kernel_file" || -e "$genes_file" ]]; then
        [[ -s "$kernel_file" && -s "$genes_file" ]] || {
            printf 'Incomplete existing K-POPS kernel output: %s\n' "$kernel_prefix" >&2
            exit 1
        }
        printf 'Retaining existing K-POPS kernel: %s\n' "$kernel_prefix"
    else
        kernel_staging="${work_directory}/kernel"
        mkdir -p "$kernel_staging"
        "$python_bin" "${kpops_destination}/prepare_kernel.py" \
            --input_prefix "$feature_prefix" \
            --kernel_type linear \
            --standardize \
            --output_prefix "${kernel_staging}/pops_features_standardized_linear"
        mkdir -p "$kernel_directory"
        mv "${kernel_staging}/pops_features_standardized_linear.bin" "$kernel_file"
        mv "${kernel_staging}/pops_features_standardized_linear.genes" "$genes_file"
    fi
    "$python_bin" - "$kernel_file" "$genes_file" <<'PY'
from pathlib import Path
import sys

kernel_path, genes_path = map(Path, sys.argv[1:])
genes = [line.strip() for line in genes_path.read_text().splitlines() if line.strip()]
expected = len(genes) * len(genes) * 4
observed = kernel_path.stat().st_size
if len(genes) != len(set(genes)):
    raise SystemExit("K-POPS kernel genes are not unique")
if observed != expected:
    raise SystemExit(f"K-POPS kernel size is {observed}; expected {expected}")
print(f"Validated K-POPS kernel dimensions for {len(genes)} genes")
PY
    kernel_sha256="$(shasum -a 256 "$kernel_file" | awk '{print $1}')"
    genes_sha256="$(shasum -a 256 "$genes_file" | awk '{print $1}')"
    kernel_manifest="${kernel_prefix}.manifest.tsv"
    printf '%s\t%s\n' \
        "field" "value" \
        "genome_build" "GRCh37" \
        "kernel_type" "linear" \
        "feature_standardization" "standardized" \
        "kpops_commit" "$KPOPS_COMMIT" \
        "gene_annotation" "$gene_annotation" \
        "gene_annotation_sha256" "$GENE_ANNOTATION_SHA256" \
        "feature_prefix" "$feature_prefix" \
        "kernel_sha256" "$kernel_sha256" \
        "genes_sha256" "$genes_sha256" \
        >"$kernel_manifest"
fi

printf '%s\n' \
    "K-POPS script: ${kpops_destination}/k-pops.py" \
    "CALDERA repository: ${caldera_destination}"
if [[ "$software_only" == false ]]; then
    printf '%s\n' \
        "K-POPS annotation: ${gene_annotation}" \
        "K-POPS feature prefix: ${feature_prefix}"
fi
if [[ "$prepare_kernel" == true ]]; then
    printf 'K-POPS kernel prefix: %s\n' "$kernel_prefix"
fi
