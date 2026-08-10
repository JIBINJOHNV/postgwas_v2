#!/usr/bin/env bash
set -euo pipefail

CALDERA_COMMIT="81a8a0308741ae986660f711bbf6a7abbd3bca19"

usage() {
    printf '%s\n' \
        "Usage: $0 --resource-root PATH --output-root PATH --python PATH" \
        "          [--magma-prefix PREFIX]" \
        "" \
        "Runs the bundled CALDERA example through PostGWAS. When a matching" \
        "MAGMA .genes.out/.genes.raw prefix is supplied, also runs K-POPS."
}

resource_root=""
output_root=""
magma_prefix=""
python_bin=""
while (($#)); do
    case "$1" in
        --resource-root)
            resource_root="${2:-}"
            shift 2
            ;;
        --output-root)
            output_root="${2:-}"
            shift 2
            ;;
        --magma-prefix)
            magma_prefix="${2:-}"
            shift 2
            ;;
        --python)
            python_bin="${2:-}"
            shift 2
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
if [[ -z "$resource_root" || -z "$output_root" || -z "$python_bin" ]]; then
    usage >&2
    exit 2
fi
if [[ ! -x "$python_bin" ]]; then
    printf 'Configured Python is not executable: %s\n' "$python_bin" >&2
    exit 1
fi

resource_root="$(cd "$resource_root" && pwd -P)"
mkdir -p "$output_root"
output_root="$(cd "$output_root" && pwd -P)"
caldera_repository="${resource_root}/caldera/${CALDERA_COMMIT}"
kernel_prefix="${resource_root}/kpops/kernels/GRCh37/pops_features_standardized_linear"
python_directory="$(cd "$(dirname "$python_bin")" && pwd -P)"
export PATH="${python_directory}:${PATH}"

"$python_bin" -m postgwas caldera \
    --pops-file "${caldera_repository}/example/inputs/lange2025_parkinsons_disease_pops.preds" \
    --credible-set-file "${caldera_repository}/example/inputs/lange2025_parkinsons_disease_cred_sets.tsv" \
    --caldera-genome-build GRCh37 \
    --dataset-id caldera_upstream_example \
    --output-directory "${output_root}/caldera" \
    --overwrite

if [[ -n "$magma_prefix" ]]; then
    "$python_bin" -m postgwas kpops \
        --magma-association-prefix "$magma_prefix" \
        --kpops-gene-annotation-file \
        "${resource_root}/pops/GRCh37_gene_annot_jun10.txt" \
        --kernel-matrix-prefix "$kernel_prefix" \
        --kpops-genome-build GRCh37 \
        --dataset-id kpops_resource_test \
        --output-directory "${output_root}/kpops" \
        --overwrite
else
    printf '%s\n' \
        "K-POPS execution skipped: provide --magma-prefix for matching" \
        "GRCh37 MAGMA .genes.out and .genes.raw files."
fi
