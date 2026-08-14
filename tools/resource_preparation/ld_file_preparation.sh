#!/usr/bin/env bash
set -euo pipefail

usage() {
  command_name="$(basename "$0")"
  printf '%s\n' \
    "Usage:" \
    "  ${command_name} --bfile PREFIX --output-dir DIR --population CODE --genome-build BUILD [required arguments]" \
    "" \
    "Required:" \
    "  --bfile PREFIX          PLINK 1 binary-fileset prefix (.bed/.bim/.fam)." \
    "  --output-dir DIR        Destination for LD files, indexes, and manifest." \
    "  --population CODE       Population label, for example EUR." \
    "  --genome-build BUILD    Coordinate build, GRCh37 or GRCh38." \
    "  --chromosomes LIST      Comma-separated chromosomes or 1-22." \
    "  --plink COMMAND         PLINK 1.9 executable." \
    "  --bgzip COMMAND         bgzip executable." \
    "  --tabix COMMAND         tabix executable." \
    "  --window-kb N           Largest stored LD distance in kb." \
    "  --ld-window-variants N  Largest number of variants in one PLINK LD window." \
    "  --minimum-r2 VALUE      Smallest stored r-squared." \
    "  --minimum-maf VALUE     PLINK MAF filter." \
    "  --threads N             PLINK thread budget."
}

bfile=""
output_dir=""
population=""
genome_build=""
chromosomes=""
plink_bin=""
bgzip_bin=""
tabix_bin=""
window_kb=""
ld_window_variants=""
minimum_r2=""
minimum_maf=""
threads=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) bfile="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --population) population="$2"; shift 2 ;;
    --genome-build) genome_build="$2"; shift 2 ;;
    --chromosomes) chromosomes="$2"; shift 2 ;;
    --plink) plink_bin="$2"; shift 2 ;;
    --bgzip) bgzip_bin="$2"; shift 2 ;;
    --tabix) tabix_bin="$2"; shift 2 ;;
    --window-kb) window_kb="$2"; shift 2 ;;
    --ld-window-variants) ld_window_variants="$2"; shift 2 ;;
    --minimum-r2) minimum_r2="$2"; shift 2 ;;
    --minimum-maf) minimum_maf="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

missing=()
[[ -n "$bfile" ]] || missing+=("--bfile")
[[ -n "$output_dir" ]] || missing+=("--output-dir")
[[ -n "$population" ]] || missing+=("--population")
[[ -n "$genome_build" ]] || missing+=("--genome-build")
[[ -n "$chromosomes" ]] || missing+=("--chromosomes")
[[ -n "$plink_bin" ]] || missing+=("--plink")
[[ -n "$bgzip_bin" ]] || missing+=("--bgzip")
[[ -n "$tabix_bin" ]] || missing+=("--tabix")
[[ -n "$window_kb" ]] || missing+=("--window-kb")
[[ -n "$ld_window_variants" ]] || missing+=("--ld-window-variants")
[[ -n "$minimum_r2" ]] || missing+=("--minimum-r2")
[[ -n "$minimum_maf" ]] || missing+=("--minimum-maf")
[[ -n "$threads" ]] || missing+=("--threads")
if [[ ${#missing[@]} -gt 0 ]]; then
  printf 'Required argument not provided: %s\n' "${missing[*]}" >&2
  usage >&2
  exit 2
fi
if [[ "$genome_build" != "GRCh37" && "$genome_build" != "GRCh38" ]]; then
  printf 'Invalid --genome-build: %s\n' "$genome_build" >&2
  exit 2
fi
if [[ ! "$population" =~ ^[A-Za-z][A-Za-z0-9_-]*$ ]]; then
  printf 'Invalid --population: %s\n' "$population" >&2
  exit 2
fi
if [[ "$chromosomes" != "1-22" && ! "$chromosomes" =~ ^([0-9]+|X|Y|XY|MT)(,([0-9]+|X|Y|XY|MT))*$ ]]; then
  printf 'Invalid --chromosomes: %s\n' "$chromosomes" >&2
  exit 2
fi
if [[ ! "$window_kb" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Invalid --window-kb: %s; expected a positive integer.\n' "$window_kb" >&2
  exit 2
fi
if [[ ! "$ld_window_variants" =~ ^[1-9][0-9]*$ ]] || (( ld_window_variants < 2 )); then
  printf 'Invalid --ld-window-variants: %s; expected an integer of at least 2.\n' "$ld_window_variants" >&2
  exit 2
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Invalid --threads: %s; expected a positive integer.\n' "$threads" >&2
  exit 2
fi
numeric_pattern='^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$'
if [[ ! "$minimum_r2" =~ $numeric_pattern ]] || ! awk -v value="$minimum_r2" 'BEGIN { exit !(value >= 0 && value <= 1) }'; then
  printf 'Invalid --minimum-r2: %s; expected a value in [0, 1].\n' "$minimum_r2" >&2
  exit 2
fi
if [[ ! "$minimum_maf" =~ $numeric_pattern ]] || ! awk -v value="$minimum_maf" 'BEGIN { exit !(value >= 0 && value <= 0.5) }'; then
  printf 'Invalid --minimum-maf: %s; expected a value in [0, 0.5].\n' "$minimum_maf" >&2
  exit 2
fi
for suffix in bed bim fam; do
  if [[ ! -s "${bfile}.${suffix}" ]]; then
    printf 'Missing or empty PLINK input: %s\n' "${bfile}.${suffix}" >&2
    exit 2
  fi
done
for executable in "$plink_bin" "$bgzip_bin" "$tabix_bin"; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    printf 'Executable not found: %s\n' "$executable" >&2
    exit 2
  fi
done

plink_version="$($plink_bin --version 2>&1)"
plink_version="${plink_version%%$'\n'*}"
if [[ ! "$plink_version" =~ PLINK[[:space:]]+v1[.]9 ]]; then
  printf 'Invalid --plink executable: expected PLINK 1.9, observed %s\n' "$plink_version" >&2
  exit 2
fi

mkdir -p "$output_dir"
chromosome_values="${chromosomes//,/ }"
if [[ "$chromosomes" == "1-22" ]]; then
  chromosome_values="$(seq 1 22)"
fi

for chromosome in $chromosome_values; do
  prefix="${output_dir}/${population}_chr${chromosome}"
  "$plink_bin" \
    --bfile "$bfile" \
    --chr "$chromosome" \
    --keep-allele-order \
    --r2 gz \
    --ld-window "$ld_window_variants" \
    --ld-window-kb "$window_kb" \
    --ld-window-r2 "$minimum_r2" \
    --maf "$minimum_maf" \
    --threads "$threads" \
    --out "$prefix"

  symmetric="${prefix}.symmetric.ld.gz"
  gzip -cd "${prefix}.ld.gz" |
    awk 'NR > 1 {
      $1=$1
      print $1, $2, $3, $4, $5, $6, $7
      if (!($1 == $4 && $2 == $5 && $3 == $6)) {
        print $4, $5, $6, $1, $2, $3, $7
      }
    }' OFS='\t' |
    LC_ALL=C sort -k1,1V -k2,2n -k3,3 -k5,5n |
    "$bgzip_bin" -c > "$symmetric"
  mv "$symmetric" "${prefix}.ld.gz"
  "$tabix_bin" -f -s 1 -b 2 -e 2 "${prefix}.ld.gz"
done

manifest="${output_dir}/ld_reference.yaml"
{
  printf 'format_version: 1\n'
  printf 'genome_build: %s\n' "$genome_build"
  printf 'populations: [%s]\n' "$population"
  printf 'orientation: symmetric_first_endpoint\n'
  printf 'file_pattern: "{population}_chr{chromosome}.ld.gz"\n'
  printf 'columns: [chromosome_a, position_a, variant_a, chromosome_b, position_b, variant_b, r2]\n'
  printf 'window_kb: %s\n' "$window_kb"
  printf 'minimum_r2: %s\n' "$minimum_r2"
  printf 'allele_order_preserved: true\n'
  printf 'plink_version: "%s"\n' "${plink_version//\"/\\\"}"
} > "$manifest"

printf 'Prepared symmetric LD reference: %s\n' "$output_dir"
printf 'Manifest: %s\n' "$manifest"
