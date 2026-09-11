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
    "  --gzip COMMAND          gzip-compatible decompressor, such as pigz." \
    "  --bgzip COMMAND         bgzip executable." \
    "  --tabix COMMAND         tabix executable." \
    "  --window-kb N           Largest stored LD distance in kb." \
    "  --ld-window-variants N  Largest number of variants in one PLINK LD window." \
    "  --minimum-r2 VALUE      Smallest stored r-squared." \
    "  --minimum-maf VALUE     PLINK MAF filter." \
    "  --threads N             Total compute thread budget." \
    "  --chromosome-workers N  Concurrent chromosome jobs within the thread budget." \
    "  --orientation VALUE     symmetric_first_endpoint or upper_triangle_dual_index." \
    "" \
    "Optional:" \
    "  --reuse-existing-ld     Keep existing indexed forward LD files and build only inventory/reverse resources."
}

bfile=""
output_dir=""
population=""
genome_build=""
chromosomes=""
plink_bin=""
gzip_bin=""
bgzip_bin=""
tabix_bin=""
window_kb=""
ld_window_variants=""
minimum_r2=""
minimum_maf=""
threads=""
chromosome_workers=""
orientation=""
reuse_existing_ld=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) bfile="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --population) population="$2"; shift 2 ;;
    --genome-build) genome_build="$2"; shift 2 ;;
    --chromosomes) chromosomes="$2"; shift 2 ;;
    --plink) plink_bin="$2"; shift 2 ;;
    --gzip) gzip_bin="$2"; shift 2 ;;
    --bgzip) bgzip_bin="$2"; shift 2 ;;
    --tabix) tabix_bin="$2"; shift 2 ;;
    --window-kb) window_kb="$2"; shift 2 ;;
    --ld-window-variants) ld_window_variants="$2"; shift 2 ;;
    --minimum-r2) minimum_r2="$2"; shift 2 ;;
    --minimum-maf) minimum_maf="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --chromosome-workers) chromosome_workers="$2"; shift 2 ;;
    --orientation) orientation="$2"; shift 2 ;;
    --reuse-existing-ld) reuse_existing_ld=true; shift ;;
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
[[ -n "$gzip_bin" ]] || missing+=("--gzip")
[[ -n "$bgzip_bin" ]] || missing+=("--bgzip")
[[ -n "$tabix_bin" ]] || missing+=("--tabix")
[[ -n "$window_kb" ]] || missing+=("--window-kb")
[[ -n "$ld_window_variants" ]] || missing+=("--ld-window-variants")
[[ -n "$minimum_r2" ]] || missing+=("--minimum-r2")
[[ -n "$minimum_maf" ]] || missing+=("--minimum-maf")
[[ -n "$threads" ]] || missing+=("--threads")
[[ -n "$chromosome_workers" ]] || missing+=("--chromosome-workers")
[[ -n "$orientation" ]] || missing+=("--orientation")
if [[ ${#missing[@]} -gt 0 ]]; then
  printf 'Required argument not provided: %s\n' "${missing[*]}" >&2
  usage >&2
  exit 2
fi
if [[ "$genome_build" != "GRCh37" && "$genome_build" != "GRCh38" ]]; then
  printf 'Invalid --genome-build: %s\n' "$genome_build" >&2
  exit 2
fi
if [[ "$orientation" != "symmetric_first_endpoint" && "$orientation" != "upper_triangle_dual_index" ]]; then
  printf 'Invalid --orientation: %s\n' "$orientation" >&2
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
if [[ ! "$chromosome_workers" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Invalid --chromosome-workers: %s; expected a positive integer.\n' "$chromosome_workers" >&2
  exit 2
fi
if (( chromosome_workers > threads )); then
  printf 'Invalid compute allocation: --chromosome-workers (%s) cannot exceed --threads (%s).\n' \
    "$chromosome_workers" "$threads" >&2
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
for executable in "$plink_bin" "$gzip_bin" "$bgzip_bin" "$tabix_bin"; do
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
if [[ "$reuse_existing_ld" == true ]]; then
  for chromosome in $chromosome_values; do
    reusable_prefix="${output_dir}/${population}_chr${chromosome}"
    if [[ ! -s "${reusable_prefix}.ld.gz" || ! -s "${reusable_prefix}.ld.gz.tbi" ]]; then
      printf 'Missing reusable LD file or index: %s\n' "${reusable_prefix}.ld.gz" >&2
      exit 2
    fi
  done
fi

# Variant membership must come from the same allele-preserving PLINK fileset as
# the LD pairs. Calculate frequencies once for the complete requested fileset:
# running --freq separately for every chromosome repeatedly scans the same BIM
# and genotype data and dominates migration time on full 1000 Genomes panels.
inventory_working="${output_dir}/.${population}.inventory_working"
inventory_chromosomes="${chromosome_values//$'\n'/,}"
inventory_chromosomes="${inventory_chromosomes// /,}"
worker_thread_budget=$(( threads / chromosome_workers ))
# Each active reversal also runs one decompressor, awk, and sort process. Reserve
# those three single-threaded slots before assigning threads to BGZF compression.
compression_threads=$(( worker_thread_budget - 3 ))
if (( compression_threads < 1 )); then
  compression_threads=1
fi
active_pids=()
active_chromosomes=()
for chromosome in $chromosome_values; do
  rm -f "${output_dir}/${population}_chr${chromosome}.variants.unsorted.tsv"
done
cleanup_inventory_working() {
  for active_pid in "${active_pids[@]:-}"; do
    if [[ -n "$active_pid" ]]; then
      kill "$active_pid" 2>/dev/null || true
    fi
  done
  for active_pid in "${active_pids[@]:-}"; do
    if [[ -n "$active_pid" ]]; then
      wait "$active_pid" 2>/dev/null || true
    fi
  done
  rm -f "${inventory_working}.frq" \
    "${inventory_working}.log" \
    "${inventory_working}.nosex"
  for cleanup_chromosome in $chromosome_values; do
    rm -f "${output_dir}/${population}_chr${cleanup_chromosome}.variants.unsorted.tsv"
  done
  if [[ -n "${manifest_working:-}" ]]; then
    rm -f "$manifest_working"
  fi
}
trap cleanup_inventory_working EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

printf 'Preparing allele-aware variant inventory with %s\n' "$plink_version"
"$plink_bin" \
  --bfile "$bfile" \
  --keep-allele-order \
  --maf "$minimum_maf" \
  --freq \
  --threads "$threads" \
  --out "$inventory_working"
awk \
  -v allowed_csv="$inventory_chromosomes" \
  -v output_dir="$output_dir" \
  -v population="$population" '
    BEGIN {
      count=split(allowed_csv, requested, ",")
      for (item=1; item <= count; item++) allowed[requested[item]]=1
    }
    NR == FNR { chrom[$2]=$1; position[$2]=$4; next }
    FNR > 1 && ($2 in position) && (chrom[$2] in allowed) {
      if ($5 !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ ||
          ($5 + 0) < 0 || ($5 + 0) > 1) {
        print "Invalid PLINK allele frequency for variant " $2 ": " $5 > "/dev/stderr"
        invalid_frequency=1
        next
      }
      allele_1=toupper($3); allele_2=toupper($4)
      if (allele_1 <= allele_2) { low=allele_1; high=allele_2 }
      else { low=allele_2; high=allele_1 }
      canonical=chrom[$2] "_" position[$2] "_" low "_" high
      destination=output_dir "/" population "_chr" chrom[$2] ".variants.unsorted.tsv"
      allele_frequency=$5 + 0
      minor_allele_frequency=(allele_frequency <= 0.5) ? $5 : 1 - allele_frequency
      print chrom[$2], position[$2], $2, canonical, allele_1, allele_2, minor_allele_frequency > destination
      destinations[destination]=1
    }
    END {
      for (destination in destinations) close(destination)
      if (invalid_frequency) exit 2
    }
  ' OFS='\t' "${bfile}.bim" "${inventory_working}.frq"

process_chromosome() {
  trap - EXIT INT TERM
  local chromosome="$1"
  local prefix="${output_dir}/${population}_chr${chromosome}"
  local inventory_unsorted
  local reverse
  local forward
  local symmetric
  printf 'Preparing chromosome %s reference sidecars\n' "$chromosome"
  if [[ "$reuse_existing_ld" == false ]]; then
    "$plink_bin" \
      --bfile "$bfile" \
      --chr "$chromosome" \
      --keep-allele-order \
      --r2 gz \
      --ld-window "$ld_window_variants" \
      --ld-window-kb "$window_kb" \
      --ld-window-r2 "$minimum_r2" \
      --maf "$minimum_maf" \
      --threads "$worker_thread_budget" \
      --out "$prefix"
  fi

  inventory_unsorted="${prefix}.variants.unsorted.tsv"
  if [[ ! -s "$inventory_unsorted" ]]; then
    printf 'No variants passed the inventory contract for chromosome %s.\n' "$chromosome" >&2
    exit 2
  fi
  LC_ALL=C sort -k1,1V -k2,2n -k3,3 "$inventory_unsorted" |
    "$bgzip_bin" -@ "$compression_threads" -c > "${prefix}.variants.tsv.gz"
  "$tabix_bin" -f -s 1 -b 2 -e 2 "${prefix}.variants.tsv.gz"
  rm -f "$inventory_unsorted"

  if [[ "$orientation" == "upper_triangle_dual_index" ]]; then
    reverse="${prefix}.reverse.ld.gz"
    "$gzip_bin" -cd "${prefix}.ld.gz" |
      awk '$2 ~ /^[0-9]+$/ {
        print $4, $5, $6, $1, $2, $3, $7
      }' OFS='\t' |
      LC_ALL=C sort -k1,1V -k2,2n -k3,3 -k5,5n |
      "$bgzip_bin" -@ "$compression_threads" -c > "$reverse"
    "$tabix_bin" -f -s 1 -b 2 -e 2 "$reverse"
    if [[ "$reuse_existing_ld" == false ]]; then
      forward="${prefix}.forward.ld.gz"
      "$gzip_bin" -cd "${prefix}.ld.gz" |
        awk '$2 ~ /^[0-9]+$/ { print $1, $2, $3, $4, $5, $6, $7 }' OFS='\t' |
        "$bgzip_bin" -@ "$compression_threads" -c > "$forward"
      mv "$forward" "${prefix}.ld.gz"
      "$tabix_bin" -f -s 1 -b 2 -e 2 "${prefix}.ld.gz"
    fi
  else
    symmetric="${prefix}.symmetric.ld.gz"
    "$gzip_bin" -cd "${prefix}.ld.gz" |
      awk '$2 ~ /^[0-9]+$/ {
        print $1, $2, $3, $4, $5, $6, $7
        if (!($1 == $4 && $2 == $5 && $3 == $6)) {
          print $4, $5, $6, $1, $2, $3, $7
        }
      }' OFS='\t' |
      LC_ALL=C sort -u -k1,1V -k2,2n -k3,3 -k5,5n |
      "$bgzip_bin" -@ "$compression_threads" -c > "$symmetric"
    mv "$symmetric" "${prefix}.ld.gz"
    "$tabix_bin" -f -s 1 -b 2 -e 2 "${prefix}.ld.gz"
  fi
  printf 'Prepared chromosome %s reference sidecars\n' "$chromosome"
}

stop_active_jobs() {
  local pid
  for pid in "${active_pids[@]}"; do
    kill "$pid" 2>/dev/null || true
  done
  for pid in "${active_pids[@]}"; do
    wait "$pid" 2>/dev/null || true
  done
  active_pids=()
  active_chromosomes=()
}

wait_for_oldest_job() {
  local pid="${active_pids[0]}"
  local chromosome="${active_chromosomes[0]}"
  local status=0
  if wait "$pid"; then
    status=0
  else
    status=$?
  fi
  active_pids=("${active_pids[@]:1}")
  active_chromosomes=("${active_chromosomes[@]:1}")
  if (( status != 0 )); then
    printf 'Chromosome %s sidecar preparation failed with exit code %s.\n' \
      "$chromosome" "$status" >&2
    stop_active_jobs
    return "$status"
  fi
}

printf 'Chromosome workers: %s; thread budget per worker: %s; compression threads per worker: %s\n' \
  "$chromosome_workers" "$worker_thread_budget" "$compression_threads"
for chromosome in $chromosome_values; do
  (process_chromosome "$chromosome") &
  active_pids+=("$!")
  active_chromosomes+=("$chromosome")
  if (( ${#active_pids[@]} >= chromosome_workers )); then
    wait_for_oldest_job
  fi
done
while (( ${#active_pids[@]} > 0 )); do
  wait_for_oldest_job
done

manifest="${output_dir}/ld_reference.yaml"
manifest_working="${manifest}.working.$$"
{
  printf 'format_version: 2\n'
  printf 'genome_build: %s\n' "$genome_build"
  printf 'populations: [%s]\n' "$population"
  printf 'orientation: %s\n' "$orientation"
  printf 'file_pattern: "{population}_chr{chromosome}.ld.gz"\n'
  printf 'reverse_file_pattern: "{population}_chr{chromosome}.reverse.ld.gz"\n'
  printf 'variant_inventory_pattern: "{population}_chr{chromosome}.variants.tsv.gz"\n'
  printf 'columns: [chromosome_a, position_a, variant_a, chromosome_b, position_b, variant_b, r2]\n'
  printf 'variant_inventory_columns: [chromosome, position, reference_id, canonical_id, allele_1, allele_2, minor_allele_frequency]\n'
  printf 'window_kb: %s\n' "$window_kb"
  printf 'minimum_r2: %s\n' "$minimum_r2"
  printf 'minimum_maf: %s\n' "$minimum_maf"
  printf 'allele_order_preserved: true\n'
  printf 'plink_version: "%s"\n' "${plink_version//\"/\\\"}"
} > "$manifest_working"
mv "$manifest_working" "$manifest"

printf 'Prepared %s LD reference: %s\n' "$orientation" "$output_dir"
printf 'Manifest: %s\n' "$manifest"
