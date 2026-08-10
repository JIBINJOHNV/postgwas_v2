


#!/bin/bash

# --- Paths ---
resource_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37"
plink_file="EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001"
output_folder="${resource_folder}/LD_ref_EUR"

mkdir -p ${output_folder}

# --- Process Chromosome-wise ---
for chrom in {9..12}; do
  echo "------------------------------------------"
  echo "Step 1: Running PLINK for Chromosome ${chrom}"
  echo "------------------------------------------"

  docker run --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
    -it jibinjv/postgwas:1.3 plink \
          --bfile ${resource_folder}/${plink_file} \
          --chr ${chrom} \
          --keep-allele-order \
          --r2 gz \
          --ld-window 99999 \
          --ld-window-r2 0.05 \
          --maf 0.00001 \
          --out ${output_folder}/EUR_chr${chrom}

        #           --ld-window-kb 1000 \

  echo "Step 2: Sorting and Indexing for Tabix..."
  # 1. Remove PLINK header (sed '1d')
  # 2. Sort by Position (Column 2) numerically
  # 3. Compress with bgzip
 echo "Step 2: Sorting and Indexing for Tabix..."

  # Standardizing the format specifically for Tabix parsing
  zcat ${output_folder}/EUR_chr${chrom}.ld.gz | sed '1d' | \
  awk '{$1=$1; print}' OFS='\t' | \
  sort -k2,2n | bgzip > ${output_folder}/EUR_chr${chrom}.sorted.ld.gz

  mv ${output_folder}/EUR_chr${chrom}.sorted.ld.gz ${output_folder}/EUR_chr${chrom}.ld.gz
  # Indexing: Col 1=Chr, Col 2=Pos, Col 2=End (Point query)
  tabix -f -s 1 -b 2 -e 2 ${output_folder}/EUR_chr${chrom}.ld.gz
  # Clean up the large unsorted file
  echo "Finished Chromosome ${chrom}."
done
