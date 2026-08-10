


## pops

https://github.com/FinucaneLab/pops

## feature file downloaded from here
https://www.dropbox.com/scl/fo/ne7xhxkt4dwhvd52a59ub/AFKkJu7ACaun1uuE99kmTkc/data/PoPS.features.txt.gz?rlkey=ltdbcld1enyr1zefg1lfqm61i&e=1&dl=0

##  gene_annot_jun10.txt downloaded from https://github.com/FinucaneLab/pops/tree/master/example/data/utils

## converted



## 1000 genome
for chr in {1..23} ; do
    os.system(f"wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz &")
    os.system(f"wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi &")

done

wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi

rename.txt
INFO/EAS_AF INFO/EAS
INFO/AMR_AF INFO/AMR
INFO/AFR_AF INFO/AFR
INFO/EUR_AF INFO/EUR
INFO/SAS_AF INFO/SAS

for chr in {3..22}; do
    # Step 1: clean and split
    bcftools view \
        --threads 8 \
        -G \
        -v snps \
        -e 'REF="." || ALT="."' \
        ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        -Ou | \
    bcftools annotate --threads 8 -x INFO/SVLEN -Ou | \
    bcftools norm -m -any --threads 8 \
        -Oz -o tmp_chr${chr}.vcf.gz
    bcftools index tmp_chr${chr}.vcf.gz
    # Step 2: add rsID and clean INFO
    bcftools annotate \
        --threads 8 \
        -a ~/Documents/software_resources/resourses/postgwas/gwas2vcf/GRCh37/dbSNP/vcf_files/GRCh37_dbSNP157_chr${chr}.vcf.gz \
        -c ID \
        tmp_chr${chr}.vcf.gz \
        -Ou | \
    bcftools annotate \
        --threads 8 \
        -x ^INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/EUR_AF,INFO/SAS_AF \
        --rename-annots rename.txt \
        -Oz -o GRCh37_1000G_freq_chr${chr}.vcf.gz
    bcftools index GRCh37_1000G_freq_chr${chr}.vcf.gz
    rm tmp_chr${chr}.vcf.gz tmp_chr${chr}.vcf.gz.csi
done


vcf_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/one_kg_vcf"
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/"
for chr in {1..9}; do
    docker run --platform=linux/amd64 \
        -v /Users/JJOHN41/:/Users/JJOHN41/ \
        jibinjv/postgwas:1.3 \
        bash -c "
    bcftools +liftover ${vcf_folder}/GRCh37_1000G_freq_chr${chr}.vcf.gz \
        -Ou --no-version -- \
        --src-fasta-ref ${resourse_folder}/GRCh37/fasta_files/GRCh37_chr${chr}.fa \
        --fasta-ref ${resourse_folder}/GRCh38/fasta_files/GRCh38_chr${chr}.fa \
        --chain ${resourse_folder}/chain_files/GRCh37_to_GRCh38.chain \
        --af-tags INFO/EAS,INFO/AMR,INFO/AFR,INFO/EUR,INFO/SAS \
        --reject ${vcf_folder}/GRCh37_to_38_rejected_1000G_freq_chr${chr}.vcf.gz \
    | bcftools view -e 'INFO/SWAP==1 || INFO/SWAP==-1' \
    | bcftools norm -m -any -d exact \
    | bcftools sort -Oz \
        -o ${vcf_folder}/GRCh38_1000G_freq_chr${chr}.vcf.gz \
        --write-index=tbi
    "
done




## flames
https://zenodo.org/records/12635505
## cadd score
https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz
https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz.tbi
## vep
https://ftp.ensembl.org/pub/grch37/release-113/variation/vep/homo_sapiens_vep_113_GRCh37.tar.gz







## need to remove ALT="." from vcf files before indexing with tabix




# Loop through all VCF.GZ files in the current folder
for FILE in *.vcf.gz; do
    echo "Processing $FILE ..."

    # 1. Run bcftools to remove lines where ALT is "."
    # We check for success (&&) before overwriting to prevent data loss
    if bcftools view -e 'ALT="."' "$FILE" -Oz -o "${FILE}.clean.vcf.gz"; then

        # 2. Replace the original file with the clean one
        mv "${FILE}.clean.vcf.gz" "$FILE"

        # 3. Index using tabix (Force overwrite with -f)
        tabix -p vcf -f "$FILE"

        echo "Successfully cleaned and indexed $FILE"
    else
        echo "ERROR processing $FILE. Skipping."
        # Clean up the partial file if it exists
        rm -f "${FILE}.clean.vcf.gz"
    fi

done




##--------------------------------------------------- MiXeR --------------------------------------------------
https://github.com/comorment/mixer
wget -O folder_name.zip "https://www.dropbox.com/scl/fo/y5yl2bd5mgplsjwwzsx77/AIFIhSJkzJTFIYhR95TwRVc?rlkey=eydtbzwva5294snzgf6lz0g5f&dl=1"


docker pull ghcr.io/precimed/gsa-mixer:2.2.1   # replace 2.2.1 with latest tag without "v"
# -u $(id -u):$(id -g)  -> Runs container as your current user/group
# -e HOME=/tmp          -> Sets temp directory as home to prevent permission errors with Python caches
export DOCKER_RUN="docker run -u $(id -u):$(id -g) -e HOME=/tmp -v /mnt/disks/sdd/:/mnt/disks/sdd/ -w /mnt/disks/sdd/"
export MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest python /tools/mixer/precimed/mixer.py"
base_dir="/mnt/disks/sdd/resourses/mixer_reference/mixer_hello_world/"


## reference
# https://github.com/precimed/mixer/blob/master/usecases/mixer_real.md
# https://github.com/precimed/mixer/blob/master/scripts/GSA_MIXER.job

THREADS=6
OUT_FOLDER=/mnt/disks/sdd/resourses/mixer_reference/sumstats/analysis/
SUMSTATS_FILE1=/mnt/disks/sdd/resourses/mixer_reference/sumstats/INT.sumstats.gz
SUMSTATS_NAME="INT"

REFERENCE_FOLDER="/mnt/disks/sdd/resourses/mixer_reference/"
BIM_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
LOADLIB_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin
ANNOT_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz


for chri in {21..22}
    do
        ${MIXER_PY} ld \
            --bfile ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chri} \
            --r2min 0.05 \
            --ldscore-r2min 0.01 \
            --out ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chri}.ld \
            --ld-window-kb 10000
    done

# generate .bin file for --loadlib-file argument
${MIXER_PY} plsa \
      --bim-file ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --use-complete-tag-indices \
      --chr2use 21-22 \
      --exclude-ranges [] \
      --savelib-file ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
      --out ${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@





## fingen info and af
gsutil -m cp -r "gs://finngen-public-data-r12/annotations" .




## ukbiobank whole genome summstat freq file preparation ;
    # first downloaded GCST90473552_ , postgwas harmonsation step
for i in {1..22} X ; do
    echo "Processing chromosome ${i}..."
    # Generate the file with header and data
    (echo -e "CHROM\tPOS\tREF\tALT\tEUR"; \
     bcftools query -r ${i} -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' GCST90473552_GRCh38_merged.vcf.gz) \
     | bgzip -c > GRCh38_wgs_ukb_freq_chr${i}.tsv.gz
    # Index the file
    tabix -s1 -b2 -e2 -f GRCh38_wgs_ukb_freq_chr${i}.tsv.gz
done


for i in {1..22} X ; do
    echo "Processing chromosome ${i}..."
    # Generate the file with header and data
    (echo -e "CHROM\tPOS\tREF\tALT\tEUR"; \
     bcftools query -r ${i} -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' GCST90473552_GRCh37_merged.vcf.gz) \
     | bgzip -c > GRCh37_wgs_ukb_freq_chr${i}.tsv.gz
    # Index the file
    tabix -s1 -b2 -e2 -f GRCh37_wgs_ukb_freq_chr${i}.tsv.gz
done


gsutil -m cp GRCh37*_wgs_ukb_freq* "gs://software_resourse/postgwas/gwas2vcf/GRCh37/default_af/tab_files"

gsutil -m cp GRCh38_wgs_ukb_freq* "gs://software_resourse/postgwas/gwas2vcf/GRCh38/default_af/tab_files"


gsutil cp  "gs://software_resourse/postgwas/gwas2vcf/GRCh37/default_af/tab_files/GRCh37_wgs_ukb_freq_chr*" .

gsutil cp  "gs://software_resourse/postgwas/gwas2vcf/GRCh38/default_af/tab_files/GRCh38_wgs_ukb_freq_chr*" .





### info score from ukb and fingen merged
import pandas as pd
import os

# Create list: numbers 1 to 22, plus 'X'
chrom_list = list(range(1, 23)) + ['X']

for i in chrom_list:
    print(f"Processing Chr {i}...")
    # Read files
    df_ukb = pd.read_csv(f'GRCh38_panukb_infoscore_chr{i}.tsv.gz', sep="\t", dtype={'CHROM': str})
    df_fin = pd.read_csv(f'GRCh38_fingen_infoscore_chr{i}.tsv.gz', sep="\t", dtype={'CHROM': str})
    # Merge
    merged = pd.merge(df_ukb, df_fin, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer', suffixes=('_ukb', '_fin'))
    # Fill INFO: Use UKB first; if missing (NaN), use Fingen
    merged['INFO'] = merged['INFO_ukb'].fillna(merged['INFO_fin'])
    # Filter and Sort (Header remains "CHROM")
    final_df = merged[['CHROM', 'POS', 'REF', 'ALT', 'INFO']].sort_values(['CHROM', 'POS'])
    # --- SAVE AND INDEX ---
    # Define filenames
    temp_tsv = f'temp_chr{i}.tsv'
    final_gz = f'GRCh38_ukb_fingen_infoscore_chr{i}.tsv.gz'
    # 1. Save as uncompressed TSV first
    final_df.to_csv(temp_tsv, sep="\t", index=False)
    # 2. Compress with bgzip
    os.system(f"bgzip -c {temp_tsv} > {final_gz}")
    # 3. Index with tabix
    # -S1: Skip the first line (Header) because it doesn't start with '#'
    # -s1, -b2, -e2: Define columns for Chrom, Start, End
    # -f: Force overwrite
    os.system(f"tabix -S1 -s1 -b2 -e2 -f {final_gz}")
    # 4. Cleanup
    os.remove(temp_tsv)

print("All chromosomes completed.")



gsutil cp GRCh37_ukb_fingen_infoscore_chr* "gs://software_resourse/postgwas/gwas2vcf/GRCh37/default_infoscore/tab_files"
gsutil cp GRCh38_ukb_fingen_infoscore_chr* "gs://software_resourse/postgwas/gwas2vcf/GRCh38/default_infoscore/tab_files"