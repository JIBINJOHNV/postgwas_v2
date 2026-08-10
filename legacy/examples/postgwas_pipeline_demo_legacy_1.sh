
docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 python /opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py --help 





## convert summstat to vcf and harmonise
docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.4 postgwas harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Documents/developing_software/postgwas/tests/example_input_file.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml




docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 python /opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py \
  --input /Users/JJOHN41/Downloads/ukb_ppp_testing/raw_data/ \
  --output-path /Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/ukb_ppp_testing_input_susmstat.csv \
  --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
  --harmonisation-output-path /Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/





docker run --platform=linux/amd64 \
  -v /mnt/fast/:/mnt/fast/ \
  -it jibinjv/postgwas:1.3 bash 

sample_id="ukb_ppp_testing"
base_dir='/Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"

## convert summstat to vcf and harmonise
docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
        --nthreads 10 \
        --max-mem 50G \
        --config /Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/ukb_ppp_testing_input_susmstat.csv \
        --defaults /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml




docker run --platform=linux/amd64 \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 postgwas sumstat_filter \
  --vcf /Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/ITGB2_P05107_OID20315_v1_Cardiometabolic/00_harmonised_sumstat/ITGB2_P05107_OID20315_v1_Cardiometabolic_GRCh37_merged.vcf.gz \
  --sample_id ITGB2_P05107_OID20315_v1_Cardiometabolic \
  --outdir /Users/JJOHN41/Downloads/ukb_ppp_testing/vcf_files/ITGB2_P05107_OID20315_v1_Cardiometabolic/01_filtered_vcf \
  --maf-cutoff 0.01 \
  --nthreads 10 \
  --external-af-name EUR \
  --allelefreq-diff-cutoff 0.2 \
  --info-cutoff 0.7 \
  --exclude-palindromic \
  --remove-mhc \
  --include-indels


docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -it jibinjv/postgwas:1.3 \
  postgwas pipeline \
        --modules flames \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --apply-manhattan \
        --apply-filter \
        --apply-imputation \
        --heritability \
        --vcf ${base_dir}/${sample_id}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/${sample_id}/ \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --ld_clump_population EUR \
        --finemap_method finemap \
        --sss \
        --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
        --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
        --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
        --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \
        --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
        --maf 0.001 \
        --merge-alleles /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
        --ref-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --w-ld-chr /Users/JJOHN41/Documents/software_resources/resourses/postgwas/1000GP_Phase3/eur_w_ld_chr/ \
        --samp-prev 0.5 \
        --pop-prev 0.01 \
        --heritability_info-min 0.7 \
        --heritability_maf-min 0.01 \
        --heritability_tool ldsc \
        --heritability \
        --finemap_method finemap \
        --population EUR \
        --corr_method pearson \
        --r2threshold 0.8 \
        --gwas2vcf_resource /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
        --gwas2vcf_default_config /Users/JJOHN41/Documents/developing_software/postgwas/tests/harmonisation.yaml \
        --imputation_tool pred_ld \
         --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ 



sample_id="PGC3_SCZ_wave3"
base_dir='/Users/JJOHN41/Downloads/sumstat_harmonisation/raw_sumstat/Harmonised_results/PGC3_SCZ_wave3/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"




docker run --platform=linux/amd64  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 postgwas pipeline \
  --modules magma \
  --nthreads 10 \
  --max-mem 60G \
  --vcf ${base_dir}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
  --sample_id ${sample_id} \
  --outdir ${base_dir}/ \
  --gene_loc_file  ${resourse_folder}/magma/gene_loc/NCBI37.3/NCBI37.3.gene_withGeneSymbol.loc \
  --geneset_file ${resourse_folder}/msigdb_v2025.1.Hs_GMTs/msigdb.v2025.1.Hs.symbols.gmt \
  --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \















docker run --platform=linux/amd64   -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/Documents:/Users/JJOHN41/Documents \
  -it jibinjv/postgwas:1.3 postgwas pathway_enrichmnet \
  --sample_id ${sample_id} \
  --outdir ${base_dir}/enrichement_analysis \
  --gene_inputfile /Users/JJOHN41/Documents/developing_software/data/oudir/gene_list.tsv \
  --biogrid-key 86cf1c2c7c8b2972cc5e6b31c45e8584 \
  --david-email jjohn41@northwell.edu \
  --dsigdb-gmt /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/enrichment/DSigDB_All.gmt 



## GCP 

for file in *_GRCh37_merged.vcf.gz; do
  sample_id="${file%_GRCh37_merged.vcf.gz}"
  base_dir='/mnt/disks/sdd/postgwas_analysis/'
  resourse_folder="/mnt/disks/sdd/resourses/postgwas/"
  genome_version="GRCh37"

  echo "Now runing ${sample_id}"

  docker run --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /mnt/disks/sdd/:/mnt/disks/sdd/ \
    --rm jibinjv/postgwas:1.3 \
    postgwas pipeline \
          --modules flames \
          --nthreads 10 \
          --max-mem 60G \
          --seed 10 \
          --apply-manhattan \
          --apply-filter \
          --heritability \
          --vcf ${base_dir}/GRCh37_vcf_files/${sample_id}_GRCh37_merged.vcf.gz \
          --sample_id ${sample_id} \
          --outdir ${base_dir}/analysis_results/${sample_id}/ \
          --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
          --ld-region-dir ${resourse_folder}ld_blocks/ \
          --ld_clump_population EUR \
          --finemap_method finemap \
          --sss \
          --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
          --gene_loc_file  ${resourse_folder}/pops/GRCh37_gene_annot_jun10.loc \
          --covariates  ${resourse_folder}/magma/covar_files/gtex_v8_ts_avg_log2TPM.txt \
          --feature_mat_prefix ${resourse_folder}/pops/features_munged/pops_features \
          --pops_gene_loc_file ${resourse_folder}/pops/GRCh37_gene_annot_jun10.txt \
          --flames_annot_dir ${resourse_folder}/flames/Annotation_data/ \
          --merge-alleles ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/w_hm3.snplist \
          --ref-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
          --w-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
          --samp-prev 0.5 \
          --pop-prev 0.01 \
          --heritability_info-min 0.7 \
          --heritability_maf-min 0.01 \
          --heritability_tool ldsc \
          --heritability \
          --finemap_method finemap \
          --ref-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
          --w-ld-chr ${resourse_folder}/1000GP_Phase3/eur_w_ld_chr/ \
          --vep_cache ${resourse_folder}/ \
          --cmd_vep /opt/conda/envs/vep/bin/vep \
          --tabix /opt/conda/envs/postgwas/bin/tabix \
          --CADD_file ${resourse_folder}/whole_genome_SNVs_GRCh37.tsv.gz 

done 





          --apply-imputation \
          --imputation_tool pred_ld \
          --ref_ld /Users/JJOHN41/Documents/software_resources/resourses/postgwas/imputation/pred-ld/ref/ \
          --r2threshold 0.8 \
          --maf 0.001 \
          --population EUR \
          --corr_method pearson \
          --gwas2vcf_resource ${resourse_folder}/gwas2vcf/ \
          --gwas2vcf_default_config /mnt/disks/sdd/postgwas_analysis/harmonisation.yaml \









