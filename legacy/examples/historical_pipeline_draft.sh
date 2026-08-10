python /Users/JJOHN41/Documents/developing_software/postgwas/src/postgwas/finemap/cli.py \
    --sample_id  AllOA \
    --outdir  /Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/AllOAAllOA/00_harmonised_sumstat/ \
    --finemap_method susie \
    --locus_file AllOA_GenomicRiskLoci_Summary.txt \
    --locus_type point \
    --window_kb 500 \
    --susie_input_file AllOA_finemap_susie.tsv \
    --finemap_ld_ref  /Users/JJOHN41/Documents/software_resources/resourses/postgwas/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 


docker run --platform=linux/amd64 \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 bash  




sample_id="AllOA"
base_dir='/Users/JJOHN41/Downloads/osteoarthritis_erectile_dysfunction/harmonisation/'
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"



docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
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
        --ld-folder ${resourse_folder}/onekg_plinkfiles/GRCh37/LD_ref_EUR/ \
        --finemap_method finemap \
        --sss \
        --locus_type point \
        --window_kb 500 \
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








docker run --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 \
  postgwas pipeline \
        --modules finemap \
        --nthreads 10 \
        --max-mem 60G \
        --seed 10 \
        --vcf ${base_dir}/${sample_id}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/${sample_id}/ \
        --finemap_ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
        --ld-region-dir ${resourse_folder}ld_blocks/ \
        --ld_clump_population EUR \
        --ld-folder ${resourse_folder}/onekg_plinkfiles/GRCh37/LD_ref_EUR/ \
        --locus_type point \
        --window_kb 500 \
        --finemap_method finemap \
        --sss 


