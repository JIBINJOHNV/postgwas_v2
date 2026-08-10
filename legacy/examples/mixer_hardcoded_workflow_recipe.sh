


docker pull ghcr.io/precimed/gsa-mixer:2.2.1   # replace 2.2.1 with latest tag without "v"
# -u $(id -u):$(id -g)  -> Runs container as your current user/group
# -e HOME=/tmp          -> Sets temp directory as home to prevent permission errors with Python caches
export DOCKER_RUN="docker run -u $(id -u):$(id -g) -e HOME=/tmp -v /mnt/disks/sdd/:/mnt/disks/sdd/ -w /mnt/disks/sdd/"
export MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest python /tools/mixer/precimed/mixer.py"
base_dir="/mnt/disks/sdd/resourses/mixer_reference/mixer_hello_world/"


## reference 
# https://github.com/precimed/mixer/blob/master/usecases/mixer_real.md
# https://github.com/precimed/mixer/blob/master/scripts/GSA_MIXER.job
# https://github.com/precimed/mixer/tree/master/usecases 


THREADS=6
OUT_FOLDER=/mnt/disks/sdd/resourses/mixer_reference/sumstats/analysis/
SUMSTATS_FILE1=/mnt/disks/sdd/resourses/mixer_reference/sumstats/INT.sumstats.gz
SUMSTATS_NAME1="INT"

SUMSTATS_FILE2=/mnt/disks/sdd/resourses/mixer_reference/sumstats/SCZ.sumstats.gz
SUMSTATS_NAME2="SCZ"

REFERENCE_FOLDER="/mnt/disks/sdd/resourses/mixer_reference/"
BIM_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
LOADLIB_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin
mixxer_ld_file=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
ANNOT_FILE=${REFERENCE_FOLDER}/ldsc/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz

mixxer_extract=${REFERENCE_FOLDER}/ldsc/w_hm3.justrs



EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.$REP.snps"

COMMON_FLAGS=""
COMMON_FLAGS="${COMMON_FLAGS} --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
COMMON_FLAGS="${COMMON_FLAGS} --ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
COMMON_FLAGS="${COMMON_FLAGS} --threads ${SLURM_CPUS_PER_TASK}"
COMMON_FLAGS="${COMMON_FLAGS} --exclude-ranges MHC"

PYTHON="singularity exec --home=pwd:/home $GITHUB/precimed/gsa-mixer/containers/latest/mixer.sif python"
MIXER_PY="$PYTHON         /tools/mixer/precimed/mixer.py"
MIXER_FIGURES_PY="$PYTHON /tools/mixer/precimed/mixer_figures.py"


mkdir -p ${OUT_FOLDER}
# GSA-MiXeR analysis - split sumstats per chromosome
  ${MIXER_PY} split_sumstats \
      --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE1} \
      --out ${OUT_FOLDER}/${SUMSTATS_NAME}.chr@.sumstats.gz \
      --log ${OUT_FOLDER}/${SUMSTATS_NAME} \
      --chr2use 1-22 

  ${MIXER_PY} split_sumstats \
      --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE2} \
      --out ${OUT_FOLDER}/${SUMSTATS_NAME2}.chr@.sumstats.gz \
      --log ${OUT_FOLDER}/${SUMSTATS_NAME2} \
      --chr2use 1-22 



# GSA-MiXeR analysis - baseline
${MIXER_PY} plsa --gsa-base \
        --trait1-file ${OUT_FOLDER}/${SUMSTATS_NAME}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_NAME}_base \
        --use-complete-tag-indices \
        --bim-file ${BIM_FILE} \
        --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --chr2use 1-22  \
        --exclude-ranges chr6:26-34MB \
        --seed 1000 \
        --hardprune-r2 0.6 \
        --threads ${THREADS} \
        --adam-epoch 3 3 \
        --adam-step 0.064 0.032

${MIXER_PY} plsa --gsa-base \
        --trait1-file ${OUT_FOLDER}/${SUMSTATS_NAME2}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_NAME2}_base \
        --use-complete-tag-indices \
        --bim-file ${BIM_FILE} \
        --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --chr2use 1-22  \
        --exclude-ranges chr6:26-34MB \
        --seed 1000 \
        --hardprune-r2 0.6 \
        --threads ${THREADS} \
        --adam-epoch 3 3 \
        --adam-step 0.064 0.032


${MIXER_PY} plsa --gsa-full \
        --trait1-file ${OUT_FOLDER}/${SUMSTATS_NAME}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_NAME}_full \
        --use-complete-tag-indices \
        --bim-file ${BIM_FILE} \
        --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-annot_10mar2023.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --load-params-file ${OUT_FOLDER}/${SUMSTATS_NAME}_base.json \
        --chr2use 1-22  \
        --exclude-ranges chr6:26-34MB \
        --seed 1000 \
        --hardprune-r2 0.6 \
        --threads ${THREADS} \
        --adam-epoch 3 3 \
        --adam-step 0.064 0.032


${MIXER_PY} plsa --gsa-full \
        --trait1-file ${OUT_FOLDER}/${SUMSTATS_NAME2}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_NAME2}_full \
        --use-complete-tag-indices \
        --bim-file ${BIM_FILE} \
        --loadlib-file ${LOADLIB_FILE} \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_10mar2023.csv \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-annot_10mar2023.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_10mar2023.csv \
        --annot-file ${ANNOT_FILE} \
        --load-params-file ${OUT_FOLDER}/${SUMSTATS_NAME2}_base.json \
        --chr2use 1-22  \
        --exclude-ranges chr6:26-34MB \
        --seed 1000 \
        --hardprune-r2 0.6 \
        --threads ${THREADS} \
        --adam-epoch 3 3 \
        --adam-step 0.064 0.032




## https://github.com/precimed/mixer/blob/master/usecases/mixer_real/MIXER_REAL.job
## https://github.com/precimed/mixer/blob/master/scripts/MIXER.job 

${MIXER_PY} fit1 \
    --trait1-file ${SUMSTATS_FILE1} \
    --out ${OUT_FOLDER}/${SUMSTATS_NAME1}_univar \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --threads ${THREADS} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 

${MIXER_PY} fit1 \
    --trait1-file ${SUMSTATS_FILE2} \
    --out ${OUT_FOLDER}/${SUMSTATS_NAME2}_univar \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --threads ${THREADS} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 

${MIXER_PY} test1 \
    --trait1-file ${SUMSTATS_FILE1} \
    --out ${OUT_FOLDER}/${SUMSTATS_NAME1}_univar_test \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 \
    --load-params-file ${OUT_FOLDER}/${SUMSTATS_NAME1}_univar.json

${MIXER_PY} test1 \
    --trait1-file ${SUMSTATS_FILE2} \
    --out ${OUT_FOLDER}/${SUMSTATS_NAME2}_univar_test \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 \
    --load-params-file ${OUT_FOLDER}/${SUMSTATS_NAME2}_univar.json


${MIXER_PY} fit2 \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 \
    --trait1-file SCZ.sumstats.gz \
    --trait2-file INT.sumstats.gz \
    --trait1-params SCZ.fit.$REP.json \
    --trait2-params INT.fit.$REP.json \
    --out SCZ_vs_INT.fit.$REP

$MIXER_PY test2 \
    --bim-file ${BIM_FILE} \
    --ld-file ${mixxer_ld_file} \
    --extract ${mixxer_extract} \
    --chr2use 1-22 \
    --trait1-file SCZ.sumstats.gz \
    --trait2-file INT.sumstats.gz \
    --load-params SCZ_vs_INT.fit.$REP.json \
    --out SCZ_vs_INT.test.$REP

















# fit baseline model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-base \
    --trait1-file ${baser_dir}/trait1.chr@.sumstats.gz \
    --use-complete-tag-indices \
    --bim-file ${baser_dir}/g1000_eur_hm3_chr@.bim \
    --loadlib-file ${baser_dir}/g1000_eur_hm3_chr@.bin \
    --annot-file ${baser_dir}/g1000_eur_hm3_chr@.annot.gz \
    --go-file ${baser_dir}/go-file-baseline.csv \
    --exclude-ranges chr21:20-21MB chr22:19100-19900KB \
    --chr2use 21-22 \
    --seed 123 \
    --adam-epoch 3 3 \
    --adam-step 0.064 0.032 \
    --out ${baser_dir}/plsa_base




for chri in {21..22}
    do ${MIXER_PY} ld \
        --bfile ${baser_dir}/g1000_eur_hm3_chr$chri \
        --r2min 0.05 \
        --ldscore-r2min 0.01 \
        --out ${baser_dir}/g1000_eur_hm3_chr$chri.ld \
        --ld-window-kb 10000
    done  



# split summary statistics into one file per chromosome
${MIXER_PY} split_sumstats \
    --trait1-file ${baser_dir}/trait1.sumstats.gz \
    --out ${baser_dir}/trait1.chr@.sumstats.gz \
    --chr2use 21-22


# generate .bin file for --loadlib-file argument
${MIXER_PY} plsa \
      --bim-file ${baser_dir}/g1000_eur_hm3_chr@.bim \
      --ld-file ${baser_dir}/g1000_eur_hm3_chr@.ld \
      --use-complete-tag-indices \
      --chr2use 21-22 \
      --exclude-ranges [] \
      --savelib-file ${baser_dir}/g1000_eur_hm3_chr@.bin \
      --out ${baser_dir}/g1000_eur_hm3_chr@


# fit baseline model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-base \
    --trait1-file ${baser_dir}/trait1.chr@.sumstats.gz \
    --use-complete-tag-indices \
    --bim-file ${baser_dir}/g1000_eur_hm3_chr@.bim \
    --loadlib-file ${baser_dir}/g1000_eur_hm3_chr@.bin \
    --annot-file ${baser_dir}/g1000_eur_hm3_chr@.annot.gz \
    --go-file ${baser_dir}/go-file-baseline.csv \
    --exclude-ranges chr21:20-21MB chr22:19100-19900KB \
    --chr2use 21-22 \
    --seed 123 \
    --adam-epoch 3 3 \
    --adam-step 0.064 0.032 \
    --out ${baser_dir}/plsa_base


${MIXER_PY} plsa --gsa-full \
    --trait1-file ${baser_dir}/trait1.chr@.sumstats.gz \
    --use-complete-tag-indices \
    --bim-file ${baser_dir}/g1000_eur_hm3_chr@.bim \
    --loadlib-file ${baser_dir}/g1000_eur_hm3_chr@.bin \
    --annot-file ${baser_dir}/g1000_eur_hm3_chr@.annot.gz \
    --go-file ${baser_dir}/go-file-gene.csv \
    --go-file-test ${baser_dir}/go-file-geneset.csv \
    --load-params-file ${baser_dir}/plsa_base.json \
    --exclude-ranges chr21:20-21MB chr22:19100-19900KB \
    --chr2use 21-22 \
    --seed 123 \
    --adam-epoch 3 3 \
    --adam-step 0.064 0.032 \
    --out ${baser_dir}/plsa_full 


${MIXER_PY} fit1 \
    --trait1-file ${base_dir}/trait1.sumstats.gz \
    --out ${base_dir}/trait1_name_univar \
    --bim-file ${base_dir}/g1000_eur_hm3_chr@.bim \
    --ld-file ${base_dir}/g1000_eur_hm3_chr@.ld \
    --chr2use 21-22 

${MIXER_PY} test1 \
    --trait1-file ${base_dir}/trait1.sumstats.gz \
    --out ${base_dir}/trait1_name_univar_test \
    --bim-file ${base_dir}/g1000_eur_hm3_chr@.bim \
    --ld-file ${base_dir}/g1000_eur_hm3_chr@.ld \
    --chr2use 21-22 \
    --load-params-file ${base_dir}/trait1_name_univar.json



${MIXER_PY} test1 \
    --trait1-file ${base_dir}/trait1.sumstats.gz \
    --out ${base_dir}/trait1_name_univar_test \
    --bim-file ${base_dir}/g1000_eur_hm3_chr@.bim \
    --ld-file ${base_dir}/g1000_eur_hm3_chr@.ld \
    --chr2use 21-22 \
    --load-params-file ${base_dir}/trait1_name_univar.json \
    --power-curve \
    --qq-plots


# Locate the figures script (usually in the same folder as mixer.py)
$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest \
    python /tools/mixer/precimed/mixer_figures.py one \
    --json ${base_dir}/trait1_name_univar_test.json \
    --out ${base_dir}/trait1_name_univar_plot \
    --trait1 "My Trait"