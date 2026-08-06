# PostGWAS

PostGWAS is a command-line toolkit that provides a collection of modules commonly used in post-GWAS analyses.
It supports tasks such as LD-block annotation, fine-mapping, summary-statistic filtering, imputation, MAGMA gene analysis,
PoPS scoring, and more. Each module is available as a standalone command, and a unified **pipeline** interface allows
chaining multiple steps together.

## Available Modules

```
┏━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Command               ┃ Description                               ┃
┡━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ annot_ldblock         │ LD-block annotation                       │
│ finemap               │ Fine-mapping (SuSiE / FINEMAP)            │
│ flames                │ FLAMES integrative scoring                │
│ formatter             │ Convert VCF to tool-specific formats      │
│ harmonisation         │ Harmonisation of input summary statistics │
│ heritability          │ LDSC-based heritability estimation        │
│ imputation            │ Summary-statistic imputation              │
│ ld_clump              │ LD clumping                               │
│ magma                 │ MAGMA gene/pathway analysis               │
│ magmacovar            │ MAGMA gene-property model                 │
│ manhattan             │ Manhattan/QQ plot generation              │
│ pipeline              │ Multi-step summary-statistic pipeline     │
│ pops                  │ PoPS gene-prioritisation                  │
│ qc                    │ QC summary reporting                      │
│ sumstat_filter        │ Summary-statistic filtering               │  
| pathway_enrichmnet    │ Run pathway enrichment module             │
└───────────────────────┴───────────────────────────────────────────┘
```

## Installation

Highly recomend using the provided Docker image for ease of use and reproducibility.

```sh
docker pull jibinjv/postgwas:1.3

docker run --rm --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -it jibinjv/postgwas:1.3 postgwas --help
```

Display help for any module:

```sh
docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -it jibinjv/postgwas:1.3 postgwas finemap --help
```

Display help for pipeline module:

```sh
docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \ 
    -it jibinjv/postgwas:1.3 postgwas pipeline --modules flames --help

```

Along with main pipleine module users can use ptional flags:

- --apply-filter
- --apply-imputation
- --apply-manhattan
- --heritability

Example:

```sh
docker run --rm  --platform=linux/amd64 \
     -u $(id -u):$(id -g) \
     jibinjv/postgwas:1.3 postgwas pipeline \
    --modules flames \
    --apply-filter \
    --heritability \
    --apply-imputation --help
```

# Test Run 

- This run make sure the docker

```sh
git clone https://github.com/JIBINJOHNV/postgwas.git
cd postgwas/tests/ 

python test_run.py \
    --input-dir `pwd` \
    --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/  


```

# Tutorials

- It is **highly recommended** to use the Docker image instead of a local installation.  
  Installing all dependencies manually is complex and error-prone, and is likely to fail.

- Please mount a local directory into the Docker container so that PostGWAS can read input files from, and write output files to, your local system.  
  Refer to the official Docker documentation for bind mounts:  
  https://docs.docker.com/engine/storage/bind-mounts/

- When mounting a local directory, make sure that the user inside the Docker container has permission to read from and write to the mounted directory.  
  Otherwise, you may encounter **permission denied** errors during execution.

- To avoid path confusion, mount the local directory into the container using the **same absolute path structure** as on your local system.  
  For example, if your local directory is:
    `/home/user/data/postgwas/`, mount it into the container as: /Users/username/data/  
- If  not mounted all path , file not found error may occur.

#### <span style="color: #5fa8ff;"><strong>Step 1: Prepare input config file for the PostGWAS harmonisation module</strong></span>

- Test files can be downloded from the `tests` folder in this repository: https://github.com/JIBINJOHNV/postgwas/tree/main/tests

```sh
docker run --rm --platform=linux/amd64 \
  -u $(id -u):$(id -g) \
  -v /Users/JJOHN41/:/Users/JJOHN41/ \
  -it jibinjv/postgwas:1.3 python /opt/postgwas/src/postgwas/scripts/create_sumstat_map_pl.py \
  --input /Users/JJOHN41/Downloads/postgwas/tests \
  --resource-folder /Users/JJOHN41/Documents/software_resources/resourses/postgwas/gwas2vcf/ \
  --output-path /Users/JJOHN41/Downloads/postgwas/tests/harmonisatio_example_input_file.csv \
  --harmonisation-output-path /Users/JJOHN41/Downloads/postgwas/tests/


# --input Full Path to folder containing summary statistics files
# --resource-folder Path to postgwas resource folder containing reference data files
# --output-path Full path to output file  with file name
# --harmonisation-output-path Full path to the base directory for harmonised outputs (a separate subfolder will be created for each GWAS).
```

####  It is highly recommended to open and review the generated output file `harmonisatio_example_input_file.csv` to ensure that all GWAS summary-statistic fields are correctly mapped before proceeding to Step 2.

- For each GWAS summary-statistic file in the input folder, the following information will be automatically detected and mapped into the output CSV file:


---

### Column Descriptions

| Column Name        | Description |
|--------------------|-------------|
| `sumstat_file`     | Full path to the original GWAS summary-statistics file. |
| `gwas_outputname`  | Base name used to label outputs for this GWAS dataset. |
| `chr_col`          | Column name in the sumstats file representing chromosome. |
| `pos_col`          | Column name representing genomic position. |
| `snp_id_col`       | Column name representing SNP identifier (e.g., rsID). |
| `ea_col`           | Column name for the effect allele. |
| `oa_col`           | Column name for the other (non-effect) allele. |
| `eaf_col`          | Column name for effect allele frequency (EAF). |
| `beta_or_col`      | Column name for effect size (beta or odds ratio). |
| `se_col`           | Column name for standard error of the effect size. |
| `imp_z_col`        | Column name for imputed Z-score (if available). |
| `pval_col`         | Column name for p-value. |
| `ncontrol_col`    | Column name for number of controls (if present in file). |
| `ncase_col`       | Column name for number of cases (if present in file). |
| `ncontrol`        | Fixed number of controls (used if not present in the file). |
| `ncase`           | Fixed number of cases (used if not present in the file). |
| `imp_info_col`    | Column name for imputation quality score (INFO). |
| `infofile`        | External file containing INFO scores (if not in sumstats). |
| `infocolumn`      | Column name in the external INFO file. |
| `eaffile`         | External file containing EAF values (if not in sumstats). |
| `eafcolumn`       | Column name in the external EAF file. |
| `liftover`        | Whether genomic liftover should be applied (`Yes` / `No`). |
| `chr_pos_col`     | Column representing combined chromosome:position (if present). |
| `resource_folder`| Path to PostGWAS reference resources (e.g., gwas2vcf, liftover chains). |
| `output_folder`   | Output directory where harmonised files for this GWAS will be written. |

---

#### After verifying and (if needed) editing the values in `harmonisatio_example_input_file.csv`, proceed to Step 2 to run the PostGWAS harmonisation module.

If your GWAS summary-statistics file does **not** contain the columns `ncase_col` and `ncontrol_col`:

- **For a binary (case–control) phenotype**  
  Open the file `harmonisatio_example_input_file.csv` and:
  - Enter the **total number of cases** in the `ncase` column  
  - Enter the **total number of controls** in the `ncontrol` column  

- **For a quantitative (continuous) phenotype**  
  Open the file `harmonisatio_example_input_file.csv` and:
  - Leave the `ncase` column as `NA`  
  - Enter the **total sample size** in the `ncontrol` column  

---

If your GWAS summary-statistics file does **not** contain an effect allele frequency column (`eaf_col`):

- Open `harmonisatio_example_input_file.csv` and manually fill:
  - `eaffile` → Path to an external allele-frequency file  
  - `eafcolumn` → Name of the column in that file containing allele frequencies  

- The external EAF file **must** contain the following columns:  
  `CHROM`, `POS`, `REF`, `ALT`  
  (plus one additional column holding the allele frequency, e.g., `EAF`)

- **If the EAF data are stored in a single file for all chromosomes:**  
  - Provide the **full file path** in `eaffile`  

- **If the EAF data are stored in separate per-chromosome files:**  
  - Provide the **base file path without the chromosome suffix**  
  - The harmonisation module will automatically append the chromosome-specific suffix  
    (for example: `_chr15.tsv.gz`) when processing each chromosome  
  - ⚠️ **Important:**  
  - Per-chromosome files **must** follow this exact naming pattern.  
  - Only files ending in `_chr<chromosome>.tsv.gz` are supported.  
  - Any other naming convention (e.g., `.chr15.txt`, `_15.tsv.gz`, `chr15.txt.gz`) will **not** be recognized by the harmonisation module.
---

If your GWAS summary-statistics file does **not** contain an imputation INFO column (`imp_info_col`):

- Open `harmonisatio_example_input_file.csv` and manually fill:
  - `infofile` → Path to an external INFO-score file  
  - `infocolumn` → Name of the column in that file containing INFO scores  

- The external INFO file **must** contain the following columns:  
  `CHROM`, `POS`, `REF`, `ALT`  
  (plus one additional column holding the INFO score, e.g., `INFO`)

- **If the INFO data are stored in a single file for all chromosomes:**  
  - Provide the **full file path** in `infofile`  

- **If the INFO data are stored in separate per-chromosome files:**  
  - Provide the **base file path without the chromosome suffix**  
  - The harmonisation module will automatically append the chromosome-specific suffix  
  - ⚠️ **Important:**  
  - Per-chromosome files **must** follow this exact naming pattern.  
  - Only files ending in `_chr<chromosome>.tsv.gz` are supported.  
  - Any other naming convention (e.g., `.chr15.txt`, `_15.tsv.gz`, `chr15.txt.gz`) will **not** be recognized by the harmonisation module.

---

After reviewing and (if needed) editing the values in `harmonisatio_example_input_file.csv`, proceed to **Step 2** to run the PostGWAS harmonisation module.

#### <span style="color: #5fa8ff; font-weight: bold;">Step 2: Run the PostGWAS harmonisation module</span>
```
The harmonisation.yaml file contains several reference file paths that must be modified to reflect the correct locations on your system.

```

```sh

docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas harmonisation \
    --nthreads 10 \
    --max-mem 50G \
    --config /Users/JJOHN41/Downloads/postgwas/tests/harmonisatio_example_input_file.csv \
    --defaults /Users/JJOHN41/Downloads/postgwas/tests/harmonisation.yaml


# --nthreads Number of threads to use
# --max-mem Maximum memory to use
# --config Full path to the input config file generated in Step 1
# --defaults Full path to a YAML file containing default parameters for the harmonisation module , example available in the tests folder. https://github.com/JIBINJOHNV/postgwas/tree/main/tests
```

#### <span style="color: #5fa8ff; font-weight: bold;">Step 3: Run PostGWAS Downstream modules individually</span>

Once the harmonisation step is successfully completed, you can proceed to run other PostGWAS modules independently using the harmonised output files as input.

- Please make sure if you are using imputation module use Grch38 sumstat vcf file generated from harmonisation step. 
  - After imputation postgwas harmonsation steps will run automatically and generate the harmonised vcf file for downstream analysis.
  - All other downstream modules will you can use the harmonised vcf file generated from imputation step.
  
- All other downstream modules can use the harmonised Grch37 VCF file generated from the harmonisation step. 
  - Make sure the genome build of the input summary-statistics VCF file matches the genome build of all reference files using in the annot_ldblock magma,pops, flames modules.
- If you are using magma output for pops, flames modules please make sure the gene loc file should contain Ensembl gene ids instead of gene symbols.
- Please ensure that the gene locus file and the gene set file use the same gene identifiers or gene symbols.


Each module follows a predefined dependency flow, as shown below.  
Modules marked as **(opt)** are optional and can be skipped depending on your analysis needs.

---

### 📌 Workflow Targets and Dependency Flow

| Target Module        | Dependency Flow |
|----------------------|------------------|
| **annot_ldblock**   | harmonisation → annot_ldblock  |
| **sumstat_filter** | harmonisation → sumstat_filter |
| **imputation**      | harmonisation → sumstat_filter(opt) → formatter → imputation |
| **manhattan**       | harmonisation → sumstat_filter(opt) → imputation(opt) → manhattan |
| **qc**               | harmonisation → sumstat_filter(opt) → imputation(opt) → qc |
| **formatter**        | harmonisation → sumstat_filter(opt) → imputation(opt) → formatter |
| **ld_clump**        | harmonisation → imputation(opt) → annot_ldblock → sumstat_filter(opt) → ld_clump |
| **heritability**    | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → heritability |
| **magma**            | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → formatter → magma |
| **magmacovar**       | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → formatter → magma → magmacovar|
| **pops**             | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → formatter → magma → pops |
| **finemap**          | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → annot_ldblock → ld_clump → formatter → finemap |
| **flames**           | harmonisation → sumstat_filter(opt) → imputation(opt) → sumstat_filter(opt) → annot_ldblock → ld_clump → formatter → magma → magmacovar→ pops → finemap → flames |
| **pathway_enrichment** | Not depend on any other postgwas module, it require a file that contain gene symbol in the first column |

---

### 🧠 Notes

- **harmonisation** is always the required first step for all downstream modules.  
- **sumstat_filter (opt)** removes low-quality or extreme variants before downstream analysis.  
- **imputation (opt)** is only needed when number of markers for the analyis is less.  
- You may run any module independently as long as all upstream dependencies are satisfied.  

---

After completing Step 3, you can chain modules programmatically or execute them individually depending on your analysis workflow. Or you can use the main `pipeline` module to run multiple steps in a single command. 

- To know all available module options and their usage, run:

```
    docker run --platform=linux/amd64 \
        -v /Users/JJOHN41/:/Users/JJOHN41/ \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -it jibinjv/postgwas:1.3 postgwas --help   
```

```sh
base_dir="/Users/JJOHN41/Downloads/postgwas/tests/"
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"
sample_id="ADHD2022_iPSYCH_deCODE_PGC"

## Please ensure that the genome build of the input summary-statistics VCF file matches the genome build of all reference files.
## Please ensure that the gene locus file and the gene set file use the same gene identifiers or gene symbols.

docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas formatter \
    --nthreads 2 \
    --max-mem 16G \
    --vcf ${base_dir}/${sample_id}/00_harmonised_sumstat/${sample_id}_GRCh37_merged.vcf.gz \
    --format magma \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/${sample_id}/magma_inputs/


docker run --rm  --platform=linux/amd64 \
    -u $(id -u):$(id -g) \
    -v /Users/JJOHN41/:/Users/JJOHN41/ \
    -it jibinjv/postgwas:1.3 postgwas magma \
    --nthreads 2 \
    --max-mem 16G \
    --sample_id ${sample_id} \
    --outdir ${base_dir}/${sample_id}/magma_analysis/ \
    --snp_loc_file ${base_dir}/${sample_id}/${sample_id}_magma_snp_loc.tsv \
    --pval_file ${base_dir}/${sample_id}/${sample_id}_magma_P_val.tsv \
    --ld_ref ${resourse_folder}/onekg_plinkfiles/GRCh37/EUR.chr1_22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_multiallele_uniqid_Grch37_maf0001 \
    --gene_loc_file ${resourse_folder}/magma/gene_loc/NCBI37.3/NCBI37.3.gene_withGeneSymbol.loc \
    --geneset_file ${resourse_folder}/msigdb_v2025.1.Hs_GMTs/msigdb.v2025.1.Hs.symbols.gmt \
    --window_upstream 35 \
    --window_downstream 10 \
    --gene_model snp-wise=mean \
    --n_sample_col N_COL



```

#### <span style="color: #5fa8ff; font-weight: bold;">Step 4: Run PostGWAS pipeline module </span>

```sh
base_dir="/Users/JJOHN41/Downloads/postgwas/tests/"
resourse_folder="/Users/JJOHN41/Documents/software_resources/resourses/postgwas/"
genome_version="GRCh37"
sample_id="ADHD2022_iPSYCH_deCODE_PGC"

docker run --rm --platform=linux/amd64 \
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
        --vcf ${base_dir}/${sample_id}/00_harmonised_sumstat/${sample_id}_GRCh38_merged.vcf.gz \
        --sample_id ${sample_id} \
        --outdir ${base_dir}/${sample_id}/flames_pipeline/ \
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
        --population EUR \
        --corr_method pearson \
        --r2threshold 0.8 \
        --gwas2vcf_resource ${resourse_folder}/gwas2vcf/ \
        --gwas2vcf_default_config ${base_dir}/harmonisation.yaml \
        --imputation_tool pred_ld \
        --ref_ld ${resourse_folder}/imputation/pred-ld/ref/ 
```

## Notes

- For research use only.
- Ensure genome builds match reference data.