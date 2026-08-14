# Fine-mapping examples

These examples use SuSiE-RSS. Replace `STUDY` and the example paths with your
files, but do not mix genome builds, populations, or variant-ID conventions.
The PLINK LD panel must represent the ancestry of the GWAS and its `.bed`,
`.bim`, and `.fam` files must share one prefix.

## Direct mode

Prepare a tab-separated SuSiE table with exactly these required columns:

| Column | Meaning |
|---|---|
| `SNP` | Variant identifier that exists in the PLINK `.bim`. |
| `CHR`, `BP` | Genome-build-matched chromosome and 1-based position. |
| `REF`, `ALT` | Non-effect and effect alleles; `EZ` is for `ALT`. |
| `EZ`, `LP` | Z score and `-log10(P)`. |
| `NEF` | Positive effective sample size. |

The locus TSV must contain `CHR`, `START`, and `END`; `GenomicLocus` is
recommended for stable provenance. The configuration uses `locus_window_kb: 0`
because these are already the intended fixed model boundaries.

```bash
postgwas finemap \
  --run-config examples/fine_mapping/direct_susie.yaml \
  --finemap-method susie \
  --susie-input-file formatted/STUDY_susie.tsv.gz \
  --locus-file loci.tsv \
  --finemap-ld-reference reference/1000G_EUR \
  --plink plink \
  --dataset-id STUDY \
  --output-directory results/direct_susie
```

## Pipeline mode

Pipeline mode starts from a harmonised, bgzip-compressed, tabix-indexed
GWAS-VCF and executes LD-block annotation, formatter export, standard LD
clumping, and fine-mapping. The GWAS-VCF must contain `ES`, `SE`, `EZ`, `LP`,
`AF`, and `NEF` FORMAT values for one study sample.

The example expects:

- `reference/ld_blocks/GRCh37_EUR_ldetect.bed.gz` plus its tabix index. The BED
  columns are chromosome, start, end, and block identifier.
- `reference/pairwise_ld/EUR_chr1.ld.gz` through the analysed chromosomes, each
  with a tabix index and seven columns: chromosome A, position A, variant A,
  chromosome B, position B, variant B, and r-squared. The current clumper
  queries the tabix index by the variant-A coordinate. Every potential index
  SNP must therefore occur as variant A, for example by storing symmetric
  pairwise rows. Otherwise that SNP is explicitly logged as self-only and its
  inferred locus may be narrower than intended.
- `reference/1000G_EUR.{bed,bim,fam}` for locus LD construction.
- One consistent unique identifier convention across VCF, pairwise-LD tables,
  and PLINK BIM. This command uses
  `{chromosome}_{position}_{reference_allele}_{alternate_allele}`; use
  `--variant-id-type rsid` instead when every relevant file uses rsIDs.

```bash
postgwas pipeline \
  --modules finemap \
  --run-config examples/fine_mapping/pipeline_susie.yaml \
  --vcf study_GRCh37.vcf.gz \
  --genome-build GRCh37 \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --ld-folder reference/pairwise_ld \
  --population EUR \
  --variant-id-type unique \
  --finemap-method susie \
  --finemap-ld-reference reference/1000G_EUR \
  --plink plink \
  --bcftools bcftools \
  --dataset-id STUDY \
  --output-directory results/pipeline_susie
```

Validate either YAML before a long run:

```bash
postgwas config validate --config examples/fine_mapping/direct_susie.yaml
postgwas config validate --config examples/fine_mapping/pipeline_susie.yaml
```

Successful runs finish with `overall status: success`. Inspect
`results/combined_results/final_combined_credible_sets.tsv`,
`quality_control/overlap_resolution_summary.tsv`, the engine QC summary under
`quality_control/`, and the resolved configuration under `run_metadata/`
before downstream analysis. The authoritative FLAMES handoff is
`downstream_inputs/flames/`. A converged run with no credible set is a valid
analysis outcome, not a successful FLAMES handoff.
