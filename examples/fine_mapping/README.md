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
GWAS-VCF and executes standard LD clumping, formatter export, and fine-mapping.
The GWAS-VCF must contain `ES`, `SE`, `EZ`, `LP`, `AF`, and `NEF` FORMAT values
for one study sample. Standard clumping does not require LD-block annotation;
that preceding step is added only when `region` is also selected.

The example expects:

- `reference/pairwise_ld/ld_reference.yaml` and, for every chromosome containing
  a genome-wide-significant variant, the format-version-2 files produced by
  `tools/resource_preparation/ld_file_preparation.sh`: forward
  `EUR_chr<CHR>.ld.gz` plus `.tbi`, reverse
  `EUR_chr<CHR>.reverse.ld.gz` plus `.tbi`, and allele-aware inventory
  `EUR_chr<CHR>.variants.tsv.gz` plus `.tbi`. The forward and reverse indexes
  let the clumper find a variant at either endpoint of PLINK's upper-triangle
  output; do not replace them with the old one-file layout.
- `reference/1000G_EUR.{bed,bim,fam}` for locus LD construction.
- One consistent unique identifier convention across VCF, pairwise-LD tables,
  and PLINK BIM. This command uses
  `{chromosome}_{position}_{reference_allele}_{alternate_allele}`; use
  `--variant-id-type rsid` instead when every relevant file uses rsIDs.

```bash
postgwas pipeline \
  --modules finemap \
  --clumping-methods standard \
  --run-config examples/fine_mapping/pipeline_susie.yaml \
  --vcf study_GRCh37.vcf.gz \
  --genome-build GRCh37 \
  --ld-folder reference/pairwise_ld \
  --population EUR \
  --variant-id-type unique \
  --finemap-method susie \
  --finemap-ld-reference reference/1000G_EUR \
  --plink plink \
  --dataset-id STUDY \
  --output-directory results/pipeline_susie
```

Validate either YAML before a long run:

```bash
postgwas config validate --config examples/fine_mapping/direct_susie.yaml
postgwas config validate --config examples/fine_mapping/pipeline_susie.yaml
```

Successful runs finish with `overall status: success`. Inspect
`results/STUDY_fine_mapping_report.html`,
`results/combined_results/final_combined_credible_sets.tsv`,
`quality_control/overlap_resolution_summary.tsv`, the engine QC summary under
`quality_control/`, and the resolved configuration under `run_metadata/`
before downstream analysis. The authoritative FLAMES handoff is
`downstream_inputs/flames/`. A converged run with no credible set is a valid
analysis outcome, not a successful FLAMES handoff. The HTML report contains
all final credible-set members and the primary-locus warnings/failures; its
small per-set bar display is visual orientation only and never filters the
complete table.
