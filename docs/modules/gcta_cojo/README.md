# GCTA-COJO conditional and joint analysis

The module has two distinct entry paths:

- `postgwas gcta_cojo` is direct mode. It consumes an existing GCTA `.ma` file
  and PLINK LD reference without running the formatter.
- `postgwas pipeline --modules gcta_cojo` is VCF pipeline mode. It inspects the
  PLINK BIM, selects its homogeneous `rsid` or configured
  coordinate-and-allele identifier convention, and asks the existing formatter
  to write the shared GCTA `.ma` file with that exact identifier type before
  running COJO.

For pipeline-created `.ma` files, the effect allele is VCF `ALT` (`A1`) and
`freq` is its effect-allele frequency.

## Analysis modes

- `slct` is the YAML default. It runs `--cojo-slct` when no conditioning list is
  supplied and selects independently associated signals by stepwise regression.
- `top_snps` runs `--cojo-top-SNPs N` without a p-value stopping threshold.
- `joint` requires `joint_snps` and fits that extracted set together with
  `--cojo-joint`.
- `cond` requires `condition_snps` and runs `--cojo-cond` for all included SNPs
  conditional on that list. Supplying `--condition-snps` on the CLI implies
  `cond` mode when `--cojo-mode` is omitted; any explicit mode must be `cond`.

These modes condition on variants within one GWAS. They are not mtCOJO, which
conditions one trait on genetically predicted effects of other traits, and they
are not COJO-SBLUP.

## SNP-list options

Each list is a plain-text file containing one SNP ID per line. The IDs must use
the same convention as the PLINK BIM and the `.ma` `SNP` column.

- `condition_snps` identifies known variants to adjust for. GCTA tests the other
  included variants after accounting for the listed variants.
- `joint_snps` identifies variants whose effects are estimated together in one
  multiple-SNP model. It defines the analysis set in `joint` mode.
- `cojo_extract` limits candidate variants to the listed IDs. GCTA still reads
  the full `.ma` file because unfiltered summary statistics are needed for its
  variance calculations. It cannot be combined with `joint` mode.
- `cojo_exclude` removes the listed IDs from the candidate variants.

For example, use `--condition-snps known_leads.txt` to adjust for established
lead variants while testing other eligible variants. Use
`--joint-snps model_snps.txt` when the objective is instead to estimate the
effects of a specified SNP set simultaneously.

## Scientific requirements

The GWAS and LD reference must use the same genome build, ancestry, variant IDs,
and compatible alleles. PostGWAS requires the user to declare build and reference
population, validates exact `.ma`/BIM ID overlap and allele compatibility, and
requires every conditioning, joint, or extract SNP to exist in both inputs. It
does not infer ancestry from filenames.

GCTA advises using all available GWAS summary variants because COJO uses them to
estimate phenotypic variance. Restrict the tested region with `cojo-extract` or
`cojo-chromosome` instead of pre-filtering the `.ma` input. GCTA also recommends
an LD reference drawn from the GWAS cohort or a large participating cohort and
suggests more than 4,000 reference samples; it does not recommend HapMap or 1000
Genomes panels for COJO because their sample sizes are small. PostGWAS reports a
warning, rather than silently changing the analysis, when the configured sample
size threshold is not met.

For conditional output, GCTA documents `NA` conditional estimates when the
multivariate correlation with the conditioning variants exceeds its collinearity
rule. PostGWAS preserves these rows, labels them
`not_estimable_collinearity`, and reports their count.

## Execution

Direct mode:

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Conditional mode:

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --condition-snps lead_snps.txt \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Pipeline mode starts from VCF and prepares the `.ma` input:

```bash
postgwas pipeline \
  --modules gcta_cojo \
  --vcf study.vcf.gz \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

The formatter runs once and passes its checked `.ma` artifact to the COJO step.

Every default is owned by
`src/postgwas/config/defaults/modules/gcta_cojo.yaml`. Export a reusable file
with:

```bash
postgwas config export --module gcta_cojo --style full --output gcta_cojo.yaml
```

Each run writes the resolved configuration, a checksummed completion manifest,
the normalized TSV, GCTA raw results and GCTA log, and a PostGWAS log. The final
terminal block reports mode, reference population and sample size, BIM ID type,
input/reference overlap, main findings, warnings, output paths, and any failure
reason.

## Sources

- [Official GCTA documentation and COJO option
  contract](https://yanglab.westlake.edu.cn/software/gcta/)
- [Yang et al. 2012, conditional and joint analysis](https://doi.org/10.1038/ng.2213)
