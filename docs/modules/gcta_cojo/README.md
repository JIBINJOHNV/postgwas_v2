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

In every `.ma` file, `A1` is the effect allele, `A2` is the other allele,
`freq` is the `A1` frequency, and `b` is the regression coefficient. For a
case-control GWAS, `b` must be `log(OR)`, not the odds ratio itself. Direct mode
cannot infer whether `A1` is REF or ALT and records it only as the input effect
allele.

## Analysis modes

- `slct` is the YAML default. It runs `--cojo-slct` when no conditioning list is
  supplied and selects independently associated signals by stepwise regression.
- `top_snps` runs `--cojo-top-SNPs N` without a p-value stopping threshold.
  Supplying `--cojo-top-snps` implies this mode when `--cojo-mode` is omitted.
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
  variance calculations. In conditional mode it must include every
  conditioning SNP. It cannot be combined with `joint` mode.
- `cojo_exclude` removes the listed IDs from the candidate variants.

A conditioning or joint-model SNP must not also occur in `cojo_exclude`.
PostGWAS verifies the SNPs written by GCTA after fitting: the `.given.cojo` SNP
set must exactly equal `condition_snps`, and the joint `.jma.cojo` SNP set must
exactly equal `joint_snps`. This catches silent loss caused by MAF, frequency,
chromosome, extract/exclude, allele, or LD-reference filters.

For example, use `--condition-snps known_leads.txt` to adjust for established
lead variants while testing other eligible variants. Use
`--joint-snps model_snps.txt` when the objective is instead to estimate the
effects of a specified SNP set simultaneously.

## Input requirements

The GWAS and LD reference must use the same genome build, ancestry, variant IDs,
and compatible alleles. PostGWAS requires the user to declare build and reference
population and does not infer ancestry from filenames.

Reference compatibility is assessed across all `.ma` rows before GCTA starts. A
variant is **GCTA-usable** only when its SNP ID occurs in the BIM and its complete
unordered allele pair matches the BIM allele pair. The YAML default requires at
least 70% usable variants (`input_validation.minimum_reference_overlap_fraction:
0.70`); a lower fraction stops the run with separate counts for absent IDs and
shared IDs with incompatible alleles. This threshold is a PostGWAS validation
policy, not a native GCTA option.

When the usable fraction is at least the configured threshold, absent BIM IDs
are reported and cannot enter COJO. Shared IDs with incompatible allele pairs
are also reported, written to a checksummed exclusion artifact, merged with any
user `cojo_exclude` list, and supplied through one GCTA `--exclude` argument.
The original complete `.ma` file remains unchanged because GCTA uses the full
summary-statistic set to estimate phenotypic variance. Every explicitly
requested conditioning, joint, or extract SNP is stricter: it must be
GCTA-usable even when the overall usable fraction passes.

PostGWAS does not strand-complement or relabel alleles during this check. The
pipeline contract is a harmonised GWAS-VCF and a build-matched, standardised
reference; widespread strand disagreement therefore indicates incompatible
inputs and should fail the usable-overlap threshold rather than be guessed away.
GCTA's upstream COJO reader checks the summary-statistic A1 against the genotype
alleles and skips A1 failures, but it does not validate the complete A1/A2 set.
The explicit PostGWAS full-pair check and exclusion therefore prevents a
partial allele match from being accepted ambiguously.

Direct mode accepts either GCTA's official `.ma` header
`SNP A1 A2 freq b se p N` or the canonical formatter header recorded in the
resolved configuration. It rejects non-finite statistics, frequencies outside
`(0, 1)`, non-positive standard errors, p-values outside `[0, 1]`, and sample
sizes below the configured threshold before starting the COJO analysis. That
threshold defaults to GCTA's minimum of 10 and cannot be configured below 10.
Per-SNP effective sample sizes may be non-integer because GCTA reads `N` as a
numeric value.

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

`--cojo-gc` is opt-in; omitting it leaves genomic control disabled. When it is
enabled, PostGWAS validates and preserves GCTA's native `p_GC`, `pJ_GC`, and
`pC_GC` column names instead of relabelling them. GCTA's supported parameter
ranges are schema-validated, including `--cojo-p` `(0, 0.05]`,
`--cojo-top-snps` `1..10000`, `--cojo-window-kb` `1..100000`,
`--cojo-collinear` `[0.01, 0.99]`, an explicit genomic-control lambda
`[1, 10]`, `--cojo-diff-freq` `[0, 1]`, `--cojo-maf` `[0, 0.5]`, and
numeric chromosome codes `1..100`.

Genomic inflation can reflect polygenicity as well as confounding, so
PostGWAS does not treat genomic control as routinely required. Enable it only
when it is justified for the supplied GWAS. The generic `--seed` is retained
for the common PostGWAS interface, but native COJO is deterministic for fixed
inputs, parameters, LD reference, and GCTA version and does not consume it.

## Execution

### Automatic chromosome execution for stepwise selection

Genome-wide `slct` runs use chromosome execution by default. This applies only
when `analysis.chromosome` is unset and the GCTA-usable `.ma`/BIM overlap
contains more than one unique positive numeric chromosome. `top_snps`, `joint`,
`cond`, an explicitly selected chromosome, and non-numeric or mixed chromosome
labels retain the single-command path. The safe fallback is deliberate: it
prevents an unrecognised chromosome label from being silently omitted.

Every chromosome worker receives the complete `.ma` file and the same PLINK
reference, exclusions, thresholds and scientific settings; only GCTA's `--chr`
candidate restriction differs. This follows GCTA's recommendation to retain
all available summary variants for phenotypic-variance estimation while doing
COJO selection chromosome by chromosome. Results are merged in numeric
chromosome order. Per-chromosome LD matrices are combined as a block-diagonal
matrix because variants on different chromosomes do not share physical LD.

PostGWAS starts the chromosome with the largest GCTA-usable variant count alone
as a memory pilot and samples the GCTA process plus descendant processes when
the operating system exposes them. The planned memory per worker is the larger
of the configured floor and the measured peak multiplied by the configured
safety factor. The log distinguishes complete process-tree samples from
root-process-only samples. Worker count is then the smallest limit imposed by:

- the number of remaining chromosomes;
- `execution.threads / chromosome_execution.threads_per_worker`;
- available `execution.memory_gb` after subtracting the PostGWAS parent process;
  and
- `chromosome_execution.max_workers`, unless it is `auto`.

If neither the GCTA root process nor its tree can be measured, execution
continues conservatively with one chromosome worker at a time. A failed
chromosome is retried without rerunning successful chromosomes. Retries are bounded by
`failed_chromosome_retries` and `retry_max_workers`; an exhausted chromosome
fails the run, and staged files are not published as completed outputs.

The chromosome `.jma.cojo` files contain the selected variants and GCTA's joint
`bJ`, `bJ_se`, and `pJ` estimates. Their `.ldr.cojo` files are merged as a
block-diagonal matrix. Standalone `postgwas gcta_cojo` then runs one
genome-wide `--cojo-cond` analysis using the merged selected-SNP list so its
public output contract also includes the separate `.cma.cojo` conditional
table, including tested variants on chromosomes where no signal was selected.

An internal consumer may request the narrower `selection_only` output contract
when it uses only the stepwise-selection results. LD clumping does this by
default because it reads the normalized `.jma.cojo` to obtain selected signals,
joint statistics, and physical loci; it does not consume `.cma.cojo`. Skipping
that separate `--cojo-cond` analysis cannot alter `.jma.cojo`, the selected
signal set, or downstream locus grouping. The contract, reconstruction status,
commands, memory measurements, worker limits, attempts, and retry decisions are
recorded in the canonical log and completion manifest, so a selection-only run
cannot be mistaken for a complete conditional-output run.

The canonical controls and shipped defaults are:

```yaml
chromosome_execution:
  enabled: true
  max_workers: auto
  threads_per_worker: 1
  minimum_memory_gb_per_worker: 1.0
  memory_safety_factor: 1.25
  memory_poll_interval_seconds: 0.1
  failed_chromosome_retries: 1
  retry_max_workers: 1
```

`enabled: false` restores one genome-wide GCTA command. Increasing
`threads_per_worker` can reduce the number of simultaneous chromosomes because
the scheduler always respects the resolved total thread budget. Raising the
memory floor or safety factor is conservative; lowering either allows more
workers but increases the risk of memory pressure. When the complete output
contract is requested, the final genome-wide conditional reconstruction is
deliberately sequential.

## Direct mode

Each example consumes the complete `.ma` input and a matching PLINK
BED/BIM/FAM prefix. Separate output directories keep the four analyses distinct.

### Stepwise selection

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_slct
```

### Fixed-count selection

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --cojo-mode top_snps \
  --cojo-top-snps 10 \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_top_snps
```

Here `10` is the requested maximum number of selected signals, not a
significance threshold; choose a count justified by the analysis.

### Joint model for a specified SNP list

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --cojo-mode joint \
  --joint-snps model_snps.txt \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_joint
```

### Conditional association

```bash
postgwas gcta_cojo \
  --cojo-file study.ma \
  --cojo-mode cond \
  --condition-snps lead_snps.txt \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_cond
```

## Pipeline mode

Supply an indexed, single-sample PostGWAS-harmonised VCF, the compatible PLINK
reference, and any analysis-specific SNP list. The pipeline inspects BIM IDs
and prepares the `.ma` input; do not pass `--cojo-file` in this mode.

### Stepwise selection

```bash
postgwas pipeline \
  --modules gcta_cojo \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_slct_pipeline
```

### Fixed-count selection

```bash
postgwas pipeline \
  --modules gcta_cojo \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --cojo-mode top_snps \
  --cojo-top-snps 10 \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_top_snps_pipeline
```

### Joint model for a specified SNP list

```bash
postgwas pipeline \
  --modules gcta_cojo \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --cojo-mode joint \
  --joint-snps model_snps.txt \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_joint_pipeline
```

### Conditional association

```bash
postgwas pipeline \
  --modules gcta_cojo \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --cojo-mode cond \
  --condition-snps lead_snps.txt \
  --cojo-reference-prefix reference/EUR_ld \
  --genome-build GRCh37 \
  --cojo-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_cond_pipeline
```

The formatter runs once and passes its checked `.ma` artifact to the COJO step.
Because the `.ma` schema is shared across GCTA consumers and COJO modes, its
columns do not change by analysis. In pipeline mode, the progress title,
formatter table, terminal summary, canonical log, and HTML report nevertheless
identify the resolved downstream analysis—for example, `GCTA-COJO top-SNP
selection (10 SNPs requested)`—instead of displaying only the generic GCTA
schema name.

## Configuration

Every default is owned by
`src/postgwas/config/defaults/modules/gcta_cojo.yaml`. Export a reusable file
with:

```bash
postgwas config export --module gcta_cojo --style full --output gcta_cojo.yaml
```

## Outputs

Each run writes the resolved configuration, a checksummed completion manifest,
the normalized TSV, GCTA raw results and GCTA log, and a PostGWAS log. Successful
standalone stepwise and fixed-count runs retain `.jma.cojo`, `.ldr.cojo`, and
`.cma.cojo`; an internal selection-only stepwise run intentionally omits the
separate `.cma.cojo`; joint runs retain `.jma.cojo` and `.ldr.cojo`; conditional runs retain
`.cma.cojo` and `.given.cojo`. When allele-incompatible shared IDs are found, the
run also retains the effective `.allele_mismatch.exclude.txt` file passed to
GCTA; this file is the deterministic union of those IDs and any user exclusion
list. A stepwise or fixed-count run for which GCTA
explicitly reports that no SNP was selected is a valid zero-row result, not a
missing-output failure. Its normalized table, log, warning, metrics, and
completion manifest remain resumable. Completion identity includes the GCTA
version and checksummed executable plus the PostGWAS, Python, and platform
identity.

The final terminal block reports mode, reference population and sample size,
BIM ID type, input/reference overlap, main findings, warnings, output paths, and
any failure reason.

## Sources

- [Official GCTA documentation and COJO option
  contract](https://yanglab.westlake.edu.cn/software/gcta/)
- [GCTA upstream option parsing and validation](https://github.com/JianYang-Lab/GCTA/blob/main/main/option.cpp)
- [GCTA upstream COJO output implementation](https://github.com/JianYang-Lab/GCTA/blob/main/main/joint_meta.cpp)
- [Yang et al. 2012, conditional and joint analysis](https://doi.org/10.1038/ng.2213)
