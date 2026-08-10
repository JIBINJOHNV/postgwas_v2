# Scientific Considerations

PostGWAS coordinates methods with different scientific assumptions. A command
that completes successfully is not sufficient evidence that its inputs and
reference data were compatible.

## Genome build

Study coordinates, reference FASTA, dbSNP, LD regions, LD genotypes, gene
annotations, and tool-specific resources must describe the same genome build.
Do not relabel coordinates to change build. Use a validated liftover procedure
and retain rejected or ambiguous variants.

## Population and LD reference

LD clumping, fine-mapping, imputation, LDSC, MAGMA, PoPS, MiXeR, and related
methods rely on population-specific LD or allele-frequency information. Select
a panel that is appropriate for the analyzed sample. Record the panel name,
version, population, build, and preparation procedure.

## Allele orientation

The harmonised GWAS-VCF is the allele-aware boundary between preparation and
analysis. Effect estimates and frequencies must remain tied to the declared
effect allele. Palindromic variants require special care because strand cannot
always be resolved from alleles alone.

## Sample size and trait type

Quantitative and case-control studies have different required metadata.
Per-variant sample size, total sample size, case count, and control count are
not universally interchangeable. LDSC and other downstream methods may also
require sample and population prevalence. Use the definitions expected by the
specific method and record them with the run.

## Effect statistics and P values

Beta coefficients, odds ratios, Z scores, standard errors, raw P values,
`-log10(P)`, and `-ln(P)` are distinct representations. Declare the input type;
do not rely on magnitude-based guessing. Treat transformations and reconstructed
statistics as scientific decisions that must appear in QC and provenance.

## Filtering and missingness

MAF, INFO, P-value, MHC, indel, palindromic-variant, and allele-frequency
discordance filters can change the scientific population represented by the
data. Use the resolved YAML values, inspect removal counts by reason, and do not
describe a virtual QC subset as a newly filtered file.

## Interpreting downstream results

Clumps are LD-defined signal groupings, not necessarily causal variants.
Fine-mapping posterior probabilities depend on the locus, LD matrix, model, and
variant set. Gene and pathway results depend on gene mapping, covariates, gene
sets, LD reference, and multiple-testing procedure. Cross-trait estimates can
be unstable when heritability, overlap assumptions, or trait metadata are weak.

