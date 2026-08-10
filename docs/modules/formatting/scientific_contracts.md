# Formatter scientific contracts

Reviewed against the upstream documentation and the PostGWAS consumer code on
2026-08-06. “Required” means consumed by the named external tool. A PostGWAS
extension is retained only when a later PostGWAS stage needs it.

| Target | Required external input | Frequency requirement | Sample-size requirement | PostGWAS decision |
|---|---|---|---|---|
| MAGMA | SNP-location: headerless `SNP CHR BP`; p-values: headed `SNP P` plus per-variant `N_COL` | None | Total sample size, including case-control results | Pipeline mode scans BIM field 2 and selects the MAGMA target's `rsid` or configured `unique` convention. Formatter retains `REF ALT` for optional exact-ID/coordinate/allele intersection; the MAGMA runner writes only the required headerless first three columns. `N_COL <- FORMAT/SS`; never use `NEF` here |
| GCTA fastBAT | `SNP P` | None | None | Bound `10^-LP` using the configured minimum raw p-value |
| GCTA mBAT-combo | `SNP A1 A2 freq BETA SE P N` | Effect-allele frequency | Total sample size | `A1 <- ALT`, `A2 <- REF`, `freq <- AF`, and `N <- SS` |
| FINEMAP | `rsid chromosome position allele1 allele2 maf beta se` | Minor allele frequency | `n_samples` in the master file | `maf <- min(EAF, 1-EAF)`; retain `NEF` so the PostGWAS adapter can create locus-specific `n_samples` |
| SuSiE-RSS | Z scores, an allele-aligned LD matrix, and the configured PostGWAS identity/coordinate fields | Optional in `susie_rss`; used only when a MAF threshold is requested | `n` is recommended | The current PostGWAS engine uses `EZ`, LD, and median `NEF`; it does not pass MAF, so `AF` is not exported |
| PRED-LD | `snp chr pos A1 A2 beta SE`; `A1=ALT`, `A2=REF` | Study AF is not an input; `--maf` filters the reference resource | None | `NC SS AF LP SI` remain carry-through metadata because PostGWAS re-harmonises the imputed results |
| CBIIT LDSC | `SNP A1 A2 P`, a signed statistic, and sample size | `FRQ` is optional QC input; `munge_sumstats.py` converts allele frequency to MAF for filtering | Binary: `N_CAS` and `N_CON`; quantitative: `N` | Infer binary when any `FORMAT/NC` value exists; infer quantitative only when `NC` is entirely missing and `FORMAT/NCO` is present. For quantitative traits, write `N <- NCO` |
| MiXeR `fit1`/`test1` | `SNP CHR BP A1 A2 N Z` | Not required | Binary: effective `N = 4/(1/Ncase + 1/Ncontrol)`; quantitative: total `N` | Write harmonised `ALT` as effect allele `A1`, `REF` as other allele `A2`, and `FORMAT/NEF` as `N`; apply configured INFO, sample-size, and SNP checks before export |

Allele direction is invariant across every effect-based output: GWAS-VCF
`ALT` is the effect allele, so it becomes FINEMAP `allele1`, PRED-LD `A1`, and
LDSC `A1`. GWAS-VCF `REF` becomes the corresponding other/non-effect allele.

## Configuration source of truth

All formatter extraction and output schemas are declared in
`src/postgwas/config/defaults/modules/formatting.yaml` and validated by the
typed configuration model. The mappings have two explicit stages:

1. `vcf_fields`: canonical formatter-table column → bcftools query expression.
2. `variant_identifiers`: general default and per-target ID selection, extraction pattern, and unique-ID template.
3. `exports.<target>`: canonical formatter-table column → downstream tool column.

The same YAML also owns output filename patterns, named transformations,
numeric parsing, validation columns, chromosome-label normalization,
variant-ID selection, study-design count columns, chromosome
selection, pipeline-to-format dependencies (`module_formats`), stable execution
order (`format_order`), input null markers, table delimiter, output null value,
temporary-table and atomic-output naming, I/O buffer size, and MiXeR QC
parameters. Canonical
table names are not fixed in Python: `vcf_fields`, semantic roles, validation
lists, and export mappings are validated together. Exporter code applies these
declarations; it does not maintain a second set of column aliases.

## Corrected formatter mappings

| Target | Previous mapping | Reviewed mapping |
|---|---|---|
| MAGMA `N_COL` | `FORMAT/NEF` | `FORMAT/SS` (total sample size) |
| FINEMAP `maf` | `FORMAT/AF` was written directly and mislabeled | `min(FORMAT/AF, 1-FORMAT/AF)` |
| SuSiE `AF` | Exported but not read | Removed from the configured engine input |
| PRED-LD `NC` | `FORMAT/NC` (cases), although the PostGWAS handoff treated it as controls | `FORMAT/NCO` (controls) |
| LDSC sample size | Case/control columns when available | Binary: `N_CAS <- FORMAT/NC`, `N_CON <- FORMAT/NCO`; quantitative: `N <- FORMAT/NCO` |

Study type is inferred only from the sample-count values; no separate VCF
metadata declaration is expected or written. `NCO` must contain at least one
value. If `NC` is absent or entirely missing, the trait is quantitative; if any
`NC` value is present, it is binary.

## Primary sources

- [Official MAGMA software, documentation, and reference resources](https://cncr.nl/research/magma/)
- [Official PLINK BIM format](https://www.cog-genomics.org/plink/1.9/formats#bim)
- [VCF 4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [GCTA fastBAT and mBAT-combo documentation](https://yanglab.westlake.edu.cn/software/gcta/)
- [FINEMAP documentation](https://christianbenner.com/)
- [susieR `susie_rss` documentation](https://stephenslab.github.io/susieR/reference/susie_rss.html)
- [PRED-LD command-line documentation](https://github.com/pbagos/PRED-LD)
- [CBIIT LDSC `munge_sumstats.py`](https://github.com/CBIIT/ldsc/blob/ldsc39/munge_sumstats.py)
- [Official MiXeR GWAS input and univariate workflow](https://github.com/precimed/mixer#gwas-summary-statistics-format)
