# Formatter scientific contracts

Reviewed against the upstream documentation and the PostGWAS consumer code on
2026-08-13. “Required” means consumed by the named external tool. A PostGWAS
extension is retained only when a later PostGWAS stage needs it.

| Target | Required external input | Frequency requirement | Sample-size requirement | PostGWAS decision |
|---|---|---|---|---|
| MAGMA | SNP-location: headerless `SNP CHR BP`; p-values: headed `SNP P` plus per-variant `N_COL` | None | Total sample size, including case-control results | Pipeline mode scans BIM field 2 and selects the MAGMA target's `rsid` or configured `unique` convention. The shared duplicate policy defaults to excluding every conflicting selected-ID group; P-value ranking requires an explicit override. Pipeline MAGMA consumes these exact paired formatter artifacts. Formatter retains `REF ALT` for optional exact-ID/coordinate/allele intersection; the MAGMA runner writes only the required headerless first three columns. `N_COL <- FORMAT/SS`; never use `NEF` here |
| GCTA fastBAT | `SNP P` | None | None | Bound `10^-LP` using the configured minimum raw p-value |
| GCTA COJO / mBAT-combo | `SNP A1 A2 freq BETA SE P N` | Effect-allele frequency | Total sample size | `A1 <- ALT`, `A2 <- REF`, `freq <- AF`, and `N <- SS`. The `.ma` schema is shared; pipeline reports identify the resolved downstream COJO mode and configured fixed SNP count without claiming that formatting differs by analysis mode |
| FINEMAP | `rsid chromosome position allele1 allele2 maf beta se` | Minor allele frequency | `n_samples` in the master file | `maf <- min(EAF, 1-EAF)`; retain `NEF` so the PostGWAS adapter can create locus-specific `n_samples` |
| SuSiE-RSS | Z scores, an allele-aligned LD matrix, and the configured PostGWAS identity/coordinate fields | Optional in `susie_rss`; used only when a MAF threshold is requested | `n` is recommended | The current PostGWAS engine uses `EZ`, LD, and median `NEF`; it does not pass MAF, so `AF` is not exported |
| PRED-LD | `snp chr pos A1 A2 beta SE`; `A1=ALT`, `A2=REF` | Study AF is not an input; `--maf` filters the reference resource | None | `NC SS AF LP SI` remain carry-through metadata because PostGWAS re-harmonises the imputed results. Only canonical chromosomes in the configured `chromosomes` list are published; excluded labels and row counts are recorded explicitly |
| CBIIT LDSC | `SNP A1 A2 P`, a signed statistic, and sample size | `FRQ` is optional QC input; `munge_sumstats.py` converts allele frequency to MAF for filtering | Binary: `N_CAS` and `N_CON`; quantitative: `N` | Infer binary when any `FORMAT/NC` value exists; infer quantitative only when `NC` is entirely missing and `FORMAT/NCO` is present. For quantitative traits, write `N <- NCO`. For binary traits, calculate each valid variant's `NC/(NC+NCO)` and reduce it with configured `median` (default) or `mean`; never sum repeated variant-level counts. When `--merge-alleles` is supplied, select by rsID plus LDSC-compatible allele order/strand before duplicate-ID validation; the heritability pipeline reuses the same file during munging |
| MiXeR `fit1`/`test1` | `SNP CHR BP A1 A2 N Z` | Not required | Binary: effective `N = 4/(1/Ncase + 1/Ncontrol)`; quantitative: total `N` | Write harmonised `ALT` as effect allele `A1`, `REF` as other allele `A2`, and `FORMAT/NEF` as `N`; apply configured INFO, sample-size, and SNP checks before export |

Allele direction is invariant across every effect-based output: GWAS-VCF
`ALT` is the effect allele, so it becomes FINEMAP `allele1`, PRED-LD `A1`, and
LDSC `A1`. GWAS-VCF `REF` becomes the corresponding other/non-effect allele.

## Configuration source of truth

All formatter extraction and output schemas are declared in
`src/postgwas/config/defaults/modules/formatting.yaml` and validated by the
typed configuration model. The mappings and reporting semantics have these
explicit configuration layers:

1. `input_contract`: accepted PostGWAS provenance and genome-build metadata,
   supported builds, and canonical chromosome/allele patterns. These values are
   validated; the formatter never repairs or re-harmonises VCF labels or alleles.
2. `vcf_fields`: canonical formatter-table column → bcftools query expression.
3. `canonical_columns`: semantic roles, including the canonical `-log10(P)`,
   effect-allele-frequency, and four sample-size sources used to report exact
   statistical representations and saved headers.
4. `sample_size_reporting`: trait-aware meanings and target-use notes for the
   total, effective, case, and control/quantitative-total sample-size roles.
5. `variant_identifiers`: general default and per-target ID type and duplicate
   policy, extraction pattern, and unique-ID template.
6. `ldsc_reference`: optional merge-alleles path, input columns, delimiter
   detection, and table-reading constraints for LDSC preselection.
7. `ldsc_sample_prevalence`: the schema-limited `median` or `mean` reduction of
   valid per-variant case fractions for binary-trait pipeline handoff.
8. `exports.<target>`: canonical formatter-table column → downstream tool column.
9. `custom_output.field_contracts`: canonical source, transformation, and
   validity rule used when direct CLI custom columns are requested.

The same YAML also owns output filename patterns, named transformations,
numeric parsing, validation columns, external-reference chromosome-label
normalization,
variant-ID selection, study-design count columns, chromosome
selection, pipeline-to-format dependencies (`module_formats`), stable execution
order (`format_order`), input null markers, table delimiter, output null value,
temporary-table and atomic-output naming, I/O buffer size, and MiXeR QC
parameters. Canonical
table names are not fixed in Python: `vcf_fields`, semantic roles, validation
lists, and export mappings are validated together. Exporter code applies these
declarations; it does not maintain a second set of column aliases.
The runtime schema report reads these same resolved mappings and transformations
to log every source-to-saved column name and identify raw P, `-log10(P)`, EAF,
derived MAF, and each sample-size source, numerical meaning, trait-specific
formula, saved header, and downstream use without maintaining a parallel Python
scientific-description registry. These reporting descriptions are retained in
resolved provenance but excluded from the scientific-output checkpoint digest,
so a wording-only change does not regenerate unchanged formatter artifacts.

## LDSC sample-prevalence reduction

Official LDSC accepts one binary-trait sample prevalence through `--samp-prev`
and requires it together with population prevalence for liability-scale
conversion ([LDSC command](https://github.com/bulik/ldsc/blob/master/ldsc.py),
[heritability documentation](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation)).
It does not prescribe how to reduce variant-level case/control counts. PostGWAS
therefore calculates `NC/(NC+NCO)` independently for every valid written LDSC
variant and applies the schema-validated `median` or `mean`. `median` is the
robust default. `sum` is deliberately unsupported because VCF variant rows
repeatedly describe the study sample rather than independent participant
groups. The reduction, row count, range, and result are logged and persisted.

## Custom CLI output contract

The optional custom table is additive and does not change a built-in tool
schema or downstream handoff. Users activate it with `--custom-output`, name
the identifier column with `--id`, and optionally request other fields by
supplying their output header names. The CLI option order is the output column
order. The identifier uses the existing configured `rsid` or `unique` policy.

All requested sources are required per written row. Missing or non-finite
values, nonpositive position/sample-size/standard-error values, negative LP,
and EAF/INFO outside `[0,1]` are excluded and counted. `ALT` remains the effect
allele and `REF` the other allele. `--p` applies the same configured
`10^-LP` conversion and minimum numeric bound as built-in outputs; `--lp`
retains LP. `--maf` applies `min(EAF, 1-EAF)` while `--eaf` retains effect-allele
frequency. Both members of either pair may be requested because each output
expression is represented independently even though it shares a canonical
source. The custom table does not infer study design and does not claim
compatibility with an unspecified external tool.

## Shared duplicate-identifier contract

Every built-in format uses `variant_identifiers.default_duplicate_policy`, with
an optional value in `target_duplicate_policies`. The packaged default is
`exclude_all`. Exact repeated extracted records are collapsed before conflicting
groups are classified. Supported formatter-side reference matching (currently
LDSC merge alleles) runs before the fallback policy. `error` stops;
`exclude_all` removes the complete group;
`most_significant`, `highest_maf`, and `highest_info` retain a record only when
one valid rank is strictly greatest. Tied maxima and missing rank values remove
the complete group rather than falling back to row order. The configured
absolute `duplicate_rank_tolerance` prevents floating-point representations of
equivalent ranks from being treated as a strict winner.

P-value, MAF, and INFO ranks do not prove biological identity. P-value ranking
can additionally create ascertainment bias and winner's-curse effects, so none
of the ranked policies is the default. Missing-ID, exact-repeat, conflicting
group, ranked-winner, unresolved-group, duplicate-row, and required-value counts
are recorded separately in the terminal summary, canonical log, result metadata,
and completion manifest. Multi-table targets such as MAGMA receive one resolved
frame, preserving identical unique-ID sets and row order.

Targets with the same resolved `(identifier type, duplicate policy)` reuse one
validated selection from the common canonical frame. The cache stores only
selections needed by multiple targets and releases each frame after its final
consumer. LDSC selection with `--merge-alleles` remains outside this cache
because its reference join changes the candidate set. Before any exporter is
called, the formatter independently verifies that the final candidate contains
no missing, empty, or duplicated selected IDs. The validation result and whether
the selection was reused are recorded for every target; a failed invariant stops
before a scientific artifact is written.

In pipeline mode, `run_magma_runner` passes these exact filtered formatter
artifacts to the MAGMA service. The service still validates identifier,
coordinate, and allele consistency and invokes MAGMA with `duplicate=error` as
a downstream safety boundary. Its configured `lowest_p` consolidation remains
relevant only when a user runs the MAGMA module directly with independently
prepared tables that repeat an otherwise coordinate/allele-consistent variant;
it is normally inactive for formatter-produced inputs.

## LDSC HapMap3 selection contract

The direct formatter does not require a HapMap3 list. If `--merge-alleles` is
omitted, the shared duplicate policy is applied directly. If the option is
supplied, the configured `SNP`, `A1`, and `A2` reference columns are validated
first: rsIDs must be unique, and alleles must be single-base, non-palindromic
A/C/G/T pairs accepted by the official LDSC `munge_sumstats.py` implementation.

An input record is compatible when its effect/other alleles equal the reference
pair in either order, directly or after strand complementation. Records absent
from the reference and records with incompatible alleles are excluded and
counted separately. A duplicated input rsID is reference-resolved when exactly
one record is compatible. More than one compatible record remains ambiguous and
is passed to the shared configured fallback; it is never selected by input
order. The `heritability` pipeline requires `--merge-alleles`, supplies it
to this formatter selection, and then supplies the identical file to
`munge_sumstats.py`.

This policy follows the CBIIT LDSC `VALID_SNPS` and `MATCH_ALLELES` definitions,
which accept non-palindromic allele pairs under allele-order and strand changes.
The formatter completion manifest fingerprints the reference because its
contents determine the selected output variants.

## Corrected formatter mappings

| Target | Previous mapping | Reviewed mapping |
|---|---|---|
| MAGMA `N_COL` | `FORMAT/NEF` | `FORMAT/SS` (total sample size) |
| FINEMAP `maf` | `FORMAT/AF` was written directly and mislabeled | `min(FORMAT/AF, 1-FORMAT/AF)` |
| SuSiE `AF` | Exported but not read | Removed from the configured engine input |
| PRED-LD `NC` | `FORMAT/NC` (cases), although the PostGWAS handoff treated it as controls | `FORMAT/NCO` (binary controls; quantitative total N), explicitly labelled as carry-through metadata |
| LDSC sample size | Case/control columns when available | Binary: `N_CAS <- FORMAT/NC`, `N_CON <- FORMAT/NCO`; quantitative: `N <- FORMAT/NCO` |

Study type is inferred only when LDSC or MiXeR is selected, because those are
the only formatter contracts with trait-specific sample-size interpretation.
MAGMA, GCTA, SuSiE, FINEMAP, and PRED-LD validate their own configured fields
without triggering inference. For LDSC or MiXeR, inference uses only sample-count
values; no separate VCF metadata declaration is expected or written. `NCO` must
contain at least one value. If `NC` is absent or entirely missing, the trait is
quantitative; if any `NC` value is present, it is binary.

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
- [Winner's-curse ascertainment bias in GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC3500533/)
- [PLINK 2 duplicate-ID handling](https://www.cog-genomics.org/plink/2.0/filter)
