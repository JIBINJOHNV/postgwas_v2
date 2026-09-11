# Input and output contracts

Files passed between PostGWAS stages must retain their scientific context.
Correct downstream interpretation depends on genome build, ancestry, dataset
identity, schema, and effect-allele convention as well as the filename.

## Input boundary

Harmonisation accepts raw GWAS summary statistics together with an explicit
sample sheet, run configuration, resources, and output location. Column roles,
study design, sample-size information, and permitted transformations must be
declared and validated rather than inferred silently.

Downstream commands commonly require:

- a harmonised GWAS-VCF;
- a stable dataset identifier;
- an output directory;
- a resolved run configuration;
- module-specific reference data compatible with the input genome build and
  population requirements.

Use standalone command help to identify inputs that must be supplied directly.
Pipeline help hides artifacts created by preceding steps and continues to show
external resources that the workflow cannot create.

### PostGWAS-origin check for study VCFs

Public commands accepting a study GWAS-VCF require PostGWAS harmonisation
metadata before analysing its records. This applies to filtering, QC, LD-block
annotation, LD clumping, Manhattan plots, formatting, standalone concordance
(`postgwas --validate`), and every pipeline, including a single-target pipeline
for MAGMA, PoPS or FLAMES.

The shared check requires each of these headers exactly once, with a non-empty
text value:

| Header (default name) | Meaning |
|---|---|
| `postgwas_version` | Version that produced the VCF; it need not equal the currently installed version |
| `postgwas_dataset_id` | Dataset identity recorded by harmonisation |
| `postgwas_vcf_status` | VCF creation/structural-validation status recorded by harmonisation |

Names are resolved from the existing canonical
`modules.formatting.input_contract.provenance_headers` configuration. No
separate per-module origin policy is introduced. The existing genome-build,
sample, index and required-field checks still apply wherever required by the
consumer; the origin check does not add formatter-specific fields to other
direct commands.

Missing, empty, duplicate or malformed origin headers stop the command with
the affected filename and an instruction to rerun harmonisation. Older VCFs
without this metadata must be regenerated from the original summary statistics;
do not add headers manually to bypass the check. A different output dataset
label does not rewrite the recorded source identity.

The check reuses the header already read for validation. It neither scans nor
changes variant values. External reference VCFs (such as dbSNP and frequency
panels), raw harmonisation inputs, and direct-module table inputs are exempt.
Internal QC assessment utilities remain usable for their existing workflows.

These are **declared provenance markers, not a cryptographic signature**. They
do not establish authenticity or prove that every variant passed scientific
QC; record-level validation and interpretation remain necessary.

## Information that must remain consistent

Check the following information whenever a file is passed to another module:

- file type and path;
- genome build and chromosome convention;
- population or ancestry context;
- sample or dataset identifier;
- schema or format version;
- allele and effect-statistic conventions;
- the input and settings used to create it.

The absence of an automatic compatibility error does not prove that two files
are scientifically compatible.

## Pipeline validation evidence

The pipeline records input and resource checks in a common, atomic YAML audit.
Each record identifies the file or cross-file requirement, its consuming
modules, the checks actually performed, observed metrics, and its outcome.
An availability pass does not establish correct contents; a header/index pass
does not establish valid values in every variant record.

External resources are checked at startup under the selected modules' current
preflight contracts. Generated tables and modified VCFs must still pass the
checks required by their consumers after creation. Startup deferrals are retained
as such, rather than marked passed merely because another step succeeded.

For PLINK references, the BIM identifier convention can guide formatting, but
the subsequently generated table must still match the reference. BIM allele
order does not establish reference-genome REF/ALT or effect orientation.
See [Pipeline Input Validation](pipeline-input-validation.md) for coverage,
configuration, evidence reuse, and remaining limitations.

## Output boundary

Treat only validated final files as analysis results. Temporary or incomplete
work must remain separate and every reused output should be traceable to its
inputs and configuration.

Before reusing an output, verify:

1. the producing stage completed successfully;
2. expected files exist and pass their format validation;
3. genome build and population assumptions match the consumer;
4. variant identifiers and allele conventions are compatible;
5. sample-size fields have the interpretation required by the method;
6. the output belongs to the intended dataset and run.

## Standalone and pipeline consistency

For standalone commands, provide every required input explicitly. In pipeline
mode, PostGWAS supplies files created by preceding stages; external references
and study-specific inputs must still be configured by the user.

## Related pages

- [Configuration](configuration.md)
- [Pipeline Workflow](pipeline-workflow.md)
- [Pipeline Input Validation](pipeline-input-validation.md)
- [Logging and Reproducibility](logging-and-reproducibility.md)
