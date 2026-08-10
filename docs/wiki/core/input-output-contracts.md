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
- [Logging and Reproducibility](logging-and-reproducibility.md)
