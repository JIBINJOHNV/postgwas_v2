# Optional EAF replacement source

Use the paired options below when a study supplies an internal frequency column
that may be minor-allele frequency (MAF):

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --resource-directory resources \
  --output-directory results \
  --eaf-replacement-file 'panel_{build}_{chromosome}.tsv.gz' \
  --eaf-replacement-column EUR
```

The equivalent full-run YAML is:

```yaml
modules:
  harmonisation:
    eaf_replacement_file: 'panel_{build}_{chromosome}.tsv.gz'
    eaf_replacement_column: EUR
```

Both canonical values default to `null`. Set both together; a missing partner
reports its public option and canonical YAML path. Explicit command-line values
override YAML. Export these settings through the normal configuration command:

```console
postgwas config export --module harmonisation --style full --output harmonisation.yaml
```

The file may contain all chromosomes or use `{build}` and `{chromosome}` path
placeholders. Relative CLI and run-YAML paths resolve from the current working
directory, matching existing run resource settings; sample-sheet row paths
continue to resolve from the sample-sheet directory. Other placeholders, format specifiers and filename prefixes are
invalid. The replacement value column must contain frequencies for the allele
declared by `external_eaf_mapping.effect_allele` (ALT in the packaged mapping).
The structural chromosome, position and allele columns and delimiter come from
`external_eaf_mapping`; no sample-sheet column is added.

PostGWAS first tests the selected internal frequency against `default_eaf`.
Accepted EAF stays in use. Replacement is selected only when the existing
comparison panel supports the MAF interpretation. A confirmed MAF input without
an explicit replacement stops that dataset. Insufficient or conflicting evidence
continues to stop the dataset. The replacement table supplies EAF after this
decision and never supplies evidence for its own selection or palindrome
orientation. An ordinary sample-sheet external EAF source keeps its existing
behavior.

Messages distinguish three outcomes: study EAF retained; input MAF supported
and the supplied replacement selected; or unresolved evidence and the dataset
stopped. Supported MAF without a replacement reports both required options and
their YAML keys. It stops that dataset before chromosome analysis; other
datasets continue under the existing batch failure policy.

The selected frequencies are a population-reference proxy, not recovered
study-measured EAF. Choose the intended population and source genome build.
The [GWAS-SSF field definitions](https://github.com/gwas-catalog/summary-statistics-standard/blob/main/gwas-ssf_v1.0.2.tex)
define frequency relative to the listed effect allele; the effect allele need
not be minor. [1000 Genomes](https://www.internationalgenome.org/faq/how-do-i-find-out-information-about-a-single-variant/)
reports ALT frequencies, including separate values for multiple ALT alleles.
Diagnostic folding and a swapped `1-AF` assignment require biallelic reference
sites. Known extra alleles prevent that inference; complete site membership in
the supplied reference remains a declared input assumption.

A supplied replacement specification must match at least one existing nonempty
file, even if a particular dataset later retains its internal EAF. Exact active
build/chromosome coverage and table contents are validated during dataset
preparation. A shared whole-genome replacement uses the same bounded staging
and worker memory accounting as existing external EAF tables.

Before a completed direct command can resume, its replacement template is
expanded and the resulting exact files join the durable checkpoint input
fingerprints. This applies to files outside `resources.root`, YAML settings,
`--option VALUE` and `--option=VALUE`. A changed tracked source invalidates the
completed checkpoint under the canonical restart policy. Original study MAF,
replacement source and final EAF interpretation remain distinct in provenance.
