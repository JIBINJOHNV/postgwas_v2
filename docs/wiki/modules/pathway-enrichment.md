# Pathway Enrichment

## Purpose

The standalone command queries multiple online pathway, drug-gene, and
interaction providers for one list of human gene symbols.

## What the analysis does

It reads unique symbols from the first column, dispatches automatically to the
managed isolated pathway-enrichment runtime, and independently attempts OmniPath,
DSigDB/WebGestaltR, BioGRID, DAVID, DGIdb, Enrichr, g:Profiler, ToppGene, and
STRINGdb.

## When to use it

Use it for exploratory follow-up of a prespecified gene list. Define the list
and background independently of results and acknowledge selection/multiplicity.

## Input requirements

A delimited file with a header and gene symbols in the first column; BioGRID
API key; DAVID-registered email; optional DSigDB GMT; network access; the managed
runtime's Python/R dependencies; dataset ID and output directory.

The [all-tools installer](../getting-started/installation.md) installs this
runtime below the main environment. Activate the main PostGWAS environment and
use the ordinary command; do not manually activate an environment named
`enricher`. The launcher locates
`share/postgwas/environments/enrichment/bin/python` under the active environment,
or uses an explicitly configured `POSTGWAS_ENRICHMENT_PYTHON` override. A missing
or unusable runtime fails before provider requests begin.

## Command

This is a standalone command, not a selectable pipeline target. It consumes a
gene-symbol list, not a GWAS-VCF or the complete output table of another module.
Select and document the intended genes before exporting the one-column input.

```console
postgwas pathway_enrichment --gene-input-file PATH [options]
```

## Minimal example

```console
postgwas pathway_enrichment \
  --gene-input-file genes.tsv \
  --biogrid-key YOUR_BIOGRID_KEY \
  --david-email you@example.org \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas pathway_enrichment \
  --gene-input-file genes.tsv \
  --biogrid-key YOUR_BIOGRID_KEY \
  --david-email you@example.org \
  --dsigdb-gmt DSigDB_All.gmt \
  --reference-set genome_protein-coding \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

The API key and DAVID email are required. DSigDB cannot run without a usable
GMT file. The command uses STRING score 400 and FDR 0.05 fallbacks.
The `--reference-set` default comes from canonical
`modules.enrichment.reference_set` YAML and is passed to DSigDB/WebGestaltR.

## Processing steps

Detect delimiter/header, deduplicate and sort symbols, dispatch the environment,
create one provider directory, call each provider inside an exception boundary,
save returned tables/networks, and print provider counts or tracebacks.

## Outputs

Provider folders are `omnipath/`, `dsigdb/`, `biogrid/`, `david/`, `dgidb/`,
`enrichr/`, `gprofiler/`, `toppgene/`, and `stringdb/`. Tables use provider
patterns such as `<dataset>_biogrid_interactions.tsv`, `_david_enrichment.tsv`,
`_dgidb_interactions.tsv`, `_enrichr_results.tsv`, `_gprofiler_enrichment.tsv`,
`_toppgene_enrichment.tsv`, and STRING enrichment/interaction/network files.

## QC and logs

Confirm mapped/input genes, provider release, organism, background, correction,
tested terms, API errors, and non-empty outputs separately. A final “completed”
message is not unified success.

## Interpretation

Providers use different databases, backgrounds, mappings, evidence, and
corrections. Their P values are not directly comparable, and provider/term
searches increase multiplicity.

## Common problems

Invalid credentials, missing DSigDB GMT, rate limits, API changes, missing
runtime, unmapped symbols, or one provider failing while others continue.

## Limitations

The provider sequence is fixed and failure-tolerant. The
`modules.enrichment.providers` configuration field does not select providers in
the current dispatcher. Provider exceptions are caught independently, so a
zero process exit or top-level completion does not certify that every provider
succeeded. Inspect each provider's diagnostics and output; distinguish a valid
empty result from a failed request. There is no unified provider-success
manifest. Installation cannot guarantee remote service availability or valid
credentials.

## Scientific references

- [g:Profiler](https://doi.org/10.1093/nar/gkad347)
- [DAVID](https://david.ncifcrf.gov/)
- [BioGRID](https://thebiogrid.org/)
- [STRING](https://string-db.org/)
- [DGIdb](https://www.dgidb.org/)
- [OmniPath](https://omnipathdb.org/)
