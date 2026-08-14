# Pathway Enrichment

## Purpose

The standalone command queries multiple online pathway, drug-gene, and
interaction providers for one list of human gene symbols.

## What the analysis does

It reads unique symbols from the first column, switches to a dedicated
`enricher` micromamba environment, and independently attempts OmniPath,
DSigDB/WebGestaltR, BioGRID, DAVID, DGIdb, Enrichr, g:Profiler, ToppGene, and
STRINGdb.

## When to use it

Use it for exploratory follow-up of a prespecified gene list. Define the list
and background independently of results and acknowledge selection/multiplicity.

## Input requirements

A delimited file with a header and gene symbols in the first column; BioGRID
API key; DAVID-registered email; optional DSigDB GMT; network access; an
`enricher` micromamba environment with Python/R dependencies; dataset ID and
output directory.

## Command

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
environment, unmapped symbols, or one provider failing while others continue.

## Limitations

The workflow is fixed, failure-tolerant, online, and lacks a unified manifest.
YAML provider selection is not implemented. Confirm that each required provider
returned results before treating the run as complete.

## Scientific references

- [g:Profiler](https://doi.org/10.1093/nar/gkad347)
- [DAVID](https://david.ncifcrf.gov/)
- [BioGRID](https://thebiogrid.org/)
- [STRING](https://string-db.org/)
- [DGIdb](https://www.dgidb.org/)
- [OmniPath](https://omnipathdb.org/)
