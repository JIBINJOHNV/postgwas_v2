# Command Reference

Use `postgwas --help` from the installed checkout as the authoritative command
inventory. The current public commands map to Wiki pages as follows.

| Command | Purpose | Guide |
|---|---|---|
| `harmonisation` | Raw summary statistics to GWAS-VCF | [Harmonisation](../harmonisation/overview.md) |
| `sumstat_filter` | Materialized GWAS-VCF filtering | [Filtering](../modules/filtering.md) |
| `formatter` | Tool-specific input export | [Formatting](../modules/formatting.md) |
| `imputation` | PRED-LD summary-statistic imputation | [Imputation](../modules/imputation.md) |
| `annot_ldblock` | LD-block INFO annotation | [LD Annotation](../modules/ld-annotation.md) |
| `ld_clump` | Independent signals and loci | [LD Clumping](../modules/ld-clumping.md) |
| `heritability` | Single-trait LDSC h² | [LDSC Heritability](../modules/ldsc.md) |
| `finemap` | SuSiE-RSS or FINEMAP | [Fine Mapping](../modules/fine-mapping.md) |
| `magma` | MAGMA gene/gene-set analysis | [MAGMA](../modules/magma.md) |
| `gcta_gene` | GCTA fastBAT/mBAT-combo | [GCTA Gene Analysis](../modules/gcta-gene.md) |
| `gcta_cojo` | GCTA-COJO conditional and joint analysis | `docs/modules/gcta_cojo/README.md` in the repository |
| `magmacovar` | MAGMA gene-property analysis | [MAGMAcovar](../modules/magmacovar.md) |
| `single_cell` | MAGMA cell typing, scDRS, LDSC cell-type analysis | [Single-Cell Integration](../modules/single-cell.md) |
| `pops` | PoPS gene prioritization | [PoPS](../modules/pops.md) |
| `kpops` | Kernel-based K-POPS gene prioritization | [K-POPS](../../modules/kpops.md) |
| `caldera` | PoPS plus fine-mapped credible sets | [CALDERA](../../modules/caldera.md) |
| `flames` | Validated FLAMES effector-gene prioritization | [FLAMES](../modules/flames.md) |
| `mixer` | MiXeR/GSA-MiXeR | [MiXeR](../modules/mixer.md) |
| `manhattan` | Association plots | [Manhattan Plots](../modules/manhattan.md) |
| `qc` | Compact GWAS-VCF QC | [QC Summary](../modules/qc-summary.md) |
| `pathway_enrichment` | Multi-provider gene-list enrichment | [Pathway Enrichment](../modules/pathway-enrichment.md) |
| `pipeline` | Dependency planning and execution | [Pipeline Workflow](../core/pipeline-workflow.md) |
| `config` | Validate/show/export configuration | [Configuration](../core/configuration.md) |
| `resources` | Install and revalidate pinned resource bundles | [Reference Resources](reference-resources.md) |
| `--validate` | Input-to-VCF concordance | [Validation Reference](validation.md) |

The command name is not always the configuration module name. `sumstat_filter`
resolves to `filtering`, `formatter` to `formatting`, `annot_ldblock` to
`ld_annotation`, `ld_clump` to `ld_clumping`, `finemap` to `fine_mapping`,
`heritability` to `ldsc`, `qc` to `qc_summary`, and `pathway_enrichment` to
`enrichment`. Use the configuration name with `postgwas config --module`.

For any command:

```console
postgwas COMMAND --help
```

For configuration:

```console
postgwas config validate --config run.yaml
postgwas config show --config run.yaml --module mixer
postgwas config export --module mixer --style full --output mixer.yaml
```
