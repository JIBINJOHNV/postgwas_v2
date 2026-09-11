# Frequently Asked Questions

## Can the pipeline start from a raw summary-statistics table?

No. Run standalone harmonisation first. The downstream pipeline starts from a
harmonised GWAS-VCF.

## Is `pip install` enough for every analysis?

No. It installs the Python package, not the complete external-tool stack or
reference data. Use the [complete Mamba installer](../getting-started/installation.md#install-the-complete-software-stack)
or the [local Docker build](../getting-started/installation.md#docker-installation)
for the software stack. Both routes still require analysis-specific reference
files; see the installation guide for platform support and validation limits.

## Does PostGWAS download reference data automatically?

No universal validated resource bundle is provided. Prepare and version the
resources required by each module, then verify build and ancestry compatibility.

## Should I trust defaults printed in command help?

Use the exported YAML or [Configuration Defaults](../reference/configuration-defaults.md)
for the installed version. Do not rely on copied values from an older guide.

## Can I run one module without the pipeline?

Usually, yes, but you must satisfy its upstream file and scientific contracts
yourself. See [Running Modules Independently](../core/running-modules-independently.md).

## Does a successful process exit mean the analysis is valid?

No. Review compatibility checks, variant retention, warnings, logs, and QC.
External tools cannot always detect an inappropriate LD population or a
mislabeled genome build.

## Are clumped variants causal?

No. LD clumping identifies approximately independent association signals under
chosen thresholds and reference LD. Causal interpretation requires additional
evidence and, commonly, fine-mapping or functional follow-up.

## Can QC create a cleaned replacement VCF?

Do not assume so. The standalone QC assessment reports metrics and can evaluate
a virtual filtered subset; it is not equivalent to the filtering module's
materialized output.
