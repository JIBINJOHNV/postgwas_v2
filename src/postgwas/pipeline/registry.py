"""Canonical registry of PostGWAS modules.

This is the only place where module descriptions, dependencies, parser
components, runner entry points, and required CLI inputs are declared.
References are strings so listing or planning modules does not import optional
scientific dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Callable, Iterable
import argparse

from postgwas.core.errors import PipelinePlanningError


@dataclass(frozen=True)
class RequiredOption:
    dest: str
    flag: str


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    description: str
    cli_entrypoint: str | None = None
    command_name: str | None = None
    dependencies: tuple[str, ...] = ()
    parser_factories: tuple[str, ...] = ()
    pipeline_example_factory: str | None = None
    pipeline_supplied_options: tuple[str, ...] = ()
    required_options: tuple[RequiredOption, ...] = ()
    preflight: str | None = None
    runner: str | None = None
    pipeline_enabled: bool = True
    unavailable_reason: str | None = None
    internal: bool = False
    pipeline_target: bool = True


def _options(*names: str) -> tuple[RequiredOption, ...]:
    preferred_flags = {
        "dataset_id": "--dataset-id",
        "output_directory": "--output-directory",
        "genome_build": "--genome-build",
    }
    return tuple(
        RequiredOption(
            name.replace("-", "_"),
            preferred_flags.get(name, "--" + name.replace("_", "-")),
        )
        for name in names
    )


COMMON = (
    "postgwas.cli.compute:get_compute_parser",
    "postgwas.cli.common:get_inputvcf_parser",
    "postgwas.cli.common:get_common_out_parser",
)

PIPELINE_REQUIRED_OPTIONS = _options("vcf", "dataset_id", "output_directory")


def get_qc_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("QC Summary Arguments")
    group.add_argument("--qc-output-prefix")
    return parser


MODULES = (
    ModuleSpec(
        "config",
        "Validate and inspect resolved PostGWAS configuration.",
        cli_entrypoint="postgwas.config.cli:main",
        pipeline_enabled=False,
        pipeline_target=False,
        unavailable_reason="configuration inspection is not a scientific pipeline step",
    ),
    ModuleSpec(
        "resources",
        "Install and validate pinned scientific reference bundles.",
        cli_entrypoint="postgwas.resources.cli:main",
        pipeline_enabled=False,
        pipeline_target=False,
        unavailable_reason="resource installation is not a scientific pipeline step",
    ),
    ModuleSpec(
        "harmonisation",
        "Standardise raw summary statistics and create GWAS-VCF artifacts.",
        cli_entrypoint="postgwas.modules.harmonisation.cli:main",
        pipeline_enabled=False,
        unavailable_reason="full-pipeline adapter is not implemented; use `postgwas harmonisation`",
    ),
    ModuleSpec(
        "sumstat_filter",
        "Filter summary statistics using explicit QC policies.",
        cli_entrypoint="postgwas.modules.filtering.cli:main",
        parser_factories=COMMON + (
            "postgwas.cli.common:get_genome_build_parser",
            "postgwas.cli.common:get_common_sumstat_filter_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_sumstat_filter_runner",
    ),
    ModuleSpec(
        "post_imputation_filter",
        "Apply the configured QC policy after imputation.",
        dependencies=("imputation",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_genome_build_parser",
            "postgwas.cli.common:get_common_sumstat_filter_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_sumstat_filter_runner",
        internal=True,
    ),
    ModuleSpec(
        "annot_ldblock",
        "Annotate variants with population-specific LD blocks.",
        cli_entrypoint="postgwas.modules.ld_annotation.cli:main",
        parser_factories=COMMON + (
            "postgwas.cli.common:get_genome_build_parser",
            "postgwas.cli.common:get_annot_ldblock_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_annot_ldblock_runner",
    ),
    ModuleSpec(
        "formatter",
        "Convert GWAS-VCF into validated downstream-tool artifacts.",
        cli_entrypoint="postgwas.modules.formatting.cli:main",
        parser_factories=COMMON + (
            "postgwas.cli.common:get_formatter_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_formatter_runner",
    ),
    ModuleSpec(
        "imputation",
        "Impute missing summary statistics.",
        cli_entrypoint="postgwas.modules.imputation.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_common_imputation_parser",
            "postgwas.cli.common:get_population_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_imputation_runner",
    ),
    ModuleSpec(
        "ld_clump",
        "Identify independent significant variants and genomic loci.",
        cli_entrypoint="postgwas.modules.ld_clumping.cli:main",
        dependencies=("annot_ldblock",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_genome_build_parser",
            "postgwas.cli.common:get_plink_binary_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
            "postgwas.cli.common:get_population_parser",
            "postgwas.cli.common:get_ld_clump_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_ld_clump_runner",
    ),
    ModuleSpec(
        "finemap",
        "Run SuSiE or FINEMAP fine-mapping.",
        cli_entrypoint="postgwas.modules.fine_mapping.cli:main",
        dependencies=("ld_clump", "formatter"),
        parser_factories=COMMON + (
            "postgwas.modules.fine_mapping.arguments:get_finemap_common_parser",
            "postgwas.modules.fine_mapping.arguments:get_common_susie_arguments",
            "postgwas.modules.fine_mapping.arguments:get_common_finemap_finemap_arguments",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.fine_mapping.cli:get_finemap_pipeline_examples"
        ),
        pipeline_supplied_options=("locus_file",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_finemap_runner",
    ),
    ModuleSpec(
        "magma",
        "Run MAGMA gene and gene-set association analysis.",
        cli_entrypoint="postgwas.modules.magma.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + (
            "postgwas.modules.magma.cli:get_magma_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.magma.cli:get_magma_pipeline_examples"
        ),
        pipeline_supplied_options=(
            "snp_location_file", "p_value_file", "variant_id_type",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_magma_runner",
    ),
    ModuleSpec(
        "gcta_cojo",
        "Run GCTA-COJO conditional, joint, or stepwise association analysis.",
        cli_entrypoint="postgwas.modules.gcta_cojo.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + (
            "postgwas.modules.gcta_cojo.cli:get_gcta_cojo_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.gcta_cojo.cli:get_gcta_cojo_pipeline_examples"
        ),
        pipeline_supplied_options=("gcta_cojo_input_file", "variant_id_type"),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_gcta_cojo_runner",
    ),
    ModuleSpec(
        "gcta_gene",
        "Run GCTA fastBAT gene/segment/set or mBAT-combo association analysis.",
        cli_entrypoint="postgwas.modules.gcta_gene.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + (
            "postgwas.modules.gcta_gene.cli:get_gcta_gene_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.gcta_gene.cli:get_gcta_gene_pipeline_examples"
        ),
        pipeline_supplied_options=("gcta_input_file",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_gcta_gene_runner",
    ),
    ModuleSpec(
        "magmacovar",
        "Run MAGMA gene-property analysis.",
        cli_entrypoint="postgwas.modules.magmacovar.cli:main",
        dependencies=("magma",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_common_magma_covar_parser",
        ),
        pipeline_supplied_options=("magma_gene_results_file",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight=(
            "postgwas.modules.magmacovar.service:preflight_magmacovar_pipeline"
        ),
        runner="postgwas.pipeline.runners:run_magmacovar_runner",
    ),
    ModuleSpec(
        "single_cell",
        "Identify GWAS-associated cell types using single-cell expression data.",
        cli_entrypoint="postgwas.modules.single_cell.cli:main",
        dependencies=("magma",),
        parser_factories=COMMON + (
            "postgwas.modules.magmacovar.cli:get_magma_covar_gene_results_parser",
            "postgwas.modules.single_cell.cli:get_single_cell_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.single_cell.cli:get_single_cell_pipeline_examples"
        ),
        pipeline_supplied_options=(
            "magma_gene_results_file", "scdrs_magma_gene_results_file",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight=(
            "postgwas.modules.single_cell.service:preflight_single_cell_pipeline"
        ),
        runner="postgwas.pipeline.runners:run_single_cell_runner",
    ),
    ModuleSpec(
        "pops",
        "Prioritise genes with PoPS.",
        cli_entrypoint="postgwas.modules.pops.cli:main",
        dependencies=("magma",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_common_pops_parser",
            "postgwas.modules.pops.cli:get_pops_pipeline_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.pops.cli:get_pops_pipeline_examples"
        ),
        pipeline_supplied_options=("magma_association_prefix",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight="postgwas.modules.pops.service:preflight_pops_pipeline",
        runner="postgwas.pipeline.runners:run_pops_runner",
    ),
    ModuleSpec(
        "kpops",
        "Prioritise genes with kernel-based K-POPS.",
        cli_entrypoint="postgwas.modules.kpops.cli:main",
        dependencies=("magma",),
        parser_factories=COMMON + (
            "postgwas.modules.kpops.cli:get_kpops_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.kpops.cli:get_kpops_pipeline_examples"
        ),
        pipeline_supplied_options=("magma_association_prefix",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight="postgwas.modules.kpops.service:preflight_kpops_pipeline",
        runner="postgwas.pipeline.runners:run_kpops_runner",
    ),
    ModuleSpec(
        "caldera",
        "Prioritise causal genes with PoPS and fine-mapped credible sets.",
        cli_entrypoint="postgwas.modules.caldera.cli:main",
        dependencies=("pops", "finemap"),
        parser_factories=COMMON + (
            "postgwas.modules.caldera.cli:get_caldera_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.caldera.cli:get_caldera_pipeline_examples"
        ),
        pipeline_supplied_options=("pops_file", "credible_set_file"),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight="postgwas.modules.caldera.service:preflight_caldera_pipeline",
        runner="postgwas.pipeline.runners:run_caldera_runner",
    ),
    ModuleSpec(
        "flames",
        "Integrate fine-mapping, MAGMA, gene-property, and PoPS evidence.",
        cli_entrypoint="postgwas.modules.flames.cli:main",
        dependencies=("finemap", "magmacovar", "pops"),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_flames_common_parser",
        ),
        pipeline_example_factory=(
            "postgwas.modules.flames.cli:get_flames_pipeline_examples"
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        preflight="postgwas.modules.flames.service:preflight_flames_pipeline",
        runner="postgwas.pipeline.runners:run_flames_runner",
    ),
    ModuleSpec(
        "heritability",
        "Estimate SNP heritability with LDSC.",
        cli_entrypoint="postgwas.modules.ldsc.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + ("postgwas.cli.common:get_ldsc_common_parser",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_heritability_runner",
    ),
    ModuleSpec(
        "manhattan",
        "Generate Manhattan and QQ plots.",
        cli_entrypoint="postgwas.modules.manhattan.cli:main",
        parser_factories=COMMON + ("postgwas.cli.common:get_assoc_plot_parser",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_manhattan_runner",
    ),
    ModuleSpec(
        "qc_summary",
        "Generate a GWAS-VCF QC summary.",
        cli_entrypoint="postgwas.modules.qc_summary.cli:main",
        command_name="qc",
        parser_factories=COMMON + (
            "postgwas.pipeline.registry:get_qc_parser",
            "postgwas.cli.common:sumstat_summary_arg_parser",
            "postgwas.cli.common:get_bcftools_binary_parser",
        ),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_qc_summary_runner",
    ),
    ModuleSpec(
        "mixer",
        "Run single-trait MiXeR architecture and/or GSA gene-set analysis.",
        cli_entrypoint="postgwas.modules.mixer.cli:main",
        dependencies=("formatter",),
        parser_factories=COMMON + (
            "postgwas.cli.common:get_genome_build_parser",
            "postgwas.modules.mixer.cli:get_mixer_parser",
        ),
        pipeline_supplied_options=("mixer_input_file",),
        required_options=PIPELINE_REQUIRED_OPTIONS,
        runner="postgwas.pipeline.runners:run_mixer_runner",
    ),
    ModuleSpec(
        "allele_orientation",
        "Audit and repair allele orientation against a compatible reference.",
        pipeline_enabled=False,
        unavailable_reason="the prototype is not scientifically safe for pipeline use",
    ),
    ModuleSpec(
        "enrichment",
        "Run pathway and interaction enrichment analyses.",
        cli_entrypoint="postgwas.modules.enrichment.main:main",
        command_name="pathway_enrichment",
        pipeline_enabled=False,
        unavailable_reason="enrichment currently runs as a standalone terminal analysis",
    ),
    ModuleSpec(
        "pipeline",
        "Plan and execute a validated multi-module workflow.",
        cli_entrypoint="postgwas.pipeline.cli:main",
        pipeline_enabled=False,
        pipeline_target=False,
        unavailable_reason="the pipeline orchestrator cannot be selected as its own target",
    ),
)


class ModuleRegistry:
    def __init__(self, specs: Iterable[ModuleSpec] = MODULES):
        self._specs = {}
        for spec in specs:
            if spec.name in self._specs:
                raise ValueError("Duplicate module registration: %s" % spec.name)
            self._specs[spec.name] = spec
        self._validate_dependencies()

    def _validate_dependencies(self) -> None:
        for spec in self._specs.values():
            unknown = [name for name in spec.dependencies if name not in self._specs]
            if unknown:
                raise ValueError("%s has unknown dependencies: %s" % (spec.name, ", ".join(unknown)))

    def names(self, include_internal: bool = False) -> tuple[str, ...]:
        return tuple(
            name
            for name, spec in self._specs.items()
            if spec.pipeline_target and (include_internal or not spec.internal)
        )

    def get(self, name: str) -> ModuleSpec:
        try:
            return self._specs[name]
        except KeyError as exc:
            choices = ", ".join(self.names())
            raise PipelinePlanningError(
                "Unknown module '%s'. Available modules: %s" % (name, choices)
            ) from exc

    def commands(self) -> dict[str, ModuleSpec]:
        return {
            (spec.command_name or spec.name): spec
            for spec in self._specs.values()
            if spec.cli_entrypoint and not spec.internal
        }

    def require_pipeline_enabled(self, name: str) -> ModuleSpec:
        spec = self.get(name)
        if not spec.pipeline_enabled or not spec.runner:
            reason = spec.unavailable_reason or "no pipeline runner is registered"
            raise PipelinePlanningError("Module '%s' is unavailable in pipeline mode: %s" % (name, reason))
        return spec


def resolve_reference(reference: str) -> Callable:
    module_name, separator, attribute = reference.partition(":")
    if not separator or not module_name or not attribute:
        raise ValueError("Invalid callable reference: %r" % reference)
    return getattr(import_module(module_name), attribute)


REGISTRY = ModuleRegistry()
