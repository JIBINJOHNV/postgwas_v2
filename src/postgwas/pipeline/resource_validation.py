"""Explicit adapters from module preflight evidence to the shared file report.

These functions never infer file types from suffixes or discover inputs by
walking arbitrary configuration objects. They report checks already performed
by a successful module preflight; they do not reopen resources. Availability
is deliberately distinguished from parsed content and scientific compatibility.
"""

from __future__ import annotations

from postgwas.core.input_validation import record_file_validation
from postgwas.core.preflight import PipelinePreflightEvidence


def _availability(path, role: str, *, message: str = "", metrics=None) -> None:
    if path is not None:
        record_file_validation(
            path, role, checks=("regular nonempty file",), metrics=metrics,
            message=message or "Availability checked; file contents are not validated by this check.",
        )


def _imputation(resources) -> None:
    for path in resources.reference_files:
        _availability(path, "PRED-LD chromosome reference")
    _availability(resources.pred_ld_script, "PRED-LD script")


def _heritability(resources) -> None:
    _availability(resources.reference.merge_alleles, "LDSC merge-alleles file")
    for path in resources.reference.required_files:
        _availability(path, "LDSC chromosome LD-score or SNP-count companion")


def _covariates(metrics) -> None:
    record_file_validation(
        metrics["covariates_file"], "MAGMA gene-covariate table",
        checks=(
            "all rows: consistent fields and unique gene identifiers",
            "finite numeric values or supported missing token",
            "configured minimum genes and per-property missingness",
            "nonconstant properties",
        ),
        metrics={
            "genes": metrics["covariate_genes"],
            "properties": metrics["properties"],
            "property_missingness": metrics["property_missingness"],
        },
        message="Compatibility with pipeline-generated MAGMA gene results is checked before use.",
    )


def _magmacovar(resources) -> None:
    _covariates(resources["covariates"])


def _pops(resources) -> None:
    _, features, annotation, controls = resources
    record_file_validation(
        features["rows_path"], "PoPS feature gene-order file",
        checks=("nonempty unique gene names", "matrix row-count agreement"),
        metrics={"genes": features["row_count"]},
    )
    for chunk in features["chunks"]:
        record_file_validation(
            chunk["columns_path"], "PoPS feature-column names",
            checks=("unique feature names across all chunks", "matrix column-count agreement"),
            metrics={"chunk": chunk["chunk"], "features": chunk["columns"]},
        )
        record_file_validation(
            chunk["matrix_path"], "PoPS feature matrix",
            checks=("matrix shape agrees with companion files", "all values: finite numeric data"),
            metrics={key: chunk[key] for key in ("chunk", "rows", "columns", "dtype")},
        )
    record_file_validation(
        annotation["path"], "PoPS gene annotation",
        checks=(
            "configured columns and unique gene identifiers",
            "finite nonnegative TSS and nonempty chromosome labels",
            "configured chromosome coverage and feature-gene compatibility",
        ),
        metrics={"genes": annotation["gene_count"], "chromosomes": annotation["chromosomes"]},
    )
    for key, role in (("feature_subset", "PoPS feature subset"), ("controls", "PoPS control features")):
        values = controls[key]
        if values is not None:
            record_file_validation(
                values["path"], role,
                checks=("nonempty unique names", "membership in configured feature universe"),
                metrics={"features": values["count"]},
            )


def _kpops(resources) -> None:
    _, values = resources
    record_file_validation(
        values["annotation"], "K-POPS gene annotation",
        checks=(
            "required columns and unique gene identifiers", "finite nonnegative TSS",
            "configured chromosome and anchor-gene compatibility",
        ),
        metrics={"genes": values["annotation_gene_count"]},
    )
    record_file_validation(
        values["kernel_genes"], "K-POPS kernel gene-order file",
        checks=("unique identifiers and configured minimum genes", "kernel genes present in annotation"),
        metrics={"genes": values["kernel_gene_count"]},
    )
    record_file_validation(
        values["kernel"], "K-POPS kernel matrix",
        checks=("binary dimensions agree with gene-order file", "all values finite", "configured symmetry tolerance"),
        metrics={"genes": values["kernel_gene_count"]},
        message="Symmetry is checked; positive semidefiniteness is not established by this check.",
    )


def _caldera(resources) -> None:
    _, values = resources
    for key, role in (
        ("coding_variants", "CALDERA coding variants"),
        ("gene_locations", "CALDERA gene locations"),
        ("model", "CALDERA trained model"),
        ("upstream_script", "CALDERA upstream script"),
        ("adapter", "CALDERA adapter"),
    ):
        _availability(values[key], role)


def _flames(resources) -> None:
    _, values = resources
    _availability(values["model"], "FLAMES trained model")
    _availability(values["script"], "FLAMES script")
    record_file_validation(
        values["features"], "FLAMES feature manifest",
        checks=("nonempty unique feature names",),
        metrics={"features": len(values["feature_names"])},
    )
    for path in dict.fromkeys((*values["annotation_inventory"].values(), *values["annotation_files"].values())):
        _availability(path, "FLAMES annotation resource")
    for key, role in (("cadd_file", "CADD score file"), ("cadd_index", "CADD tabix index")):
        _availability(values[key], role)
    if values["vep_cache"] is not None:
        record_file_validation(
            values["vep_cache"], "VEP cache directory",
            checks=("directory exists and contains a file",),
            message="Configured build declarations agree; cache contents and version are not established by this check.",
        )


def _ld_clump(resources) -> None:
    if resources.reference_manifest is not None:
        module = resources.configuration.modules.ld_clumping
        record_file_validation(
            resources.reference_directory / module.reference.manifest_filename,
            "LD-clumping reference manifest",
            checks=(
                "schema-valid manifest", "declared build and population compatibility",
                "configured orientation, column and allele-order contract",
                "reference window, r2 and MAF coverage declarations",
            ),
            metrics=resources.reference_manifest.model_dump(mode="json"),
            message="Manifest declarations checked; this does not prove every chromosome table's content.",
        )
        record_file_validation(
            None, "LD-clumping chromosome reference content", status="deferred",
            message="Validate candidate-dependent chromosome tables, indexes and variant inventories using the configured missing-reference policy before LD lookup.",
        )


def _mixer(resources) -> None:
    # These exact identities were captured after expansion of the configured
    # chromosome patterns. No assumption is made about an opaque LD binary.
    for identity in resources.file_identities:
        _availability(identity.path, "MiXeR input resource or executable")
    if resources.gsa_resources is not None:
        for key in ("baseline_go", "model_go", "test_go"):
            record_file_validation(
                resources.gsa_resources[key], "GSA-MiXeR gene-set file",
                checks=("configured header columns", "first data row present"),
                message="Header and first record checked; this check does not inspect all ontology rows.",
            )


def _single_cell(resources) -> None:
    for method, evidence in resources["methods"].items():
        if method == "magma_celltype":
            _covariates(evidence)
        elif method == "scdrs":
            record_file_validation(
                evidence.h5ad_file, "scDRS H5AD expression atlas",
                checks=(
                    "unique cell and gene names", "all X values finite and nonnegative",
                    "configured matrix-state and annotation requirements",
                    "configured expression-filter feasibility",
                ),
                metrics=evidence.h5ad_summary,
            )
            _availability(evidence.gene_identifier_map_file, "scDRS gene-identifier crosswalk")
            if evidence.covariate_file is not None:
                record_file_validation(
                    evidence.covariate_file, "scDRS cell covariates",
                    checks=("cell identifiers and configured covariate requirements",),
                    metrics=evidence.covariate_summary,
                )
        elif method == "ldsc_celltype":
            record_file_validation(
                evidence.ldcts_file, "LDSC cell-type manifest",
                checks=("all rows: configured manifest syntax", "unique labels and valid resource-prefix lists"),
                metrics={"cell_types": len(evidence.entries)},
            )
            for entry in evidence.reference_inventory:
                _availability(
                    entry["path"], "LDSC cell-type %s" % entry["role"],
                    metrics={"chromosome": entry["chromosome"], "resource": entry["resource"]},
                )
            _availability(evidence.merge_alleles_file, "LDSC cell-type merge-alleles file")
        else:
            raise RuntimeError("No resource-report adapter for single-cell method %r" % method)


def _already_recorded(resources) -> None:
    """File validators record their evidence directly; do not rescan or duplicate it."""


# Module names are registry protocol identities, not defaults or scientific
# parameters. Coverage is checked against all registered pipeline preflights.
RESOURCE_REPORTERS = {
    "sumstat_filter": _already_recorded,
    "post_imputation_filter": _already_recorded,
    "formatter": _already_recorded,
    "manhattan": _already_recorded,
    "qc_summary": _already_recorded,
    "annot_ldblock": _already_recorded,
    "magma": _already_recorded,
    "gcta_cojo": _already_recorded,
    "gcta_gene": _already_recorded,
    "finemap": _already_recorded,
    "imputation": _imputation,
    "heritability": _heritability,
    "ld_clump": _ld_clump,
    "magmacovar": _magmacovar,
    "pops": _pops,
    "kpops": _kpops,
    "caldera": _caldera,
    "flames": _flames,
    "mixer": _mixer,
    "single_cell": _single_cell,
}


def record_pipeline_resource_validation(
    module_name: str,
    validated: PipelinePreflightEvidence,
) -> None:
    """Publish exact already-completed file checks in the active module scope."""
    if validated.module != module_name:
        raise ValueError("Resource-report module does not match its preflight evidence")
    reporter = RESOURCE_REPORTERS.get(module_name)
    if reporter is not None:
        reporter(validated.resources)


__all__ = ["RESOURCE_REPORTERS", "record_pipeline_resource_validation"]
