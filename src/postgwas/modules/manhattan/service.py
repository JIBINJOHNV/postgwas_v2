"""Resolve canonical Manhattan settings and validate inputs and plotted outputs."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import shlex
import subprocess

from postgwas.config import load_run_configuration_for_module
from postgwas.config.cli_overrides import explicit_overrides
from postgwas.core.paths import configured_output_path, require_nonempty_file, resolve_executable
from postgwas.core.input_validation import validate_once
from postgwas.core.vcf import declared_vcf_tag_definitions, validate_indexed_vcf
from postgwas.core.ui import print_screen_block, screen_field, screen_line

_PLOT_OPTIONS = {
    "nauto": "autosome_count", "min_af": "minimum_af",
    "min_lp": "minimum_neglog10_p", "loglog_pval": "loglog_pvalue",
    "cyto_ratio": "cytoband_ratio", "max_height": "maximum_height",
    "spacing": "chromosome_spacing", "width": "width", "height": "height",
    "fontsize": "font_size", "allelic_shift": "allelic_shift", "csq": "flag_coding",
}


def resolve_plot_inputs(args):
    """Share the read-only direct/pipeline boundary without creating outputs."""
    overrides = explicit_overrides(args, {
        **_PLOT_OPTIONS, "threads": "execution.threads", "memory_gb": "execution.memory_gb",
    })
    if getattr(args, "pdf", None) or getattr(args, "png", None):
        overrides["file_format"] = "png" if getattr(args, "png", None) else "pdf"
    configuration = load_run_configuration_for_module(
        "manhattan", getattr(args, "run_config", None),
        module_overrides=overrides,
    )
    module = configuration.modules.manhattan
    bcftools = resolve_executable(configuration.resources.executables.bcftools, "bcftools")
    rscript = resolve_executable(configuration.resources.executables.rscript, "Rscript")
    script = require_nonempty_file(Path(__file__).parent / "resources" / "assoc_plot.R", "Manhattan R adapter")
    if getattr(args, "pdf", None) and getattr(args, "png", None):
        raise ValueError("--pdf and --png are mutually exclusive")
    required = ["INFO/" + module.allelic_shift_field if module.allelic_shift else "FORMAT/" + module.pvalue_field]
    if module.minimum_af > 0:
        required.append("FORMAT/" + module.frequency_field)
    if module.flag_coding:
        required.append("INFO/" + module.consequence_info_field)
    contract = configuration.modules.formatting.input_contract
    indexed = validate_indexed_vcf(
        args.vcf, args.dataset_id, bcftools,
        genome_build_header="##" + contract.genome_build_metadata,
        supported_genome_builds=contract.supported_genome_builds,
        required_fields=required,
        provenance_headers=contract.provenance_headers.model_dump(),
        expected_genome_build=getattr(args, "genome_build", None),
    )
    pheno = getattr(args, "pheno", None)
    if pheno and pheno != indexed.sample:
        raise ValueError("--pheno %r does not match the VCF sample %r" % (pheno, indexed.sample))
    for field in required:
        category, tag = field.split("/", 1)
        definition = declared_vcf_tag_definitions(indexed.header, category)[tag]
        if field == "INFO/" + module.consequence_info_field:
            continue
        # One LP/AF per association; AS is the two-count binomial input.
        expected_numbers = {"2"} if field == "INFO/" + module.allelic_shift_field else {"1", "A"}
        if definition.get("number") not in expected_numbers or definition.get("type") not in {"Float", "Integer"}:
            raise ValueError("%s must be numeric with Number in %s" % (field, sorted(expected_numbers)))
    try:
        runtime = validate_once(
            (rscript, script),
            {"validator": "manhattan_r_runtime", "environment": {
                key: os.environ.get(key) for key in (
                    "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE",
                    "R_PROFILE", "R_PROFILE_USER", "R_ENVIRON", "R_ENVIRON_USER",
                )
            }},
            lambda: subprocess.run(
                [rscript, str(script), "--check-runtime"], stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, check=True,
            ).stdout.strip(),
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError("Manhattan R runtime check failed; verify the configured Rscript and compatible R libraries: %s" % exc.stdout) from exc
    if module.flag_coding:
        try:
            columns = validate_once(
                (bcftools, indexed.vcf),
                {"validator": "manhattan_consequence_header",
                 "annotation": module.consequence_info_field,
                 "plugins": os.environ.get("BCFTOOLS_PLUGINS")},
                lambda: subprocess.run(
                    [bcftools, "+split-vep", str(indexed.vcf), "--annotation", module.consequence_info_field, "--list"],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True,
                ).stdout,
            )
        except subprocess.CalledProcessError as exc:
            raise RuntimeError("Coding highlighting requires a working bcftools split-vep plugin and structured consequence header: %s" % exc.stderr) from exc
        if module.consequence_field not in columns.split():
            raise ValueError("Configured consequence field is absent from the split-vep annotation header: %s" % module.consequence_field)
    return configuration, indexed, rscript, runtime


def run_assoc_plot_direct(args):
    """Render a checked plot and retain the exact point data and command."""
    configuration, indexed, rscript, runtime = resolve_plot_inputs(args)
    module = configuration.modules.manhattan
    file_format = module.file_format
    explicit_output = getattr(args, file_format, None)
    output = Path(explicit_output).expanduser().resolve() if explicit_output else configured_output_path(
        args.output_directory, module.output_file, dataset_id=args.dataset_id, file_format=file_format,
    )
    paths = {
        "log": configured_output_path(args.output_directory, module.log_file, dataset_id=args.dataset_id),
        "plot_data": configured_output_path(args.output_directory, module.plot_data_file, dataset_id=args.dataset_id),
    }
    if len({output, *paths.values()}) != 3 or indexed.vcf in {output, *paths.values()}:
        raise ValueError("Manhattan output, log, audit table and input paths must be distinct")
    script = Path(__file__).parent / "resources" / "assoc_plot.R"
    command = [rscript, str(require_nonempty_file(script, "Manhattan R adapter"))]
    options = {
        file_format: output, "vcf": indexed.vcf, "genome": indexed.genome_build,
        "pheno": indexed.sample, "bcftools": indexed.bcftools,
        "plot-data": paths["plot_data"], "threads": configuration.execution.threads,
        "significance": module.significance_threshold, "suggestive": module.suggestive_threshold,
        "single-height": module.single_chromosome_height, "genome-height": module.genomewide_height,
        "dpi": module.png_dpi, "lp-field": module.pvalue_field,
        "af-field": module.frequency_field, "as-field": module.allelic_shift_field,
        "consequence-field": module.consequence_field,
        "consequence-info-field": module.consequence_info_field,
        "coding-terms": ",".join(module.coding_terms),
        "chromosome-label-policy": module.chromosome_label_policy,
        "chromosome-prefix": module.chromosome_prefix,
        "axis-label-rows": module.axis_label_rows,
        "caption-fontsize": module.caption_font_size,
        "caption-template": module.caption_template,
        "input-records": indexed.variant_count,
    }
    options.update({key.replace("_", "-"): getattr(module, value) for key, value in _PLOT_OPTIONS.items() if key not in {"allelic_shift", "csq"}})
    command.extend("--%s=%s" % (key, value) for key, value in options.items() if value is not None)
    if module.allelic_shift:
        command.append("--as")
    if module.flag_coding:
        command.append("--csq")
    for path in (output, *paths.values()):
        path.parent.mkdir(parents=True, exist_ok=True)
    print_screen_block("\n".join((
        screen_line("analysis", "Generating Manhattan plot"),
        screen_field("info", "Plot", output, path_value=True),
        screen_field("genetic", "Genome build", indexed.genome_build),
        screen_field("info", "Sample", indexed.sample),
    )))
    with paths["log"].open("w", encoding="utf-8") as log:
        log.write("RESOLVED CONFIGURATION:\n" + json.dumps(module.model_dump(mode="json"), indent=2) + "\n")
        log.write("INPUT RECORDS: %s\nCOMMAND: %s\n" % (indexed.variant_count, shlex.join(command)))
        log.write("RUNTIME: %s\n" % runtime)
        log.flush()
        try:
            subprocess.run(command, stdout=log, stderr=log, text=True, check=True)
            require_nonempty_file(output, "Manhattan plot")
            # PDF and PNG magic bytes are file-format protocol invariants.
            with output.open("rb") as plot:
                signature = plot.read(8)
            if not (signature.startswith(b"%PDF-") if file_format == "pdf" else signature == b"\x89PNG\r\n\x1a\n"):
                raise ValueError("Manhattan output is not a valid %s file" % file_format.upper())
            with require_nonempty_file(paths["plot_data"], "Manhattan point audit").open() as handle:
                rows = 0
                for row in csv.DictReader(handle, delimiter="\t"):
                    values = [float(row[key]) for key in ("pos", "raw_lp", "lp", "chrompos")]
                    if not all(math.isfinite(value) and value >= 0 for value in values) or values[0] < 1:
                        raise ValueError("Invalid coordinate or association statistic in Manhattan point audit")
                    rows += 1
            if not rows or rows > indexed.variant_count:
                raise ValueError("Manhattan point count is empty or exceeds the input record count")
            log.write("VALIDATED: %s plotted records; %s excluded by configured filters or missing LP.\n" % (rows, indexed.variant_count - rows))
        except Exception as exc:
            log.write("FAILED: %s\n" % exc)
            raise RuntimeError("Manhattan plot failed; see %s: %s" % (paths["log"], exc)) from exc
    print_screen_block("\n".join((
        screen_line("success", "Manhattan plot validated"),
        screen_field("count", "Plotted records", f"{rows:,}"),
        screen_field("count", "Excluded records", f"{indexed.variant_count - rows:,}"),
        screen_field("info", "Point audit", paths["plot_data"], path_value=True),
    )))
    return {"plot": str(output), **{key: str(value) for key, value in paths.items()}}
