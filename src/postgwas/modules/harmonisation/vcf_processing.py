"""
GWAS summary statistics -> VCF, annotation and liftover, driven by bcftools.

Runtime choices are resolved from the canonical policy registry in
``harmonisation/policies.py``. Passing ``policies=None`` loads that registry's
validated defaults.

Policy keys read by this module::

    chromosome.allowed_after_split  ordered chromosome merge contract
    vcf.liftover_swap                exclude | keep
    vcf.liftover_warn_fraction       0.10
    vcf.liftover_critical_fraction   0.30
    vcf.liftover_fail_fraction       0.50
    vcf.variant_drop_warn_fraction   0.20
    vcf.min_valid_size_bytes         100
    vcf.sort_memory_mb_per_thread    512
    vcf.concat_max_workers           3
    vcf.concat_max_attempts          2
    vcf.on_merge_failure             continue | fail
"""

from datetime import datetime
from pathlib import Path
import re
import shlex
import shutil
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, Mapping, Sequence


from .policies import default_policies
from postgwas.core.paths import configured_output_path, resolve_executable
from postgwas.core.processes import run_checked_command
from postgwas.core.ui.screen import screen_line
from postgwas.core.vcf import (
    annotation_vcf_fields,
    count_indexed_vcf_records,
    declared_vcf_tag_definitions,
    read_vcf_header,
)

LIFTOVER_SWAP_CHOICES = ("exclude", "keep")
ON_MERGE_FAILURE_CHOICES = ("continue", "fail")

_LIFTOVER_SUMMARY = re.compile(
    r"^Lines\s+total/swapped/reference added/rejected:\s*"
    r"(\d+)/(\d+)/(\d+)/(\d+)\s*$",
    re.MULTILINE,
)
_NORM_SUMMARY = re.compile(
    r"^Lines\s+total/split/joined/realigned/mismatch_removed/"
    r"dup_removed/skipped:\s*"
    r"(\d+)/(\d+)/(\d+)/(\d+)/(\d+)/(\d+)/(\d+)\s*$",
    re.MULTILINE,
)


# ============================================================
# Typed errors (library code must not call sys.exit)
# ============================================================
class MissingExecutableError(RuntimeError):
    """A required external binary is not on PATH."""


class VcfMergeError(RuntimeError):
    """A required chromosome input or merged VCF failed integrity validation."""


class LiftoverFailureRateError(RuntimeError):
    """Liftover exceeded its failure limit or produced no usable variants."""


# ============================================================
# Small helpers
# ============================================================
def _policy_value(policies, key: str):
    """Read a resolved value from the supplied or canonical policy registry."""
    return (policies if policies is not None else default_policies()).get(key)


def _emit(message: str, logger=None, level: str = "info") -> None:
    """Send one audit record to the structured logger, or to stdout standalone.

    Workers must never write to stdout when a logger is available, so ``print``
    is only the last resort (the concat step runs in the parent process).
    """
    if logger is not None:
        if level == "warn":
            logger.warn(message)
        elif level == "error":
            logger.error(message)
        else:
            logger.info(message)
    if logger is None:
        print(screen_line(
            "error" if level == "error" else "warning" if level == "warn" else "info",
            message,
            indent=4,
        ))


def _vcf_metadata_line(name: str, value: Any) -> str:
    """Render one safe, quoted custom VCF metadata line."""
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", str(name)):
        raise ValueError("Invalid configured VCF metadata name: %r" % name)
    text = str(value).replace("\r", " ").replace("\n", " ").replace("\t", " ")
    text = text.replace("\\", "\\\\").replace('"', '\\"')
    return '##%s="%s"' % (name, text)


def _render_provenance_headers(
    provenance: Mapping[str, Any] | None,
    *,
    vcf_config: Mapping[str, Any],
    input_build: str,
    output_build: str,
    merge_group: str,
) -> list[str]:
    """Render configured provenance in the same pass as VCF concatenation."""
    if provenance is None:
        return []
    settings = dict(vcf_config["provenance"])
    values = dict(provenance)
    values.update({
        "vcf_created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input_genome_build": input_build,
        "output_genome_build": output_build,
        "liftover": (
            "Applied: %s to %s" % (input_build, output_build)
            if output_build != input_build
            else "Not applied; records failed target liftover"
            if merge_group == "notlifted"
            else "Not applied"
        ),
    })
    missing = str(settings["missing_value"])
    return [
        _vcf_metadata_line(header_name, values.get(logical_name, missing))
        for logical_name, header_name in dict(settings["headers"]).items()
    ]


def require_binaries(
    executables: Mapping[str, str],
    plugins: Sequence[str] = (),
    logger=None,
) -> Dict[str, str]:
    """Resolve configured programs and verify configured bcftools plugins."""
    found = {
        name: resolve_executable(
            value, "%s executable" % name, error_type=MissingExecutableError,
        )
        for name, value in executables.items()
    }
    if plugins:
        available = set(run_checked_command(
            [found["bcftools"], "plugin", "-l"],
            "Reading available bcftools plugins",
            logger=logger,
            error_type=MissingExecutableError,
        ).split())
        missing_plugins = sorted(set(plugins) - available)
        if missing_plugins:
            message = (
                "Required bcftools plugin(s) not available: %s. Install the plugins "
                "for this bcftools build and set BCFTOOLS_PLUGINS when they are not "
                "in its standard plugin directory."
                % ", ".join(missing_plugins)
            )
            if logger is not None:
                logger.error(message)
            raise MissingExecutableError(message)
    return found


# ============================================================
# HELPER: Run bcftools step
# ============================================================
def _run_bcftools_step(cmd, transcript_file, step_name, chromosome, logger=None):
    """Run one configured VCF command and return its appended transcript block."""
    transcript_path = Path(transcript_file)
    start = transcript_path.stat().st_size if transcript_path.exists() else 0
    run_checked_command(
        cmd,
        "Chromosome %s VCF step %s; transcript: %s"
        % (chromosome, step_name, transcript_file),
        logger=logger,
        error_type=RuntimeError,
        stdout_path=transcript_file,
        append_stdout=True,
        stderr_to_stdout=True,
        stdout_header=(
            "\n[chr%s] === %s ===\ncommand=%s\n"
            % (chromosome, step_name, shlex.join(str(part) for part in cmd))
        ),
    )
    try:
        with transcript_path.open("rb") as handle:
            handle.seek(start)
            return handle.read().decode("utf-8", errors="replace")
    except OSError as exc:
        raise RuntimeError(
            "Chromosome %s VCF step %s completed but its audit transcript "
            "could not be read from %s: %s"
            % (chromosome, step_name, transcript_path, exc)
        ) from exc


# ============================================================
# HELPER: Variant count (safe)
# ============================================================
def get_vcf_variant_count(vcf_path: str, bcftools: str, logger=None):
    """Return the number of records in an indexed VCF, or None.

    A failed ``bcftools index -n`` is recorded with stderr so an unreadable or
    unindexed file cannot disable count-based QC without an explanation.
    """
    try:
        return count_indexed_vcf_records(
            vcf_path,
            bcftools,
            logger=logger,
            error_type=RuntimeError,
        )
    except RuntimeError as exc:
        _emit(
            "The VCF variant count is unavailable: %s" % exc,
            logger=logger, level="warn",
        )
        return None


def _resolve_liftover_tag_arguments(
    header: str,
    vcf_config: Mapping[str, object],
) -> tuple[list[str], dict[str, list[str]]]:
    """Validate canonical allele roles and build explicit plugin tag lists.

    Population-frequency fields are derived from the annotation columns, so a
    panel field cannot be attached to the source-build VCF and then omitted
    from allele-aware liftover by a second, independently maintained list.
    """
    roles = dict(vcf_config["liftover_tag_roles"])
    external_frequency_fields = annotation_vcf_fields(
        vcf_config["external_frequency_columns"], "INFO",
    )
    frequency_fields = list(dict.fromkeys(
        [str(value) for value in roles["allele_frequency"]]
        + external_frequency_fields
    ))
    effect_fields = list(dict.fromkeys(
        str(value) for value in roles["signed_effect"]
    ))
    definitions = {
        category: declared_vcf_tag_definitions(header, category)
        for category in ("INFO", "FORMAT")
    }
    problems = []
    for role, fields in (
        ("allele-frequency", frequency_fields),
        ("signed-effect", effect_fields),
    ):
        for field in fields:
            category, tag = field.split("/", 1)
            metadata = definitions[category].get(tag)
            if metadata is None:
                problems.append("%s is not declared" % field)
                continue
            if metadata.get("number") != "A":
                problems.append(
                    "%s has Number=%s, expected Number=A for its %s role"
                    % (field, metadata.get("number", "missing"), role)
                )
            if metadata.get("type") != "Float":
                problems.append(
                    "%s has Type=%s, expected Type=Float for its %s role"
                    % (field, metadata.get("type", "missing"), role)
                )
    if problems:
        raise ValueError(
            "Annotated VCF allele-dependent fields are incompatible with safe "
            "liftover: %s." % "; ".join(problems)
        )

    def plugin_field(field: str) -> str:
        category, tag = field.split("/", 1)
        return ("FMT" if category == "FORMAT" else category) + "/" + tag

    resolved = {
        "allele_frequency": frequency_fields,
        "signed_effect": effect_fields,
    }
    return [
        "--af-tags", ",".join(plugin_field(field) for field in frequency_fields),
        "--es-tags", ",".join(plugin_field(field) for field in effect_fields),
    ], resolved


def _parse_liftover_audit(transcript: str) -> dict[str, dict[str, int]]:
    """Parse the pinned plugin and bcftools norm accounting summaries."""
    if not isinstance(transcript, str):
        raise RuntimeError("The liftover command did not return a readable transcript.")
    liftover_matches = _LIFTOVER_SUMMARY.findall(transcript)
    norm_matches = _NORM_SUMMARY.findall(transcript)
    if len(liftover_matches) != 1 or len(norm_matches) != 2:
        raise RuntimeError(
            "The liftover transcript must contain one plugin accounting summary "
            "and two normalization summaries; found %d and %d."
            % (len(liftover_matches), len(norm_matches))
        )
    total, swapped, reference_added, rejected = map(int, liftover_matches[0])
    norm_keys = (
        "total", "split", "joined", "realigned", "mismatch_removed",
        "duplicate_removed", "skipped",
    )
    normalization, deduplication = [
        dict(zip(norm_keys, map(int, match))) for match in norm_matches
    ]
    return {
        "plugin": {
            "total": total,
            "swapped": swapped,
            "reference_added": reference_added,
            "rejected": rejected,
        },
        "normalization": normalization,
        "deduplication": deduplication,
    }


# ============================================================
# HELPER: Record variant-count QC
# ============================================================
def record_vcf_variant_counts(
    step,
    current,
    previous,
    chromosome,
    drop_warn_fraction: float,
    logger=None,
):

    if current is None:
        _emit(
            "%s: the variant count could not be read, so the drop check for "
            "this step could not run." % step,
            logger=logger, level="warn",
        )
        return

    if logger is not None and previous is None:
        logger.record("OBSERVED", "vcf_variant_count", stage=step, count=current)
    elif logger is not None:
        logger.qc(
            step,
            "Variant count after the %s stage of the bcftools pipeline." % step,
            previous,
            current,
            step="16 bcftools_annotate_liftover",
            warn=bool(previous) and current < previous * (1.0 - drop_warn_fraction),
        )
    else:
        _emit(
            "chr%s %s: %s variants" % (chromosome, step, current),
            logger=None,
        )


# ============================================================
# 2. ANNOTATION + LIFTOVER PIPELINE
# ============================================================
def annotate_and_liftover_vcf(
    output_dir: str,
    gwas_outputname: str,
    chromosome: str,
    external_eaf_file: str,
    default_dbsnp_file: str,
    genome_fasta_file: str,
    target_genome_fasta_file: str,
    gff_file: str,
    chain_file: str,
    grch_version: str,
    threads: int,
    output_layout: Mapping[str, str],
    executables: Mapping[str, str],
    vcf_config: Mapping[str, object],
    policies=None,
    logger=None,
    qc_info=None,
) -> str:

    plugin = str(vcf_config["liftover_plugin"])
    binaries = dict(executables)
    bash = binaries["bash"]
    bcftools = binaries["bcftools"]
    tabix = binaries["tabix"]

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    target_builds = dict(vcf_config["target_builds"])
    if grch_version not in target_builds:
        raise ValueError("No target genome build is configured for %s" % grch_version)
    target_build = str(target_builds[grch_version])
    values = {
        "dataset_id": gwas_outputname,
        "chromosome": chromosome,
        "build": grch_version,
        "target_build": target_build,
    }
    transcript_file = configured_output_path(
        outdir, output_layout["bcftools_transcript"], **values,
    )

    # ---- policy resolution -------------------------------------------------
    liftover_swap = _policy_value(policies, "vcf.liftover_swap")
    if liftover_swap not in LIFTOVER_SWAP_CHOICES:
        raise ValueError(
            "vcf.liftover_swap must be one of %s - got %r"
            % (", ".join(LIFTOVER_SWAP_CHOICES), liftover_swap)
        )
    warn_fraction = float(
        _policy_value(policies, "vcf.liftover_warn_fraction")
    )
    critical_fraction = float(
        _policy_value(policies, "vcf.liftover_critical_fraction")
    )
    fail_fraction = float(
        _policy_value(policies, "vcf.liftover_fail_fraction")
    )
    drop_warn_fraction = float(
        _policy_value(policies, "vcf.variant_drop_warn_fraction")
    )
    sort_memory_mb = int(
        _policy_value(policies, "vcf.sort_memory_mb_per_thread")
    )

    input_vcf = configured_output_path(
        outdir, output_layout["chromosome_raw_vcf"], **values,
    )

    if not input_vcf.exists():
        raise FileNotFoundError(f"[chr{chromosome}] Missing VCF: {input_vcf}")

    norm_vcf = configured_output_path(
        outdir, output_layout["chromosome_normalized_vcf"], **values,
    )
    id_vcf = configured_output_path(
        outdir, output_layout["chromosome_id_vcf"], **values,
    )
    af_vcf = configured_output_path(
        outdir, output_layout["chromosome_frequency_vcf"], **values,
    )
    csq_vcf = configured_output_path(
        outdir, output_layout["chromosome_annotated_vcf"], **values,
    )
    target_vcf = configured_output_path(
        outdir, output_layout["chromosome_lifted_vcf"], **values,
    )
    reject_vcf = configured_output_path(
        outdir, output_layout["chromosome_not_lifted_vcf"], **values,
    )

    counts = {}
    sort_tmp_dir = configured_output_path(
        outdir, output_layout["sort_temp_directory"], **values,
    )
    sort_tmp_dir.mkdir(parents=True, exist_ok=True)

    # The sort temp dir holds multi-GB spill files in a predictable, shared
    # location; it must be removed on the failure path too, not only on success.
    try:
        # =======================================================
        # ORIGINAL COUNT
        # =======================================================
        counts["ORIGINAL"] = get_vcf_variant_count(
            str(input_vcf), bcftools, logger=logger,
        )
        record_vcf_variant_counts(
            "ORIGINAL", counts["ORIGINAL"], None, chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        # =======================================================
        # STEP 0: NORMALIZATION
        # =======================================================
        norm_commands = [
            [bcftools, "norm", "--threads", str(threads), "-Ou", "-m-any",
             "-d", "exact", "--fasta-ref", str(genome_fasta_file), str(input_vcf)],
            [bcftools, "annotate", "--threads", str(threads), "-x", "ID",
             "-Oz", "-o", str(norm_vcf), "--write-index=tbi"],
        ]
        _run_bcftools_step(
            [bash, "-c", "set -euo pipefail\n" + " | ".join(
                shlex.join(command) for command in norm_commands
            )],
            transcript_file, "STEP0_NORM", chromosome, logger=logger,
        )

        counts["NORM"] = get_vcf_variant_count(
            str(norm_vcf), bcftools, logger=logger,
        )
        record_vcf_variant_counts(
            "NORM", counts["NORM"], counts["ORIGINAL"], chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        # =======================================================
        # STEP 1: ID ANNOTATION
        # =======================================================
        id_commands = [
            [bcftools, "annotate", "--threads", str(threads),
             "--annotations", str(default_dbsnp_file), "--columns", "ID",
             "--pair-logic", "exact", "-Ou", str(norm_vcf)],
            [bcftools, "annotate", "--threads", str(threads),
             "--set-id", str(vcf_config["missing_id_format"]), "-Oz", "-o",
             str(id_vcf), "--write-index=tbi"],
        ]
        _run_bcftools_step(
            [bash, "-c", "set -euo pipefail\n" + " | ".join(
                shlex.join(command) for command in id_commands
            )],
            transcript_file, "STEP1_ID", chromosome, logger=logger,
        )

        counts["ID"] = get_vcf_variant_count(
            str(id_vcf), bcftools, logger=logger,
        )
        record_vcf_variant_counts(
            "ID", counts["ID"], counts["NORM"], chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        # =======================================================
        # STEP 2: EAF ANNOTATION
        # =======================================================
        _run_bcftools_step([
            bcftools, "annotate",
            "--threads", str(threads),
            "--annotations", str(external_eaf_file),
            "--pair-logic", "exact",
            "--columns", ",".join(vcf_config["external_frequency_columns"]),
            "-Oz", "-o", str(af_vcf),
            "--write-index=tbi",
            str(id_vcf)
        ], transcript_file, "STEP2_AF", chromosome, logger=logger)

        counts["EAF"] = get_vcf_variant_count(
            str(af_vcf), bcftools, logger=logger,
        )
        record_vcf_variant_counts(
            "EAF", counts["EAF"], counts["ID"], chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        # =======================================================
        # STEP 3: CSQ
        # =======================================================
        _run_bcftools_step([
            bcftools, "csq",
            "--fasta-ref", str(genome_fasta_file),
            "--gff-annot", str(gff_file),
            "--threads", str(threads),
            "--unify-chr-names", "-,chr,-",
            "-Oz", "-o", str(csq_vcf),
            "--write-index=tbi",
            str(af_vcf)
        ], transcript_file, "STEP3_CSQ", chromosome, logger=logger)

        counts["CSQ"] = get_vcf_variant_count(
            str(csq_vcf), bcftools, logger=logger,
        )
        record_vcf_variant_counts(
            "CSQ", counts["CSQ"], counts["EAF"], chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        csq_header = read_vcf_header(
            csq_vcf, bcftools, logger=logger, error_type=RuntimeError,
        )
        plugin_tag_args, resolved_liftover_tags = _resolve_liftover_tag_arguments(
            csq_header, vcf_config,
        )
        if logger is not None:
            logger.record(
                "PARAM", "liftover_allele_dependent_fields",
                allele_frequency=resolved_liftover_tags["allele_frequency"],
                signed_effect=resolved_liftover_tags["signed_effect"],
            )

        # =======================================================
        # STEP 4: LIFTOVER
        # =======================================================
        # vcf.liftover_swap:
        #   exclude - count and discard successfully harmonised SWAP records.
        #   keep    - retain successfully harmonised SWAP records.
        # Tag harmonisation is mandatory in both modes.
        if logger is not None:
            logger.record(
                "ACTION", "liftover_allele_swap",
                policy=liftover_swap,
                outcome=("discard" if liftover_swap == "exclude" else "keep"),
            )

        sort_mem_mb = max(1, int(threads * sort_memory_mb))

        liftover_commands = [[
            bcftools, "+%s" % plugin, str(csq_vcf), "--no-version", "-Ou", "--",
            "--src-fasta-ref", str(genome_fasta_file),
            "--fasta-ref", str(target_genome_fasta_file),
            "--chain", str(chain_file), "--reject", str(reject_vcf),
            "--reject-type", "z", *plugin_tag_args,
        ]]
        if liftover_swap == "exclude":
            # INFO/SWAP is defined by the configured liftover plugin; these are
            # the plugin's two allele-swap states, not a tunable threshold.
            liftover_commands.append([
                bcftools, "view", "--threads", str(threads), "-e",
                "INFO/SWAP==1 || INFO/SWAP==-1",
            ])
        liftover_commands.extend([
            [bcftools, "norm", "--threads", str(threads), "-Ou", "-m-any",
             "--fasta-ref", str(target_genome_fasta_file)],
            [bcftools, "sort", "-m", "%dM" % sort_mem_mb, "--temp-dir",
             str(sort_tmp_dir), "-Ou"],
            [bcftools, "norm", "--threads", str(threads), "-d", "exact",
             "-Oz", "-o", str(target_vcf), "--write-index=tbi"],
        ])
        liftover_transcript = _run_bcftools_step(
            [bash, "-c", "set -euo pipefail\n" + " | ".join(
                shlex.join(command) for command in liftover_commands
            )],
            transcript_file, "STEP4_LIFTOVER", chromosome, logger=logger,
        )
        liftover_audit = _parse_liftover_audit(liftover_transcript)

        if reject_vcf.is_file() and reject_vcf.stat().st_size > 0:
            _run_bcftools_step(
                [tabix, "-f", "-p", "vcf", str(reject_vcf)],
                transcript_file, "STEP4_INDEX_REJECT", chromosome, logger=logger,
            )

        # =======================================================
        # POST-LIFTOVER QC
        # =======================================================
        counts["LIFTED"] = get_vcf_variant_count(
            str(target_vcf), bcftools, logger=logger,
        )

        reject_path = Path(reject_vcf)
        if not reject_path.exists() or reject_path.stat().st_size == 0:
            # A chromosome where everything lifted over leaves a 0-byte reject
            # file, which is not indexable; that is zero rejects, not "unknown".
            counts["REJECTED"] = 0
        else:
            counts["REJECTED"] = get_vcf_variant_count(
                str(reject_vcf), bcftools, logger=logger,
            )

        plugin_counts = liftover_audit["plugin"]
        normalization_counts = liftover_audit["normalization"]
        deduplication_counts = liftover_audit["deduplication"]
        swap_candidates = (
            plugin_counts["swapped"] + plugin_counts["reference_added"]
        )
        swap_excluded = swap_candidates if liftover_swap == "exclude" else 0
        expected_normalization_input = (
            plugin_counts["total"] - plugin_counts["rejected"] - swap_excluded
        )
        accounting_problems = []
        if expected_normalization_input < 0:
            accounting_problems.append(
                "plugin rejection and policy-exclusion counts exceed plugin input"
            )
        if counts["CSQ"] is not None and counts["CSQ"] != plugin_counts["total"]:
            accounting_problems.append(
                "CSQ=%s but the plugin processed %s"
                % (counts["CSQ"], plugin_counts["total"])
            )
        if (
            counts["REJECTED"] is not None
            and counts["REJECTED"] != plugin_counts["rejected"]
        ):
            accounting_problems.append(
                "reject VCF=%s but the plugin reported %s"
                % (counts["REJECTED"], plugin_counts["rejected"])
            )
        if normalization_counts["total"] != expected_normalization_input:
            accounting_problems.append(
                "target normalization received %s records, expected %s after "
                "plugin rejection and swap policy"
                % (normalization_counts["total"], expected_normalization_input)
            )
        expected_lifted = (
            deduplication_counts["total"]
            - deduplication_counts["duplicate_removed"]
        )
        if counts["LIFTED"] is not None and counts["LIFTED"] != expected_lifted:
            accounting_problems.append(
                "final target VCF=%s but sorted exact deduplication predicts %s"
                % (counts["LIFTED"], expected_lifted)
            )
        if accounting_problems:
            raise RuntimeError(
                "Liftover record accounting is inconsistent: %s."
                % "; ".join(accounting_problems)
            )

        if qc_info is not None:
            qc_info.update(
                input=plugin_counts["total"],
                rejected=plugin_counts["rejected"],
                swap_excluded=swap_excluded,
                swap_policy=liftover_swap,
                final=counts["LIFTED"],
            )
        if logger is not None:
            logger.record(
                "OBSERVED", "liftover_record_accounting",
                input=plugin_counts["total"],
                rejected=plugin_counts["rejected"],
                swapped=plugin_counts["swapped"],
                reference_added=plugin_counts["reference_added"],
                swap_excluded=swap_excluded,
                normalization_split=normalization_counts["split"],
                normalization_input=normalization_counts["total"],
                normalization_realigned=normalization_counts["realigned"],
                normalization_mismatch_removed=normalization_counts[
                    "mismatch_removed"
                ],
                normalization_skipped=normalization_counts["skipped"],
                deduplication_input=deduplication_counts["total"],
                duplicates_removed=deduplication_counts["duplicate_removed"],
                final=counts["LIFTED"],
            )
        if liftover_swap == "exclude" and swap_excluded:
            _emit(
                "Liftover excluded %s successfully harmonised SWAP records "
                "under vcf.liftover_swap='exclude' (swapped=%s, "
                "reference-added=%s)."
                % (
                    swap_excluded,
                    plugin_counts["swapped"],
                    plugin_counts["reference_added"],
                ),
                logger=logger,
                level="warn",
            )
        elif liftover_swap == "keep" and swap_candidates:
            _emit(
                "Liftover retained %s successfully harmonised SWAP records "
                "(swapped=%s, reference-added=%s)."
                % (
                    swap_candidates,
                    plugin_counts["swapped"],
                    plugin_counts["reference_added"],
                ),
                logger=logger,
            )

        record_vcf_variant_counts(
            "LIFTED", counts["LIFTED"], counts["CSQ"], chromosome,
            drop_warn_fraction=drop_warn_fraction, logger=logger,
        )

        if counts["REJECTED"]:
            _emit(
                "Liftover rejected %s variants." % counts["REJECTED"],
                logger=logger,
            )

        # A required target-build chromosome containing no records is never a
        # scientifically usable liftover result. This invariant is independent
        # of the configurable plugin-rejection limit because policy exclusion,
        # normalization and deduplication are separate audited transformations.
        if (
            counts["CSQ"] is not None
            and counts["CSQ"] > 0
            and counts["LIFTED"] == 0
        ):
            raise LiftoverFailureRateError(
                "chr%s: %s annotated variants entered liftover, but the required "
                "target-build VCF contains zero variants (reported rejected count: "
                "%r). A zero-survivor liftover is invalid regardless of "
                "vcf.liftover_fail_fraction."
                % (chromosome, counts["CSQ"], counts["REJECTED"])
            )

        # 🚨 LIFTOVER PLUGIN REJECTION WARNING
        if counts["CSQ"] and counts["REJECTED"] is not None:
            failure_rate = counts["REJECTED"] / counts["CSQ"]
            _emit(
                "Liftover plugin rejection rate: %.2f%%."
                % (failure_rate * 100),
                logger=logger,
            )

            if failure_rate > critical_fraction:
                message = (
                    "[chr%s] 🚨 WARNING: >%.0f%% variants were rejected by "
                    "the liftover plugin"
                    % (chromosome, critical_fraction * 100)
                )
                _emit(message, logger=logger, level="warn")
            elif failure_rate > warn_fraction:
                message = (
                    "[chr%s] ⚠️ Notice: >%.0f%% variants were rejected by "
                    "the liftover plugin"
                    % (chromosome, warn_fraction * 100)
                )
                _emit(message, logger=logger, level="warn")

            if failure_rate > fail_fraction:
                message = (
                    "chr%s: the liftover plugin rejected %.2f%% of variants, "
                    "which is above "
                    "vcf.liftover_fail_fraction (%.2f%%)."
                    % (chromosome, failure_rate * 100, fail_fraction * 100)
                )
                raise LiftoverFailureRateError(message)
        else:
            message = (
                "[chr%s] ⚠️ Liftover plugin rejection rate could not be computed "
                "(CSQ count %r, rejected count %r), so the failure-rate gate did not run."
                % (chromosome, counts.get("CSQ"), counts.get("REJECTED"))
            )
            _emit(message, logger=logger, level="warn")

        return str(target_vcf)

    finally:
        shutil.rmtree(sort_tmp_dir, ignore_errors=True)


# ============================================================
# 3. CONCAT FUNCTION
# ============================================================
def concat_vcfs_by_build(
    output_dir: str,
    gwas_outputname: str,
    grch_version: str,
    expected_chromosomes: Sequence[str],
    threads: int,
    output_layout: Mapping[str, str],
    executables: Mapping[str, str],
    vcf_config: Mapping[str, object],
    provenance: Mapping[str, Any] | None = None,
    policies=None,
    logger=None,
):
    binaries = dict(executables)
    bash = binaries["bash"]
    bcftools = binaries["bcftools"]
    tabix = binaries["tabix"]
    outdir = Path(output_dir)
    target_builds = dict(vcf_config["target_builds"])
    if grch_version not in target_builds:
        raise ValueError("No target genome build is configured for %s" % grch_version)
    target_build = str(target_builds[grch_version])

    def report(message, level="info"):
        _emit(message, logger=logger, level=level)

    # ---- policy resolution -------------------------------------------------
    min_valid_size = int(
        _policy_value(policies, "vcf.min_valid_size_bytes")
    )
    max_workers_policy = int(
        _policy_value(policies, "vcf.concat_max_workers")
    )
    concat_max_attempts = int(
        _policy_value(policies, "vcf.concat_max_attempts")
    )
    if concat_max_attempts < 1:
        raise ValueError("vcf.concat_max_attempts must be at least 1")
    on_merge_failure = _policy_value(policies, "vcf.on_merge_failure")
    if on_merge_failure not in ON_MERGE_FAILURE_CHOICES:
        raise ValueError(
            "vcf.on_merge_failure must be one of %s - got %r"
            % (", ".join(ON_MERGE_FAILURE_CHOICES), on_merge_failure)
        )

    total_cpus = max(1, int(threads))
    max_workers = max(1, min(max_workers_policy, total_cpus))
    threads_per_job = max(1, total_cpus // max_workers)

    # -------------------------------------------------------
    # 🔥 VALIDATION (AUTO FIX + SKIP BAD FILES)
    # -------------------------------------------------------
    def validate_and_fix_vcf(vcf, *, require_records=False):
        vcf = Path(vcf)

        if not vcf.is_file():
            report("Missing chromosome VCF: %s" % vcf, level="warn")
            return None, "file is missing"

        # skip tiny files (header-only / broken)
        if vcf.stat().st_size < min_valid_size:
            report(f"⚠️ Skipping empty VCF: {vcf}", level="warn")
            return None, "file is smaller than vcf.min_valid_size_bytes"

        tbi = Path(str(vcf) + ".tbi")

        # try to auto-index if missing
        if not tbi.is_file() or tbi.stat().st_size == 0:
            report(f"⚠️ Missing index → attempting to index: {vcf}", level="warn")
            try:
                run_checked_command(
                    [tabix, "-f", "-p", "vcf", str(vcf)],
                    "Indexing VCF before concatenation: %s" % vcf,
                    logger=logger,
                    error_type=RuntimeError,
                )
            except (RuntimeError, OSError) as exc:
                report(
                    "❌ Failed to index → skipping: %s (%s)" % (vcf, exc),
                    level="warn",
                )
                return None, "index creation failed: %s" % exc

        # final check
        if not tbi.is_file() or tbi.stat().st_size == 0:
            report(f"❌ Still no usable index → skipping: {vcf}", level="warn")
            return None, "index is missing or empty"

        try:
            count_text = run_checked_command(
                [bcftools, "index", "-n", str(vcf)],
                "Validating chromosome VCF and index: %s" % vcf,
                logger=logger,
                error_type=RuntimeError,
            ).strip()
            count = int(count_text)
            if count < 0:
                raise ValueError("negative record count")
            if require_records and count == 0:
                raise ValueError("zero records in a required VCF")
        except (RuntimeError, ValueError) as exc:
            report("❌ Invalid chromosome VCF or index → skipping: %s (%s)" % (vcf, exc),
                   level="warn")
            return None, str(exc)

        return vcf, None

    configured_chromosome_order = list(dict.fromkeys(
        str(value)
        for value in _policy_value(policies, "chromosome.allowed_after_split")
    ))
    allowed_chromosomes = set(configured_chromosome_order)
    chromosomes = list(dict.fromkeys(str(value) for value in expected_chromosomes))
    unexpected = [value for value in chromosomes if value not in allowed_chromosomes]
    if unexpected:
        raise ValueError(
            "Expected merge chromosomes are not allowed by "
            "chromosome.allowed_after_split: %s" % ", ".join(unexpected)
        )
    chromosome_rank = {
        chromosome: rank
        for rank, chromosome in enumerate(configured_chromosome_order)
    }
    chromosomes.sort(key=chromosome_rank.__getitem__)
    if logger is not None:
        logger.record(
            "PARAM", "vcf_concat_chromosome_order",
            chromosomes=chromosomes,
            source="chromosome.allowed_after_split",
        )

    def chromosome_outputs(pattern_name, *, build, target):
        paths = []
        for chromosome in chromosomes:
            path = configured_output_path(
                outdir,
                output_layout[pattern_name],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
                build=build,
                target_build=target,
            )
            if path.is_file():
                paths.append(path)
        return paths

    def required_chromosome_outputs(pattern_name, *, build, target):
        valid = []
        failures = []
        if not chromosomes:
            return valid, ["no completed chromosomes were supplied to the merge"]
        for chromosome in chromosomes:
            path = configured_output_path(
                outdir,
                output_layout[pattern_name],
                dataset_id=gwas_outputname,
                chromosome=chromosome,
                build=build,
                target_build=target,
            )
            if not path.is_file():
                failures.append("chromosome %s input is missing: %s" % (chromosome, path))
                continue
            checked, validation_failure = validate_and_fix_vcf(
                path, require_records=True,
            )
            if checked is None:
                failures.append(
                    "chromosome %s input is invalid: %s (%s)"
                    % (chromosome, path, validation_failure)
                )
            else:
                valid.append(checked)
        return valid, failures

    def optional_chromosome_outputs(pattern_name, *, build, target):
        valid = []
        for path in chromosome_outputs(pattern_name, build=build, target=target):
            checked, _failure = validate_and_fix_vcf(
                path, require_records=False,
            )
            if checked is not None:
                valid.append(checked)
        return valid

    # Build membership is carried by the resolved source/target build values
    # and rendered output-layout paths. A dataset ID containing a build name
    # therefore cannot move a chromosome VCF into the wrong merge group.
    source_annotated, source_failures = required_chromosome_outputs(
        "chromosome_annotated_vcf", build=grch_version, target=target_build,
    )
    target_annotated, target_failures = required_chromosome_outputs(
        "chromosome_lifted_vcf", build=grch_version, target=target_build,
    )
    build_files = {grch_version: source_annotated, target_build: target_annotated}
    notlifted_files = optional_chromosome_outputs(
        "chromosome_not_lifted_vcf", build=grch_version, target=target_build,
    )
    raw_files, raw_failures = required_chromosome_outputs(
        "chromosome_raw_vcf", build=grch_version, target=target_build,
    )

    def unlink_vcf(path):
        for candidate in (Path(path), Path(str(path) + ".tbi"), Path(str(path) + ".csi")):
            try:
                candidate.unlink()
            except FileNotFoundError:
                pass
            except OSError as exc:
                report(f"⚠️ Could not remove {candidate}: {exc}", level="warn")

    # -------------------------------------------------------
    # MERGE FUNCTION (SAFE)
    # -------------------------------------------------------
    header_template = str(vcf_config["genome_build_header"])

    def merge(
        tag, vcf_list, out, *, genome_build, merge_group, require_records,
    ):
        if not vcf_list:
            return None, "no valid chromosome inputs", 0

        try:
            build_header = header_template.format(build=genome_build)
        except (KeyError, ValueError) as exc:
            raise ValueError(
                "vcf_processing.genome_build_header cannot be rendered for %s: %s"
                % (genome_build, exc)
            ) from exc

        provenance_headers = _render_provenance_headers(
            provenance,
            vcf_config=vcf_config,
            input_build=grch_version,
            output_build=genome_build,
            merge_group=merge_group,
        )
        declared_headers = [build_header, *provenance_headers]
        annotate_command = [bcftools, "annotate"]
        for header_line in declared_headers:
            annotate_command.extend(["--header-line", header_line])
        annotate_command.extend([
            "--output-type", "z", "--output", str(out),
            "--threads", str(threads_per_job), "--write-index=tbi",
        ])

        index_path = Path(str(out) + ".tbi")
        failures = []
        for attempt in range(1, concat_max_attempts + 1):
            unlink_vcf(out)
            try:
                commands = [
                    [
                        bcftools, "concat", "-a", "--output-type", "u",
                        *map(str, vcf_list),
                    ],
                    annotate_command,
                ]
                run_checked_command(
                    [
                        bash, "-c", "set -euo pipefail\n" + " | ".join(
                            shlex.join(command) for command in commands
                        ),
                    ],
                    "Concatenating %s VCFs and declaring genome build %s "
                    "(attempt %d/%d)"
                    % (tag, genome_build, attempt, concat_max_attempts),
                    logger=logger,
                    error_type=RuntimeError,
                    expected_outputs=[out, index_path],
                )
                if Path(out).stat().st_size < min_valid_size:
                    raise RuntimeError(
                        "merged output is smaller than vcf.min_valid_size_bytes (%d): %s"
                        % (min_valid_size, out)
                    )
                header = run_checked_command(
                    [bcftools, "view", "--header-only", str(out)],
                    "Validating genome-build header for %s (attempt %d/%d)"
                    % (tag, attempt, concat_max_attempts),
                    logger=logger,
                    error_type=RuntimeError,
                )
                build_headers = [
                    line.strip()
                    for line in header.splitlines()
                    if line.lower().startswith("##genome_build=")
                ]
                if build_headers != [build_header]:
                    raise RuntimeError(
                        "merged VCF %s must contain exactly one %r header; found %r"
                        % (out, build_header, build_headers)
                    )
                for expected_header in provenance_headers:
                    matching_headers = [
                        line.strip()
                        for line in header.splitlines()
                        if line.strip() == expected_header
                    ]
                    if matching_headers != [expected_header]:
                        raise RuntimeError(
                            "merged VCF %s must contain exactly one provenance "
                            "header %r; found %d"
                            % (out, expected_header, len(matching_headers))
                        )
                if provenance_headers:
                    definitions = {
                        category: declared_vcf_tag_definitions(header, category)
                        for category in ("INFO", "FORMAT")
                    }
                    absent_output_fields = []
                    for role, field in dict(
                        vcf_config["provenance"]["output_fields"]
                    ).items():
                        category, tag_name = str(field).split("/", 1)
                        if tag_name not in definitions[category]:
                            absent_output_fields.append("%s=%s" % (role, field))
                    if absent_output_fields:
                        raise RuntimeError(
                            "merged VCF %s is missing configured provenance output "
                            "field(s): %s"
                            % (out, ", ".join(absent_output_fields))
                        )
                count_text = run_checked_command(
                    [bcftools, "index", "-n", str(out)],
                    "Validating merged VCF and index for %s (attempt %d/%d)"
                    % (tag, attempt, concat_max_attempts),
                    logger=logger,
                    error_type=RuntimeError,
                ).strip()
                try:
                    count = int(count_text)
                except ValueError as exc:
                    raise RuntimeError(
                        "merged VCF validation did not return a record count for %s: %r"
                        % (out, count_text)
                    ) from exc
                if count < 0:
                    raise RuntimeError("merged VCF returned a negative record count: %s" % out)
                if require_records and count == 0:
                    raise RuntimeError(
                        "required merged VCF contains zero records: %s" % out
                    )
                if attempt > 1:
                    report(
                        "✅ [%s] concat recovered on attempt %d/%d."
                        % (tag, attempt, concat_max_attempts)
                    )
                if logger is not None:
                    logger.record(
                        "OUTPUT", "vcf_genome_build_header",
                        build=genome_build, header=build_header, path=str(out),
                    )
                    logger.record(
                        "OUTPUT", "vcf_postgwas_provenance_headers",
                        path=str(out), merge_group=merge_group,
                        header_count=len(provenance_headers),
                    )
                return str(out), None, attempt
            except RuntimeError as exc:
                unlink_vcf(out)
                failures.append("attempt %d/%d: %s" % (
                    attempt, concat_max_attempts, exc,
                ))
                if attempt < concat_max_attempts:
                    report(
                        "⚠️ [%s] concat attempt %d/%d failed: %s. Retrying."
                        % (tag, attempt, concat_max_attempts, exc),
                        level="warn",
                    )
                    continue
                report(
                    "❌ [%s] concat of %d files failed after %d attempt(s): %s"
                    % (tag, len(vcf_list), concat_max_attempts, exc),
                    level=(
                        "error"
                        if require_records and on_merge_failure == "fail"
                        else "warn"
                    ),
                )
        return None, " | ".join(failures), concat_max_attempts

    # -------------------------------------------------------
    # PARALLEL EXECUTION
    # -------------------------------------------------------
    # tag -> (files that were merged, output path, records' coordinate build)
    merge_inputs = {
        build.lower(): (
            files,
            configured_output_path(
                outdir,
                output_layout["merged_build_vcf"],
                dataset_id=gwas_outputname,
                build=build,
                target_build=target_build,
            ),
            build,
        )
        for build, files in build_files.items()
    }
    merge_inputs.update({
        "gwas2vcf": (
            raw_files,
            configured_output_path(
                outdir, output_layout["merged_raw_vcf"],
                dataset_id=gwas_outputname,
                build=grch_version,
                target_build=target_build,
            ),
            grch_version,
        ),
        # Rejected records never entered the target coordinate system.
        "notlifted": (
            notlifted_files,
            configured_output_path(
                outdir, output_layout["merged_not_lifted_vcf"],
                dataset_id=gwas_outputname,
                build=grch_version,
                target_build=target_build,
            ),
            grch_version,
        ),
    })
    merge_labels = {
        **{build.lower(): build for build in build_files},
        "gwas2vcf": "gwas2vcf_%s" % grch_version,
        "notlifted": "notlifted_%s" % target_build,
    }
    group_keys = {
        "input_build": grch_version.lower(),
        "target_build": target_build.lower(),
        "raw_gwas2vcf": "gwas2vcf",
    }
    required_groups = {
        group_keys[name] for name in vcf_config["required_merge_groups"]
    }
    input_failures = {
        grch_version.lower(): source_failures,
        target_build.lower(): target_failures,
        "gwas2vcf": raw_failures,
        "notlifted": [],
    }

    # Remove stale merged outputs before deciding which tasks can run. A path
    # from an earlier attempt must never satisfy downstream output validation.
    for _files, out_path, _build in merge_inputs.values():
        unlink_vcf(out_path)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for key, (files, out_path, build) in merge_inputs.items():
            blocked = bool(input_failures[key]) and on_merge_failure == "fail"
            if files and not blocked:
                futures[key] = executor.submit(
                    merge,
                    merge_labels[key],
                    files,
                    out_path,
                    genome_build=build,
                    merge_group=key,
                    require_records=key in required_groups,
                )

    # -------------------------------------------------------
    # WAIT
    # -------------------------------------------------------
    results = dict((key, None) for key in merge_inputs)
    merge_attempts = dict((key, 0) for key in merge_inputs)
    failure_details = dict(
        (key, list(input_failures[key])) for key in merge_inputs
    )
    for key, future in futures.items():
        path, failure, attempts_used = future.result()
        results[key] = path
        merge_attempts[key] = attempts_used
        if failure:
            failure_details[key].append(failure)
    for key in required_groups:
        if not merge_inputs[key][0] and not failure_details[key]:
            failure_details[key].append("no valid chromosome inputs")

    merge_failures = [key for key in merge_inputs if failure_details[key]]
    required_merge_failures = [
        key for key in merge_failures if key in required_groups
    ]
    optional_merge_failures = [
        key for key in merge_failures if key not in required_groups
    ]

    # -------------------------------------------------------
    # CLEANUP
    # -------------------------------------------------------
    # Only the inputs of a merge that produced a real output file are deleted.
    # Previously six unquoted `os.system("rm ...")` globs ran unconditionally, so
    # a failed concat destroyed the entire product of the liftover pipeline while
    # the discarded wait status hid whether the removal had even worked.
    # Pure intermediates: nothing downstream reads them, so they always go.
    for pattern_name in (
        "chromosome_normalized_vcf", "chromosome_id_vcf",
        "chromosome_frequency_vcf",
    ):
        for path in chromosome_outputs(
            pattern_name, build=grch_version, target=target_build,
        ):
            unlink_vcf(path)

    for key, (files, _out, _build) in merge_inputs.items():
        if results.get(key) and key not in merge_failures:
            for path in files:
                unlink_vcf(path)
        elif files:
            outcome = (
                "a partial merged output was written but the required chromosome set "
                "was incomplete"
                if results.get(key)
                else "the merge did not produce a validated output file"
            )
            report(
                "🛟 [%s] keeping %d per-chromosome file(s): %s, so its inputs were "
                "NOT deleted." % (merge_labels[key], len(files), outcome),
                level="warn",
            )

    if required_merge_failures and on_merge_failure == "fail":
        detail = "; ".join(
            "%s: %s" % (merge_labels[key], " | ".join(failure_details[key]))
            for key in required_merge_failures
        )
        raise VcfMergeError(
            "Required VCF merge outputs failed validation. %s. Per-chromosome VCFs "
            "were kept so the merge can be retried. Set vcf.on_merge_failure to "
            "'continue' to carry on with explicitly recorded partial results." % detail
        )

    if required_merge_failures:
        report(
            "Required VCF merge output%s incomplete under "
            "vcf.on_merge_failure='continue': %s. Dataset status is PARTIAL."
            % (
                "s are" if len(required_merge_failures) != 1 else " is",
                ", ".join(merge_labels[key] for key in required_merge_failures),
            ),
            level="warn",
        )

    results["merge_failures"] = [merge_labels[k] for k in merge_failures]
    results["required_merge_failures"] = [
        merge_labels[key] for key in required_merge_failures
    ]
    results["optional_merge_failures"] = [
        merge_labels[key] for key in optional_merge_failures
    ]
    results["merge_failure_details"] = {
        merge_labels[key]: failure_details[key] for key in merge_failures
    }
    results["merge_attempts"] = {
        merge_labels[key]: merge_attempts[key] for key in merge_inputs
    }
    results["merge_status"] = "PARTIAL" if required_merge_failures else "OK"
    return results
