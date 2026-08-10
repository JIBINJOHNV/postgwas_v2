"""Resolve overlapping primary fine-mapping results with one joint rerun.

Primary loci are always analysed before this module is called.  Connected
overlap groups are defined from the exact closed intervals used in that first
round.  Their union is rerun without applying the configured flank again.
"""

from __future__ import annotations

import copy
import logging
import shutil
from pathlib import Path

import pandas as pd

from postgwas.modules.fine_mapping.progress import write_run_configuration
from postgwas.modules.fine_mapping.output_layout import resolve_output_paths


logger = logging.getLogger("postgwas.modules.fine_mapping")


_PLAN_COLUMNS = [
    "overlap_group_id",
    "primary_genomic_locus",
    "chromosome",
    "primary_start",
    "primary_end",
    "joint_start",
    "joint_end",
    "joint_genomic_locus",
    "source_primary_genomic_loci",
    "resolution_status",
    "resolution_reason",
    "primary_analysis_status",
    "primary_warning_reason",
    "primary_failure_reason",
    "joint_analysis_status",
    "joint_warning_reason",
    "joint_failure_reason",
]

_PROVENANCE_COLUMNS = [
    "analysis_round",
    "final_selection_reason",
    "overlap_group_id",
    "final_genomic_locus",
    "source_primary_genomic_loci",
    "source_primary_credible_set_files",
    "source_primary_status_files",
    "source_primary_manifest_files",
    "source_primary_warning_reasons",
    "source_primary_failure_reasons",
    "source_result_file",
    "final_credible_set_file",
    "warning_reason",
    "failure_reason",
]


def _read_loci(path: str | Path) -> pd.DataFrame:
    loci = pd.read_csv(path, sep="\t", dtype=str)
    required = {"GenomicLocus", "chr", "start", "end"}
    missing = required - set(loci.columns)
    if missing:
        raise ValueError(
            "Fine-mapping genomic-loci output is missing columns: "
            + ", ".join(sorted(missing))
        )
    return _normalise_loci(loci)


def build_connected_overlap_plan(
    loci: pd.DataFrame,
    maximum_joint_region_kb: int | None,
    overlap_group_prefix: str,
    provenance_delimiter: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return a locus audit and one joint interval per connected overlap group."""
    records = []
    joint_records = []
    group_number = 0

    for chromosome, chromosome_loci in loci.groupby("chr", sort=False):
        connected = []
        connected_end = None
        for row in chromosome_loci.itertuples(index=False):
            if connected and int(row.start) > int(connected_end):
                group_number = _finish_group(
                    connected,
                    chromosome,
                    group_number,
                    maximum_joint_region_kb,
                    overlap_group_prefix,
                    provenance_delimiter,
                    records,
                    joint_records,
                )
                connected = []
                connected_end = None
            connected.append(row)
            connected_end = max(int(row.end), int(connected_end or row.end))
        if connected:
            group_number = _finish_group(
                connected,
                chromosome,
                group_number,
                maximum_joint_region_kb,
                overlap_group_prefix,
                provenance_delimiter,
                records,
                joint_records,
            )

    plan = pd.DataFrame(records, columns=_PLAN_COLUMNS)
    joints = pd.DataFrame(
        joint_records,
        columns=[
            "overlap_group_id", "CHR", "START", "END", "LP",
            "GenomicLocus", "source_primary_genomic_loci", "eligible_for_rerun",
            "resolution_reason",
        ],
    )
    return plan, joints


def _finish_group(
    rows,
    chromosome,
    group_number,
    maximum_joint_region_kb,
    overlap_group_prefix,
    provenance_delimiter,
    records,
    joint_records,
):
    if len(rows) == 1:
        row = rows[0]
        records.append({
            "overlap_group_id": "",
            "primary_genomic_locus": row.GenomicLocus,
            "chromosome": chromosome,
            "primary_start": int(row.start),
            "primary_end": int(row.end),
            "joint_start": pd.NA,
            "joint_end": pd.NA,
            "joint_genomic_locus": "",
            "source_primary_genomic_loci": row.GenomicLocus,
            "resolution_status": "primary_retained",
            "resolution_reason": "no_primary_result_overlap",
            "primary_analysis_status": "",
            "primary_warning_reason": "",
            "primary_failure_reason": "",
            "joint_analysis_status": "not_required",
            "joint_warning_reason": "",
            "joint_failure_reason": "",
        })
        return group_number

    group_number += 1
    group_id = f"{overlap_group_prefix}{group_number:04d}"
    joint_start = min(int(row.start) for row in rows)
    joint_end = max(int(row.end) for row in rows)
    source_loci = provenance_delimiter.join(str(row.GenomicLocus) for row in rows)
    joint_locus = f"chr{chromosome}:{joint_start}-{joint_end}"
    span_kb = (joint_end - joint_start + 1) / 1000.0
    eligible = maximum_joint_region_kb is None or span_kb <= maximum_joint_region_kb
    reason = (
        "connected_primary_boundaries_overlap"
        if eligible
        else "joint_region_exceeds_configured_maximum"
    )
    for row in rows:
        records.append({
            "overlap_group_id": group_id,
            "primary_genomic_locus": row.GenomicLocus,
            "chromosome": chromosome,
            "primary_start": int(row.start),
            "primary_end": int(row.end),
            "joint_start": joint_start,
            "joint_end": joint_end,
            "joint_genomic_locus": joint_locus,
            "source_primary_genomic_loci": source_loci,
            "resolution_status": (
                "joint_rerun_requested" if eligible else "excluded_before_joint_rerun"
            ),
            "resolution_reason": reason,
            "primary_analysis_status": "",
            "primary_warning_reason": "",
            "primary_failure_reason": "",
            "joint_analysis_status": "pending" if eligible else "not_run",
            "joint_warning_reason": "",
            "joint_failure_reason": "",
        })
    joint_records.append({
        "overlap_group_id": group_id,
        "CHR": chromosome,
        "START": joint_start,
        "END": joint_end,
        "LP": pd.NA,
        "GenomicLocus": joint_locus,
        "source_primary_genomic_loci": source_loci,
        "eligible_for_rerun": eligible,
        "resolution_reason": reason,
    })
    return group_number


def _read_index(result: dict) -> pd.DataFrame:
    columns = ["Filename", "GenomicLocus", "Annotfiles"]
    index_value = result.get("flames_index")
    if not index_value:
        return pd.DataFrame(columns=columns)
    index_file = Path(index_value)
    if not index_file.is_file():
        return pd.DataFrame(columns=columns)
    index = pd.read_csv(index_file, sep="\t", dtype=str)
    missing = set(columns) - set(index.columns)
    if missing:
        raise ValueError(
            "Fine-mapping FLAMES index is missing columns: "
            + ", ".join(sorted(missing))
        )
    index = index[columns].copy()
    index["Filename"] = index["Filename"].map(
        lambda value: str(
            (Path(value) if Path(value).is_absolute() else index_file.parent / value)
            .resolve()
        )
    )
    return index


def _read_manifests(result: dict) -> pd.DataFrame:
    frames = []
    for manifest in result.get("credible_set_manifests", []):
        path = Path(manifest)
        if not path.is_file() or path.stat().st_size == 0:
            continue
        frame = pd.read_csv(path, sep="\t", dtype=str)
        if {"Filename", "GenomicLocus"}.issubset(frame.columns):
            frame["source_result_file"] = str(path.resolve())
            frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def _status_lookup(result: dict) -> dict[str, dict]:
    status_path = result.get("locus_status")
    if not status_path or not Path(status_path).is_file():
        return {}
    status = pd.read_csv(status_path, sep="\t", dtype=str)
    key = next(
        (column for column in ("genomic_locus", "GenomicLocus") if column in status),
        None,
    )
    if key is None:
        return {}
    return {
        str(row[key]): row.to_dict()
        for _, row in status.drop_duplicates(key, keep="last").iterrows()
    }


def _read_primary_result_loci(result: dict) -> pd.DataFrame:
    """Read every successful primary locus, including converged no-set SuSiE fits."""
    status_path = result.get("locus_status")
    if status_path and Path(status_path).is_file():
        status = pd.read_csv(status_path, sep="\t", dtype=str)
        susie_columns = {
            "genomic_locus", "locus_chr", "locus_start", "locus_end", "converged",
        }
        finemap_columns = {
            "genomic_locus", "chromosome", "start", "end", "status",
        }
        if susie_columns.issubset(status.columns):
            successful = status[
                status["converged"].astype(str).str.lower().eq("true")
            ]
            loci = successful.rename(columns={
                "genomic_locus": "GenomicLocus",
                "locus_chr": "chr",
                "locus_start": "start",
                "locus_end": "end",
            })
            return _normalise_loci(loci)
        if finemap_columns.issubset(status.columns):
            successful = status[status["status"].astype(str).eq("success")]
            loci = successful.rename(columns={
                "genomic_locus": "GenomicLocus",
                "chromosome": "chr",
            })
            return _normalise_loci(loci)
    return _read_loci(result["genomic_loci"])


def _normalise_loci(loci: pd.DataFrame) -> pd.DataFrame:
    required = ["GenomicLocus", "chr", "start", "end"]
    temporary = loci[required].drop_duplicates().copy()
    temporary["chr"] = temporary["chr"].astype(str).str.replace(
        r"^chr", "", regex=True, case=False
    )
    temporary["start"] = pd.to_numeric(
        temporary["start"], errors="raise"
    ).astype("int64")
    temporary["end"] = pd.to_numeric(
        temporary["end"], errors="raise"
    ).astype("int64")
    if (temporary["start"] < 1).any() or (
        temporary["end"] < temporary["start"]
    ).any():
        raise ValueError("Fine-mapping locus status has invalid boundaries")
    if temporary["GenomicLocus"].duplicated().any():
        raise ValueError("Fine-mapping locus status has duplicate locus identifiers")
    return temporary.sort_values(
        ["chr", "start", "end"], kind="stable"
    ).reset_index(drop=True)


def _manifest_row(manifests, filename, genomic_locus):
    if manifests.empty:
        return {}
    basename = Path(filename).name
    matched = manifests[
        manifests["Filename"].map(lambda value: Path(str(value)).name).eq(basename)
        & manifests["GenomicLocus"].astype(str).eq(str(genomic_locus))
    ]
    return matched.iloc[0].to_dict() if not matched.empty else {}


def _write_final_outputs(
    settings,
    output_paths,
    primary_result,
    joint_result,
    plan,
):
    final_dir = output_paths["downstream_flames_directory"]
    if final_dir.exists():
        shutil.rmtree(final_dir)
    annotation_dir = output_paths["downstream_flames_annotations_directory"]

    primary_index = _read_index(primary_result)
    joint_index = _read_index(joint_result or {})
    overlapping = set(
        plan.loc[plan["overlap_group_id"].ne(""), "primary_genomic_locus"].astype(str)
    )
    retained_primary = primary_index[
        ~primary_index["GenomicLocus"].astype(str).isin(overlapping)
    ].copy()
    selected = [("primary", retained_primary, primary_result)]
    if joint_result:
        selected.append(("joint", joint_index, joint_result))

    primary_files_by_locus = {
        locus: settings["provenance_delimiter"].join(group["Filename"].astype(str))
        for locus, group in primary_index.groupby("GenomicLocus", sort=False)
    }
    primary_manifest_files = settings["provenance_delimiter"].join(
        str(Path(path).resolve())
        for path in primary_result.get("credible_set_manifests", [])
    )
    group_by_joint = {
        str(row.joint_genomic_locus): row
        for row in plan.loc[plan["overlap_group_id"].ne("")]
        .drop_duplicates("joint_genomic_locus")
        .itertuples(index=False)
    }
    primary_manifests = _read_manifests(primary_result)
    joint_manifests = _read_manifests(joint_result or {})
    primary_status = _status_lookup(primary_result)
    joint_status = _status_lookup(joint_result or {})
    index_rows = []
    combined_rows = []

    for round_name, index, result in selected:
        manifests = primary_manifests if round_name == "primary" else joint_manifests
        statuses = primary_status if round_name == "primary" else joint_status
        prefix = settings[f"{round_name}_file_prefix"]
        for sequence, row in enumerate(index.itertuples(index=False), start=1):
            source = Path(row.Filename)
            if not source.is_absolute():
                source = Path(result["flames_input"]) / source
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(f"Selected credible-set file is missing: {source}")
            final_dir.mkdir(parents=True, exist_ok=True)
            destination = final_dir / f"{prefix}{sequence:05d}_{source.name}"
            shutil.copy2(source, destination)
            annotation = annotation_dir / (
                f"{settings['annotation_prefix']}{destination.stem}.txt"
            )
            index_rows.append({
                "Filename": str(destination.resolve()),
                "GenomicLocus": row.GenomicLocus,
                "Annotfiles": str(annotation.resolve()),
            })

            metadata = _manifest_row(manifests, row.Filename, row.GenomicLocus)
            status = statuses.get(str(row.GenomicLocus), {})
            if round_name == "primary":
                overlap_group = ""
                source_loci = str(row.GenomicLocus)
                source_files = str(source.resolve())
                selection_reason = "non_overlapping_primary_result"
            else:
                group = group_by_joint.get(str(row.GenomicLocus))
                overlap_group = group.overlap_group_id if group else ""
                source_loci = group.source_primary_genomic_loci if group else ""
                source_files = settings["provenance_delimiter"].join(
                    value
                    for locus in source_loci.split(settings["provenance_delimiter"])
                    for value in [primary_files_by_locus.get(locus, "")]
                    if value
                )
                selection_reason = "joint_rerun_replaces_overlapping_primary_results"
            source_statuses = [
                primary_status.get(locus, {})
                for locus in source_loci.split(settings["provenance_delimiter"])
                if locus
            ]
            source_warnings = settings["provenance_delimiter"].join(
                f"{locus}={value}"
                for locus, source_status in zip(
                    source_loci.split(settings["provenance_delimiter"]),
                    source_statuses,
                )
                for value in [source_status.get("warning_reason")]
                if pd.notna(value) and str(value)
            )
            source_failures = settings["provenance_delimiter"].join(
                f"{locus}={value}"
                for locus, source_status in zip(
                    source_loci.split(settings["provenance_delimiter"]),
                    source_statuses,
                )
                for value in [source_status.get("failure_reason")]
                if pd.notna(value) and str(value)
            )
            metadata.update({
                "analysis_round": round_name,
                "final_selection_reason": selection_reason,
                "overlap_group_id": overlap_group,
                "final_genomic_locus": row.GenomicLocus,
                "source_primary_genomic_loci": source_loci,
                "source_primary_credible_set_files": source_files,
                "source_primary_status_files": str(
                    Path(primary_result["locus_status"]).resolve()
                ),
                "source_primary_manifest_files": primary_manifest_files,
                "source_primary_warning_reasons": source_warnings,
                "source_primary_failure_reasons": source_failures,
                "source_result_file": metadata.get(
                    "source_result_file", str(source.resolve())
                ),
                "final_credible_set_file": str(destination.resolve()),
                "warning_reason": status.get(
                    "warning_reason", metadata.get("warning_reason")
                ),
                "failure_reason": status.get(
                    "failure_reason", metadata.get("failure_reason")
                ),
            })
            combined_rows.append(metadata)

    index_frame = pd.DataFrame(
        index_rows, columns=["Filename", "GenomicLocus", "Annotfiles"]
    )
    if not index_frame.empty:
        index_frame.to_csv(
            final_dir / settings["index_filename"], sep="\t", index=False
        )
    combined = pd.DataFrame(combined_rows)
    for column in _PROVENANCE_COLUMNS:
        if column not in combined:
            combined[column] = pd.Series(dtype="object")
    combined_path = (
        output_paths["combined_results_directory"]
        / settings["final_combined_filename"]
    )
    combined_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(combined_path, sep="\t", index=False)

    loci_columns = [
        "GenomicLocus", "chr", "start", "end", "analysis_round",
        "overlap_group_id", "source_primary_genomic_loci",
    ]
    locus_rows = []
    for row in retained_primary.drop_duplicates("GenomicLocus").itertuples(index=False):
        source = plan.loc[plan["primary_genomic_locus"].eq(row.GenomicLocus)].iloc[0]
        locus_rows.append({
            "GenomicLocus": row.GenomicLocus,
            "chr": source.chromosome,
            "start": source.primary_start,
            "end": source.primary_end,
            "analysis_round": "primary",
            "overlap_group_id": "",
            "source_primary_genomic_loci": row.GenomicLocus,
        })
    for locus in joint_index["GenomicLocus"].drop_duplicates():
        group = group_by_joint.get(str(locus))
        if group:
            locus_rows.append({
                "GenomicLocus": locus,
                "chr": group.chromosome,
                "start": group.joint_start,
                "end": group.joint_end,
                "analysis_round": "joint",
                "overlap_group_id": group.overlap_group_id,
                "source_primary_genomic_loci": group.source_primary_genomic_loci,
            })
    if not index_frame.empty:
        pd.DataFrame(locus_rows, columns=loci_columns).to_csv(
            final_dir / settings["genomic_loci_filename"], sep="\t", index=False
        )
    return final_dir, combined


def resolve_overlapping_results(args, primary_result: dict, run_engine):
    """Run connected overlap groups once and build the authoritative handoff."""
    settings = dict(args.overlap_resolution)
    root = Path(primary_result["output_dir"]).resolve()
    output_paths = resolve_output_paths(
        root,
        args.fine_mapping_output_layout,
        args.dataset_id,
    )
    for key in ("quality_control_directory", "run_metadata_directory"):
        output_paths[key].mkdir(parents=True, exist_ok=True)
    protected_primary_directories = {
        Path(value).resolve()
        for value in [
            primary_result.get("flames_input"),
            Path(primary_result["flames_index"]).parent
            if primary_result.get("flames_index") else None,
            Path(primary_result["genomic_loci"]).parent
            if primary_result.get("genomic_loci") else None,
        ]
        if value
    }
    configured_directories = {
        "joint_rerun_directory": output_paths["joint_rerun_directory"],
        "downstream_flames_directory": output_paths[
            "downstream_flames_directory"
        ],
    }
    collisions = {
        name: path
        for name, path in configured_directories.items()
        if path in protected_primary_directories
    }
    if collisions:
        detail = ", ".join(
            f"{name}={path}" for name, path in sorted(collisions.items())
        )
        raise ValueError(
            "Overlap-resolution output directories cannot replace primary "
            f"fine-mapping artifacts: {detail}"
        )
    configuration_path = (
        output_paths["run_metadata_directory"]
        / settings["configuration_filename"]
    )
    write_run_configuration(
        configuration_path,
        {
            "policy": settings["policy"],
            "joint_failure_policy": settings["joint_failure_policy"],
            "settings": settings,
            "primary_output_directory": primary_result["output_dir"],
            "primary_locus_status": primary_result.get("locus_status"),
        },
    )
    loci = _read_primary_result_loci(primary_result)
    if loci.empty:
        logger.info("Overlap resolution skipped: primary analysis had no successful loci")
        return primary_result
    plan, joints = build_connected_overlap_plan(
        loci,
        settings["maximum_joint_region_kb"],
        settings["overlap_group_prefix"],
        settings["provenance_delimiter"],
    )
    primary_status = _status_lookup(primary_result)
    for index, row in plan.iterrows():
        status = primary_status.get(str(row["primary_genomic_locus"]), {})
        converged = str(status.get("converged", "")).lower()
        plan.at[index, "primary_analysis_status"] = status.get(
            "status",
            "converged" if converged == "true" else status.get("stage", "success"),
        )
        plan.at[index, "primary_warning_reason"] = status.get("warning_reason", "")
        plan.at[index, "primary_failure_reason"] = status.get("failure_reason", "")
    plan_path = output_paths["quality_control_directory"] / settings["plan_filename"]
    joint_locus_path = (
        output_paths["prepared_loci_directory"]
        / settings["joint_locus_filename"]
    )
    eligible = joints[joints["eligible_for_rerun"].eq(True)].copy()
    joint_result = None

    if not eligible.empty:
        eligible["LP"] = float(args.lp_threshold)
        joint_locus_path.parent.mkdir(parents=True, exist_ok=True)
        eligible[["CHR", "START", "END", "LP", "GenomicLocus"]].to_csv(
            joint_locus_path, sep="\t", index=False
        )
        joint_args = copy.copy(args)
        joint_args.locus_file = str(joint_locus_path)
        joint_args.locus_type = "range"
        joint_args.window_kb = 0
        joint_args.output_directory = str(output_paths["joint_rerun_directory"])
        if primary_result.get("summary_statistics_chromosome_manifest"):
            joint_args.summary_statistics_chromosome_manifest = primary_result[
                "summary_statistics_chromosome_manifest"
            ]
        logger.info(
            "Overlap resolution: rerunning %d connected group(s) without a second flank",
            len(eligible),
        )
        try:
            joint_result = run_engine(joint_args)
        except Exception as exc:
            logger.exception("Joint overlap rerun failed: %s", exc)
            plan.loc[
                plan["resolution_status"].eq("joint_rerun_requested"),
                ["resolution_status", "resolution_reason"],
            ] = ["joint_rerun_failed_excluded", str(exc)]
            plan.loc[
                plan["joint_analysis_status"].eq("pending"),
                ["joint_analysis_status", "joint_failure_reason"],
            ] = ["failed", str(exc)]
            joint_output_paths = resolve_output_paths(
                joint_args.output_directory,
                joint_args.fine_mapping_output_layout,
                joint_args.dataset_id,
            )
            status_key = (
                "susie_qc_file"
                if joint_args.finemap_method == "susie"
                else "finemap_locus_status_file"
            )
            partial_status_path = joint_output_paths[status_key]
            partial_joint_status = _status_lookup(
                {"locus_status": str(partial_status_path)}
                if partial_status_path and partial_status_path.is_file()
                else {}
            )
            for index, row in plan.loc[
                plan["joint_analysis_status"].eq("failed")
            ].iterrows():
                status = partial_joint_status.get(
                    str(row["joint_genomic_locus"]), {}
                )
                if not status:
                    continue
                failure_reason = status.get("failure_reason", str(exc))
                plan.at[index, "resolution_reason"] = failure_reason
                plan.at[index, "joint_warning_reason"] = status.get(
                    "warning_reason", ""
                )
                plan.at[index, "joint_failure_reason"] = failure_reason

    successful_joint_loci = set()
    if joint_result:
        if joint_result.get("genomic_loci") and Path(
            joint_result["genomic_loci"]
        ).is_file():
            successful_joint_loci = set(
                _read_loci(joint_result["genomic_loci"])["GenomicLocus"].astype(str)
            )
        joint_status = _status_lookup(joint_result)
        for index, row in plan.loc[plan["overlap_group_id"].ne("")].iterrows():
            status = joint_status.get(str(row["joint_genomic_locus"]), {})
            converged = str(status.get("converged", "")).lower()
            plan.at[index, "joint_analysis_status"] = status.get(
                "status",
                "converged" if converged == "true" else status.get("stage", "unknown"),
            )
            plan.at[index, "joint_warning_reason"] = status.get("warning_reason", "")
            plan.at[index, "joint_failure_reason"] = status.get("failure_reason", "")
    requested = plan["resolution_status"].eq("joint_rerun_requested")
    plan.loc[
        requested & plan["joint_genomic_locus"].isin(successful_joint_loci),
        ["resolution_status", "resolution_reason"],
    ] = ["replaced_by_joint_rerun", "joint_rerun_produced_credible_set"]
    plan.loc[
        requested & ~plan["joint_genomic_locus"].isin(successful_joint_loci),
        ["resolution_status", "resolution_reason"],
    ] = ["joint_rerun_no_credible_set_excluded", "joint_rerun_produced_no_credible_set"]
    plan.to_csv(plan_path, sep="\t", index=False)

    final_dir, combined = _write_final_outputs(
        settings, output_paths, primary_result, joint_result, plan
    )
    n_overlap_groups = int(
        plan.loc[
            plan["overlap_group_id"].ne(""), "overlap_group_id"
        ].nunique()
    )
    n_joint_success = int(
        plan.loc[
            plan["resolution_status"].eq("replaced_by_joint_rerun"),
            "overlap_group_id",
        ].nunique()
    )
    n_joint_no_sets = int(
        plan.loc[
            plan["resolution_status"].eq("joint_rerun_no_credible_set_excluded"),
            "overlap_group_id",
        ].nunique()
    )
    n_joint_failed = int(
        plan.loc[
            plan["resolution_status"].isin(
                ["joint_rerun_failed_excluded", "excluded_before_joint_rerun"]
            ),
            "overlap_group_id",
        ].nunique()
    )
    logger.info("Fine-mapping overlap-resolution final summary")
    logger.info("  successful primary loci: %d", len(loci))
    logger.info("  connected overlap groups: %d", n_overlap_groups)
    logger.info("  successful joint reruns: %d", n_joint_success)
    logger.info("  joint reruns without credible sets: %d", n_joint_no_sets)
    logger.info("  failed or size-excluded joint reruns: %d", n_joint_failed)
    logger.info("  final credible sets: %d", len(combined))
    logger.info("  overlap audit TSV: %s", plan_path)
    final_combined_path = (
        output_paths["combined_results_directory"]
        / settings["final_combined_filename"]
    )
    logger.info("  final combined TSV: %s", final_combined_path)

    if combined.empty and n_joint_failed:
        final_status = "failed_no_authoritative_credible_sets"
    elif combined.empty:
        final_status = "completed_no_credible_sets"
    elif n_joint_failed or n_joint_no_sets:
        final_status = "completed_with_excluded_overlap_groups"
    else:
        final_status = "success"
    logger.info("  overall status: %s", final_status)
    summary_path = (
        output_paths["quality_control_directory"] / settings["summary_filename"]
    )
    pd.DataFrame([{
        "status": final_status,
        "n_successful_primary_loci": len(loci),
        "n_overlap_groups": n_overlap_groups,
        "n_joint_rerun_successful": n_joint_success,
        "n_joint_rerun_without_credible_sets": n_joint_no_sets,
        "n_joint_rerun_failed_or_size_excluded": n_joint_failed,
        "n_final_credible_sets": len(combined),
        "overlap_audit_tsv": str(plan_path),
        "final_combined_tsv": str(final_combined_path),
    }]).to_csv(summary_path, sep="\t", index=False)
    logger.info("  overlap summary TSV: %s", summary_path)
    if final_status == "failed_no_authoritative_credible_sets":
        raise RuntimeError(
            "Overlap resolution produced no authoritative credible sets; inspect "
            f"{plan_path} and {summary_path}"
        )
    resolved = dict(primary_result)
    resolved.update({
        "status": final_status,
        "flames_input": str(final_dir) if not combined.empty else None,
        "final_combined_credible_sets": str(
            final_combined_path
        ),
        "overlap_resolution": str(plan_path),
        "overlap_resolution_summary": str(summary_path),
        "overlap_resolution_configuration": str(configuration_path),
        "n_overlap_groups": n_overlap_groups,
        "n_joint_rerun_successful": n_joint_success,
        "n_joint_rerun_failed": n_joint_failed,
        "n_joint_rerun_without_credible_sets": n_joint_no_sets,
        "n_final_credible_sets": int(len(combined)),
    })
    return resolved
