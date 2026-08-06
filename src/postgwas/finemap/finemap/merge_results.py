"""Select FINEMAP models and create scientifically valid FLAMES inputs."""

import logging
import re
import shutil
from pathlib import Path

import pandas as pd

from postgwas.finemap.defaults import (
    DEFAULT_COVERAGE_TOLERANCE,
    DEFAULT_CREDIBLE_SET_COVERAGE,
    DEFAULT_GENOME_BUILD,
)

logger = logging.getLogger("postgwas.finemap")
TARGET_CREDIBLE_SET_COVERAGE = DEFAULT_CREDIBLE_SET_COVERAGE
_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def parse_cred_header(file_path):
    """Read FINEMAP's posterior probability for the model represented by a .credK file."""
    try:
        with open(file_path, "r") as handle:
            for line_number, line in enumerate(handle):
                if "Post-Pr" in line:
                    match = re.search(rf"=\s*({_NUMBER})", line)
                    return float(match.group(1)) if match else 0.0
                if line_number >= 30 or not line.startswith("#"):
                    break
    except (OSError, ValueError):
        pass
    return 0.0


def reformat_snp_id(snp_id):
    """Backward-compatible formatter for an already encoded CHR_POS_A1_A2 ID."""
    if isinstance(snp_id, str):
        return snp_id.replace("_", ":", 2)
    return snp_id


def _model_size(filename):
    match = re.search(r"\.cred(\d+)$", filename)
    return int(match.group(1)) if match else 10**9


def select_and_copy_best_models(input_dir, output_dir, allowed_loci=None):
    """Select the most probable causal-count model for each successful locus."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    allowed = set(allowed_loci) if allowed_loci is not None else None
    all_files_tracking = []

    loci_folders = sorted(
        path for path in Path(input_dir).glob("*")
        if path.is_dir() and (allowed is None or path.name in allowed)
    )
    for folder in loci_folders:
        locus_name = folder.name
        logger.info("[STAGE] locus=%s stage=model_selection status=started", locus_name)
        for stale in output_path.glob(f"{locus_name}.cred*"):
            if stale.is_file():
                stale.unlink()
        stale_z = output_path / f"{locus_name}.z"
        if stale_z.is_file():
            stale_z.unlink()
        cred_files = sorted(
            path for path in folder.glob("*.cred*")
            if path.is_file() and re.search(r"\.cred\d+$", path.name)
        )
        if not cred_files:
            logger.warning(
                "[STAGE] locus=%s stage=model_selection status=skipped reason=no_credible_files",
                locus_name,
            )
            continue

        ranked = sorted(
            ((parse_cred_header(path), _model_size(path.name), path) for path in cred_files),
            key=lambda value: (-value[0], value[1], value[2].name),
        )
        best_file_path = ranked[0][2]
        logger.info(
            "[STAGE] locus=%s stage=model_selection status=completed model=%s posterior_probability=%.12g",
            locus_name,
            best_file_path.name,
            ranked[0][0],
        )
        for probability, k_value, path in ranked:
            record = {
                "Locus": locus_name,
                "File_Name": path.name,
                "K_Causal": k_value,
                "Posterior_Prob": probability,
                "Selected": path == best_file_path,
            }
            all_files_tracking.append(record)
            if path == best_file_path:
                shutil.copy2(path, output_path / path.name)
                z_source = folder / f"{locus_name}.z"
                if not z_source.is_file():
                    raise FileNotFoundError(
                        f"Cannot map FINEMAP credible-set IDs without {z_source}"
                    )
                shutil.copy2(z_source, output_path / z_source.name)
    return all_files_tracking


def save_summary_csv(tracking_data, output_dir):
    """Save all evaluated FINEMAP model probabilities and the selected model."""
    summary_path = Path(output_dir) / "finemap_all_models_summary.csv"
    if tracking_data:
        summary_df = pd.DataFrame(tracking_data).sort_values(
            by=["Locus", "Posterior_Prob", "K_Causal"],
            ascending=[True, False, True],
            kind="stable",
        )
    else:
        summary_df = pd.DataFrame(
            columns=["Locus", "File_Name", "K_Causal", "Posterior_Prob", "Selected"]
        )
    summary_df.to_csv(summary_path, index=False)


def filter_metadata(lines, k_val):
    """Retained for API compatibility; FLAMES exports no misleading FINEMAP comments."""
    target_idx = int(k_val) - 1
    filtered_lines = []
    for line in lines:
        if "Post-Pr" in line:
            filtered_lines.append(line)
            continue
        parts = line.split()
        value_idx = 1 + target_idx * 2
        na_idx = value_idx + 1
        filtered_lines.append(
            f"{parts[0]} {parts[value_idx]} {parts[na_idx]}\n"
            if len(parts) > na_idx else line
        )
    return filtered_lines


def _flames_chromosome(value):
    chrom = re.sub(r"^chr", "", str(value), flags=re.IGNORECASE)
    if chrom.upper() == "X":
        return "23"
    if chrom.upper() == "Y":
        return "24"
    return chrom


def _safe_filename(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "locus"


def create_flames_file(
    input_dir,
    output_dir,
    allowed_loci=None,
    locus_metadata=None,
    genome_build=DEFAULT_GENOME_BUILD,
    target_coverage=DEFAULT_CREDIBLE_SET_COVERAGE,
    coverage_tolerance=DEFAULT_COVERAGE_TOLERANCE,
):
    """Map FINEMAP rsIDs and emit one validated 95% credible set per file.

    ``target_coverage`` remains an argument for transparent configuration and
    testing, but values other than the scientifically required 0.95 are
    rejected.
    """
    if float(target_coverage) != DEFAULT_CREDIBLE_SET_COVERAGE:
        raise ValueError(
            "FINEMAP-to-FLAMES export requires target_coverage=0.95"
        )
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    annotation_dir = output_path / "annots"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    allowed = set(allowed_loci) if allowed_loci is not None else None
    metadata = locus_metadata or {}

    cred_files = sorted(
        path for path in input_path.glob("*.cred*")
        if path.is_file() and re.search(r"\.cred\d+$", path.name)
    )
    index_rows = []
    manifest_rows = []

    for file_path in cred_files:
        match = re.match(r"^(.*)\.cred(\d+)$", file_path.name)
        if not match:
            continue
        locus_name, model_k = match.group(1), int(match.group(2))
        if allowed is not None and locus_name not in allowed:
            continue
        logger.info(
            "[STAGE] locus=%s stage=flames_export status=started source=%s",
            locus_name,
            file_path.name,
        )
        z_file = input_path / f"{locus_name}.z"
        if not z_file.is_file():
            raise FileNotFoundError(f"Missing harmonised FINEMAP Z file: {z_file}")

        z_df = pd.read_csv(z_file, sep=r"\s+", dtype={"rsid": str, "chromosome": str})
        required_z = {"rsid", "chromosome", "position", "allele1", "allele2"}
        if not required_z.issubset(z_df.columns):
            raise ValueError(f"{z_file} lacks columns needed for FLAMES variant IDs")
        if z_df["rsid"].duplicated().any():
            raise ValueError(f"{z_file} contains duplicate rsIDs")
        variant_lookup = z_df.set_index("rsid")

        credible_df = pd.read_csv(file_path, comment="#", sep=r"\s+", dtype=str)
        credible_columns = sorted(
            (column for column in credible_df.columns if re.fullmatch(r"cred\d+", column)),
            key=lambda column: int(column.replace("cred", "")),
        )
        model_probability = parse_cred_header(file_path)
        for credible_column in credible_columns:
            cs_number = int(credible_column.replace("cred", ""))
            probability_column = f"prob{cs_number}"
            if probability_column not in credible_df.columns:
                raise ValueError(f"{file_path} has {credible_column} but no {probability_column}")

            subset = credible_df[[credible_column, probability_column]].copy()
            subset[probability_column] = pd.to_numeric(subset[probability_column], errors="coerce")
            subset = subset.dropna()
            subset = subset[
                ~subset[credible_column].astype(str).str.upper().isin({"NA", "NAN", "."})
            ]
            if subset.empty:
                continue
            if subset[credible_column].duplicated().any():
                raise ValueError(f"Duplicate variant in {file_path}:{credible_column}")
            probabilities = subset[probability_column]
            if ((probabilities < 0) | (probabilities > 1)).any():
                raise ValueError(f"Invalid posterior probability in {file_path}:{probability_column}")

            missing_ids = [rsid for rsid in subset[credible_column] if rsid not in variant_lookup.index]
            if missing_ids:
                preview = ", ".join(map(str, missing_ids[:5]))
                raise ValueError(
                    f"{file_path} contains credible-set IDs absent from its harmonised Z file: {preview}"
                )
            formatted_ids = []
            for rsid in subset[credible_column]:
                variant = variant_lookup.loc[rsid]
                formatted_ids.append(
                    f"{_flames_chromosome(variant['chromosome'])}:"
                    f"{int(variant['position'])}:"
                    f"{str(variant['allele1']).upper()}_{str(variant['allele2']).upper()}"
                )
            if len(set(formatted_ids)) != len(formatted_ids):
                raise ValueError(
                    f"Duplicate CHR:BP:A1_A2 identifier in {file_path}:{credible_column}"
                )

            flames_df = pd.DataFrame({
                "index": range(1, len(subset) + 1),
                "cred1": formatted_ids,
                "prob1": probabilities.astype(float).values,
            }).sort_values("prob1", ascending=False, kind="stable").reset_index(drop=True)
            flames_df["index"] = range(1, len(flames_df) + 1)
            achieved_coverage = float(flames_df["prob1"].sum())
            if achieved_coverage < float(target_coverage) - float(coverage_tolerance):
                raise ValueError(
                    f"{file_path}:{credible_column} has cumulative probability "
                    f"{achieved_coverage:.6g}, below the required 95%"
                )

            out_name = f"{file_path.name}_L{cs_number}.txt"
            out_file = output_path / out_name
            flames_df.to_csv(out_file, sep="\t", index=False, float_format="%.12g")

            locus_info = metadata.get(locus_name, {})
            genomic_locus = str(locus_info.get("GenomicLocus", locus_name))
            annotation_file = annotation_dir / f"annotated_{_safe_filename(locus_name)}_CS{cs_number}.txt"
            index_rows.append({
                "Filename": str(out_file.resolve()),
                "GenomicLocus": genomic_locus,
                "Annotfiles": str(annotation_file.resolve()),
            })
            manifest_rows.append({
                "Filename": out_name,
                "GenomicLocus": genomic_locus,
                "engine": "FINEMAP",
                "genome_build": genome_build,
                "model_k": model_k,
                "model_posterior_probability": model_probability,
                "credible_set": cs_number,
                "target_coverage": float(target_coverage),
                "achieved_coverage": achieved_coverage,
                "n_variants": len(flames_df),
                "allele1_definition": "effect_allele_aligned_to_LD_reference_A1",
                "credible_set_membership_definition": "finemap_95_percent_credible_set",
                "probability_definition": "finemap_snp_posterior_in_selected_model",
            })
            logger.info(
                "[STAGE] locus=%s stage=flames_export status=completed credible_set=%d "
                "variants=%d achieved_coverage=%.12g output=%s",
                locus_name,
                cs_number,
                len(flames_df),
                achieved_coverage,
                out_file,
            )

    if not index_rows:
        raise ValueError("No valid 95% FINEMAP credible sets were available for FLAMES export")
    index_file = output_path / "indexfile.txt"
    pd.DataFrame(index_rows).to_csv(index_file, sep="\t", index=False)
    pd.DataFrame(manifest_rows).to_csv(
        output_path / "finemap_FLAMES_manifest.tsv", sep="\t", index=False
    )
    logger.info(
        "[STAGE] stage=flames_index status=completed credible_sets=%d index=%s",
        len(index_rows),
        index_file,
    )


def process_finemap_output(
    raw_dir,
    inter_dir,
    final_dir,
    allowed_loci=None,
    locus_metadata=None,
    genome_build=DEFAULT_GENOME_BUILD,
    target_coverage=DEFAULT_CREDIBLE_SET_COVERAGE,
    coverage_tolerance=DEFAULT_COVERAGE_TOLERANCE,
):
    """Select best FINEMAP models and create validated FLAMES files and index."""
    tracking_data = select_and_copy_best_models(
        raw_dir, inter_dir, allowed_loci=allowed_loci
    )
    save_summary_csv(tracking_data, inter_dir)
    create_flames_file(
        inter_dir,
        final_dir,
        allowed_loci=allowed_loci,
        locus_metadata=locus_metadata,
        genome_build=genome_build,
        target_coverage=target_coverage,
        coverage_tolerance=coverage_tolerance,
    )
