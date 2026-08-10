"""Convert the structured harmonisation QC JSON to a tabular report."""

import json

import pandas as pd

from .policies import default_policies


def flatten_dict(values, parent_key="", separator="_"):
    """Flatten nested QC mappings while preserving their full key path."""
    flattened = {}
    for key, value in values.items():
        name = "%s%s%s" % (parent_key, separator, key) if parent_key else key
        if isinstance(value, dict):
            flattened.update(flatten_dict(value, name, separator))
        else:
            flattened[name] = value
    return flattened


def qc_results_to_dataframe(json_path=None, data=None, allowed_chromosomes=None):
    """Return chromosome metrics and global metrics in one DataFrame."""
    if data is None:
        if json_path is None:
            raise ValueError("Provide QC results as data or json_path")
        with open(json_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    if not isinstance(data, dict):
        raise TypeError("QC results must be a dictionary")

    rows = []
    global_metrics = {}
    configured_chromosomes = (
        allowed_chromosomes
        if allowed_chromosomes is not None
        else default_policies().get("chromosome.allowed")
    )
    allowed_chromosomes = {
        str(chromosome).upper() for chromosome in configured_chromosomes
    }
    for chromosome, content in data.items():
        if (
            str(chromosome).upper() not in allowed_chromosomes
            or not isinstance(content, dict)
        ):
            global_metrics[chromosome] = content
            continue
        for section, section_content in content.items():
            metrics = (
                flatten_dict(section_content)
                if isinstance(section_content, dict)
                else {section: section_content}
            )
            rows.extend(
                {
                    "section": section if isinstance(section_content, dict) else "summary",
                    "metric": metric,
                    "chromosome": "chr%s" % chromosome,
                    "value": value,
                }
                for metric, value in metrics.items()
            )

    if not rows:
        return pd.DataFrame()

    table = pd.DataFrame(rows).pivot(
        index=["section", "metric"], columns="chromosome", values="value"
    ).reset_index()
    if not global_metrics:
        return table

    global_table = pd.DataFrame(
        {
            "section": "GLOBAL",
            "metric": key,
            "GLOBAL_VALUE": value,
        }
        for key, value in global_metrics.items()
    )
    return pd.concat([table, global_table], ignore_index=True)
