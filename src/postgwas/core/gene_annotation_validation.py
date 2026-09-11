"""Read-only named gene/TSS annotation contracts for gene-score consumers."""

from copy import deepcopy

import numpy as np
import pandas as pd

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.io.tables import read_pandas_table, require_table_columns
from postgwas.core.paths import require_nonempty_file


def validate_gene_tss_annotation(
    path, *, delimiter, id_column, chromosome_column, tss_column, name_column,
    require_names, require_nonempty_identifiers, label, error_type=ValueError,
):
    """Validate a static annotation, leaving chromosome/anchor policy to callers.

    Preserve the existing pandas string-conversion and missing-name semantics;
    this extraction does not redefine identifier normalization. TSS positions
    are checked as finite/nonnegative without guessing the genome build.
    """
    path = require_nonempty_file(path, label, error_type=error_type)

    def inspect():
        table = read_pandas_table(path, delimiter, label, error_type=error_type)
        required = [id_column, chromosome_column, tss_column]
        if require_names:
            required.append(name_column)
        require_table_columns(table, required, label, error_type=error_type)
        identifiers = table[id_column].astype(str)
        chromosome = table[chromosome_column].astype(str)
        if identifiers.duplicated().any():
            raise error_type("%s contains duplicate gene identifiers." % label)
        if require_nonempty_identifiers and identifiers.str.strip().eq("").any():
            raise error_type("%s contains empty gene identifiers." % label)
        if chromosome.str.strip().eq("").any():
            raise error_type("%s contains empty chromosome labels." % label)
        tss = pd.to_numeric(table[tss_column], errors="coerce")
        if tss.isna().any() or not np.isfinite(tss.to_numpy()).all() or (tss < 0).any():
            raise error_type("%s TSS values must be finite and non-negative." % label)
        names = identifiers
        name_present = name_column in table.columns
        if name_present:
            configured_names = table[name_column].astype(str)
            if require_names and configured_names.str.strip().eq("").any():
                raise error_type("%s contains empty gene names." % label)
            names = configured_names.where(configured_names.str.strip().ne(""), identifiers)
        result = {
            "path": path, "genes": set(identifiers), "gene_count": len(identifiers),
            "chromosomes": sorted(set(chromosome)),
            "tss_minimum": float(tss.min()), "tss_maximum": float(tss.max()),
            "name_column_present": name_present,
            "gene_names": dict(zip(identifiers, names)),
            "gene_chromosomes": dict(zip(identifiers, chromosome)),
            "gene_tss": dict(zip(identifiers, tss)),
        }
        record_file_validation(
            path, "Gene TSS annotation",
            checks=("configured columns", "unique gene identifiers", "finite nonnegative TSS"),
            metrics={key: result[key] for key in ("gene_count", "chromosomes", "tss_minimum", "tss_maximum")},
            message="Genome build is a declaration, not inferred from these coordinates.",
        )
        return result

    return deepcopy(validate_once((path,), {
        "validator": "gene_tss_annotation", "delimiter": delimiter,
        "id_column": id_column, "chromosome_column": chromosome_column,
        "tss_column": tss_column, "name_column": name_column,
        "require_names": require_names,
        "require_nonempty_identifiers": require_nonempty_identifiers,
    }, inspect, error_type=error_type))
