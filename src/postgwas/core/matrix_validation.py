"""Read-only name and matrix file contracts shared by resource consumers.

The checks do not scale, reorder, impute, or filter data. Numeric tolerances,
binary layout and chunk sizes are supplied by the consuming configuration.
Only compact evidence and identifier vectors are cached, never full matrices.
"""

from pathlib import Path

import numpy as np

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.paths import require_nonempty_file


def read_unique_names(path, label, *, error_type=ValueError):
    """Read a nonempty unique NumPy text-name vector once per input identity."""
    path = require_nonempty_file(path, label, error_type=error_type)

    def inspect():
        try:
            values = np.atleast_1d(np.loadtxt(path, dtype=str)).reshape(-1)
        except (OSError, ValueError) as exc:
            raise error_type("Cannot read %s %s: %s" % (label, path, exc)) from exc
        if values.size == 0 or any(not str(value).strip() for value in values):
            raise error_type("%s contains no usable names: %s" % (label, path))
        if len(values) != len(set(values.tolist())):
            raise error_type("%s contains duplicate names: %s" % (label, path))
        values.flags.writeable = False
        record_file_validation(
            path, "Unique name vector", checks=("nonempty unique names",),
            metrics={"names": len(values)},
        )
        return values

    return validate_once((path,), {"validator": "numpy_unique_name_vector"}, inspect,
                         error_type=error_type)


def matrix_values_are_finite(matrix):
    """Scan a numeric array in buffered blocks, without a full-size Boolean mask."""
    chunks = np.nditer(
        matrix, flags=["external_loop", "buffered"],
        op_flags=[["readonly"]], order="K",
    )
    return all(np.isfinite(chunk).all() for chunk in chunks)


def validate_feature_matrix_bundle(rows_path, chunks, *, error_type=ValueError):
    """Validate NPY chunks against their explicit row/column name companions."""
    rows_path = Path(rows_path).expanduser().resolve()
    chunks = tuple(tuple(Path(path).expanduser().resolve() for path in pair) for pair in chunks)
    if not chunks or any(len(pair) != 2 for pair in chunks):
        raise error_type("Feature matrix bundle requires column-name/matrix companion pairs")

    def inspect():
        rows = read_unique_names(rows_path, "Feature rows", error_type=error_type)
        all_columns = []
        chunk_metrics = []
        for chunk, (columns_path, matrix_path) in enumerate(chunks):
            columns = read_unique_names(columns_path, "Feature columns chunk %d" % chunk, error_type=error_type)
            try:
                matrix = np.load(matrix_path, mmap_mode="r", allow_pickle=False)
            except (OSError, ValueError) as exc:
                raise error_type("Cannot read feature matrix %s: %s" % (matrix_path, exc)) from exc
            if matrix.ndim != 2 or matrix.shape != (len(rows), len(columns)):
                raise error_type(
                    "Feature matrix chunk %d has shape %s; expected (%d, %d)."
                    % (chunk, matrix.shape, len(rows), len(columns))
                )
            if not np.issubdtype(matrix.dtype, np.number):
                raise error_type("Feature matrix chunk %d must have a numeric data type; found %s." % (chunk, matrix.dtype))
            if not matrix_values_are_finite(matrix):
                raise error_type("Feature matrix chunk %d contains non-finite values: %s" % (chunk, matrix_path))
            all_columns.extend(columns.tolist())
            metrics = {
                "chunk": chunk, "rows": matrix.shape[0], "columns": matrix.shape[1],
                "columns_path": columns_path, "matrix_path": matrix_path, "dtype": str(matrix.dtype),
            }
            chunk_metrics.append(metrics)
            record_file_validation(
                matrix_path, "Feature matrix", checks=("companion dimensions", "all values finite and numeric"),
                metrics={name: metrics[name] for name in ("rows", "columns", "dtype")},
            )
        if len(all_columns) != len(set(all_columns)):
            raise error_type("Feature names must be unique across all matrix chunks.")
        return {
            "rows": rows, "rows_path": rows_path, "row_count": len(rows),
            "feature_names": tuple(all_columns), "feature_count": len(all_columns),
            "chunks": chunk_metrics,
        }

    result = validate_once(
        (rows_path, *(path for pair in chunks for path in pair)),
        {"validator": "named_npy_matrix_bundle", "chunks": chunks}, inspect,
        error_type=error_type,
    )
    # Consumers attach run-specific compatibility/control evidence; do not let
    # that mutate the cached matrix contract or another consumer's metrics.
    return {**result, "chunks": [dict(chunk) for chunk in result["chunks"]]}


def validate_binary_kernel(
    path, *, dimension, dtype, bytes_per_value, chunk_rows,
    relative_tolerance, absolute_tolerance, label, error_type=ValueError,
):
    """Check a C-order square binary kernel; symmetry does not establish PSD."""
    path = require_nonempty_file(path, label, error_type=error_type)
    dtype = np.dtype(dtype)
    if dimension <= 0 or chunk_rows <= 0 or dtype.itemsize != bytes_per_value:
        raise error_type("Binary kernel dimensions, chunk rows or dtype/byte-width contract are invalid")
    if not all(np.isfinite(value) and value >= 0 for value in (relative_tolerance, absolute_tolerance)):
        raise error_type("Binary kernel symmetry tolerances must be finite and nonnegative")

    def inspect():
        expected_bytes = dimension ** 2 * bytes_per_value
        if path.stat().st_size != expected_bytes:
            raise error_type("%s has %d bytes; expected %d for %d %s genes" % (
                label, path.stat().st_size, expected_bytes, dimension, dtype,
            ))
        kernel = np.memmap(path, dtype=dtype, mode="r", shape=(dimension, dimension), order="C")
        for start in range(0, dimension, chunk_rows):
            stop = min(dimension, start + chunk_rows)
            block = np.asarray(kernel[start:stop, :])
            if not np.isfinite(block).all():
                raise error_type("%s contains non-finite values" % label)
            if not np.allclose(block, np.asarray(kernel[:, start:stop]).T,
                               rtol=relative_tolerance, atol=absolute_tolerance):
                raise error_type("%s matrix must be symmetric" % label)
        metrics = {"dimension": dimension, "bytes": expected_bytes, "dtype": str(dtype)}
        record_file_validation(
            path, "Binary kernel", checks=("binary dimensions", "all values finite", "configured symmetry tolerance"),
            metrics=metrics, message="Positive semidefiniteness is not established by symmetry.",
        )
        return metrics

    return dict(validate_once((path,), {
        "validator": "binary_kernel", "dimension": dimension, "dtype": dtype.str,
        "bytes_per_value": bytes_per_value, "chunk_rows": chunk_rows,
        "relative_tolerance": relative_tolerance, "absolute_tolerance": absolute_tolerance,
    }, inspect, error_type=error_type))
