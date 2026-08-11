"""
rejects.py — every variant the pipeline removes, written out with all of its
original columns and the reason it went.

Plan v3, Part 4 and the Appendix reason registry.

Rules this module enforces rather than documents:

  * A variant appears exactly once.  The first rule that rejects it wins and it
    leaves the working frame immediately, so the per-reason counts sum to the
    total and the reconciliation balances.
  * The mask is applied as `mask.fill_null(True)`.  A null mask value rejects
    the row.  Without this a null allele makes the mask null, and because
    `filter` treats null as false the row is dropped by BOTH the keep filter and
    the reject filter and vanishes from every output - the live bug at
    io.py:526-534, verified as 3 rows in, 1 kept, 1 rejected, 1 gone.
  * The reason must be in REASONS.  An unregistered code raises, so a reason
    can never reach the output without a plain-English description.
  * Every rejection is logged through logger.qc(), which requires the counts
    before and after.
  * The file is flushed from a `finally`, so a chromosome that fails at step 9
    still has steps 1 to 8 on disk, and it is written even when empty (header
    only) so that an absent file means the step never ran.

Appended columns, in this order, reason LAST:

    reject_step, reject_detail, reject_reason

The stored columns come from an immutable snapshot captured immediately after
the study table is parsed and before any harmonisation. Rejected working rows
carry only a stable source-row ID into this module; the original values are
looked up from that snapshot when the reject file is materialised. The output
is written as text so per-chromosome files concatenate without dtype clashes.

This module imports nothing from postgwas.* - only the standard library and
polars - so it is importable and testable on its own.

Python 3.8 compatible.
"""

import csv
import gzip
import os
from types import MappingProxyType
from typing import Any, Dict, Iterable, List, Optional, Union

import polars as pl

__all__ = [
    "REASONS",
    "REASON_STEPS",
    "SOURCE_INPUT_ROW_COLUMN",
    "RejectCollector",
    "RejectOutputError",
    "RejectReasonError",
    "ReconciliationError",
    "concat_reject_files",
    "reconcile",
    "REJECT_COLUMNS",
]


# =============================================================================
# 1. Reason registry  (plan Appendix — reject reason registry)
# =============================================================================

_REASONS = {
    # -- step 01, reading -----------------------------------------------------
    "missing_required_columns":
        "A column the variant needs was empty or could not be parsed as its "
        "configured numeric type, so the variant could not be used.",

    # -- step 02, coordinates and alleles -------------------------------------
    "invalid_chromosome":
        "The chromosome is not one of the chromosomes being kept.",
    "invalid_position":
        "The position is missing, is not a positive number, or could not be read.",
    "non_standard_allele":
        "An allele is not written as plain A, C, G or T letters - insertions and "
        "deletions coded as I, D or - land here.",
    "null_allele":
        "The effect allele or the other allele is missing, so the variant cannot be "
        "oriented or matched to a reference allele pair.",
    "duplicate_variant":
        "Another row describes the same variant with consistent scientific values, "
        "and a deterministically better-ranked row was retained.",
    "conflicting_duplicate":
        "Rows share the configured duplicate key but disagree on one or more "
        "scientific values, so no row in the group can be selected safely.",

    # -- step 03, effect allele frequency -------------------------------------
    "eaf_out_of_range":
        "The frequency is below 0 or above 1 by more than rounding error.",
    "eaf_degenerate":
        "The frequency is exactly 0 or exactly 1, which makes the standard-error "
        "formula divide by zero.",
    "eaf_null":
        "The frequency is missing and could not be recovered from any reference.",
    "eaf_unmatched_external":
        "The variant is not in the external frequency reference.",
    "palindromic_ambiguous":
        "The variant is palindromic (A/T or C/G) and its frequency sits inside the "
        "band where both strand orientations are equally plausible.",
    "af_discordant":
        "The frequency disagrees with the reference panel under both orientations, "
        "so it cannot be trusted either way.",
    "reference_unmatched":
        "No REF/ALT allele orientation at this coordinate matches the configured "
        "chromosome reference.",
    "reference_ambiguous":
        "More than one reference allele orientation remains possible, so the "
        "effect allele cannot be selected safely.",

    # -- step 04, sample size --------------------------------------------------
    "neff_invalid":
        "The effective sample size is missing, zero, negative or infinite.",
    "sample_size_invalid":
        "The case or control count is missing or not a usable number.",

    # -- step 05, effect type --------------------------------------------------
    "effect_non_positive_or":
        "The odds ratio is zero or negative, so it has no logarithm.",

    # -- step 06, Z score ------------------------------------------------------
    "z_invalid":
        "The Z score is missing or is not a finite number.",

    # -- step 07, p-value ------------------------------------------------------
    "pval_out_of_range":
        "The p-value is outside the range from just above 0 to 1.",
    "pval_null":
        "The p-value is missing.",

    # -- step 08, standard error from p-value -------------------------------
    "zero_p_missing_se":
        "The reported p-value is zero and neither a supplied standard error nor a "
        "usable Z score is available, so an exact standard error cannot be "
        "reconstructed.",

    # -- step 10, effect statistics validation ---------------------------------
    "se_null":
        "The standard error is missing.",
    "se_non_finite":
        "The standard error is infinite or not a number. This is checked before the "
        "test for being positive, because infinity counts as greater than zero.",
    "se_non_positive":
        "The standard error is zero or negative.",
    "se_from_clipped_pval":
        "The standard error was worked out from a p-value that had to be clipped, "
        "so the exact significance was lost even if the derived value is positive "
        "and finite.",
    "beta_invalid":
        "The effect size is missing or is not a finite number.",
    "beta_zero":
        "The effect size is exactly zero.",
    "beta_se_z_discordant":
        "The supplied Z score does not agree with the effect size divided by its "
        "standard error within the configured rounding tolerances.",
    "z_pval_discordant":
        "The Z score and the p-value do not agree with each other.",

    # -- step 11, imputation quality -------------------------------------------
    "info_out_of_range":
        "The imputation quality is outside the allowed window.",
    "info_missing":
        "The imputation quality is missing.",

    # -- step 13, final completeness check --------------------------------------
    "final_missing_chr":
        "The chromosome is missing at the final check, immediately before export.",
    "final_missing_pos":
        "The position is missing at the final check, immediately before export.",
    "final_missing_snp":
        "The variant identifier is missing at the final check, immediately before export.",
    "final_missing_ea":
        "The effect allele is missing at the final check, immediately before export.",
    "final_missing_oa":
        "The other allele is missing at the final check, immediately before export.",
    "final_missing_eaf":
        "The effect allele frequency is missing at the final check, immediately "
        "before export.",
    "final_missing_beta":
        "The effect size is missing at the final check, immediately before export.",
    "final_missing_se":
        "The standard error is missing at the final check, immediately before export.",
    "final_missing_zscore":
        "The Z score is missing at the final check, immediately before export.",
    "final_missing_pval":
        "The p-value is missing at the final check, immediately before export.",
    "final_missing_info":
        "The imputation quality is missing at the final check when explicitly required.",
    "final_missing_n":
        "The configured sample size is missing at the final check, immediately before export.",
}

_REASON_STEPS = {
    "missing_required_columns": "01",
    "invalid_chromosome": "02",
    "invalid_position": "02",
    "non_standard_allele": "02",
    "null_allele": "02",
    "duplicate_variant": "02",
    "conflicting_duplicate": "02",
    "eaf_out_of_range": "03",
    "eaf_degenerate": "03",
    "eaf_null": "03",
    "eaf_unmatched_external": "03",
    "palindromic_ambiguous": "03",
    "af_discordant": "03",
    "reference_unmatched": "03",
    "reference_ambiguous": "03",
    "neff_invalid": "04",
    "sample_size_invalid": "04",
    "effect_non_positive_or": "05",
    "z_invalid": "06",
    "pval_out_of_range": "07",
    "pval_null": "07",
    "zero_p_missing_se": "08",
    "se_null": "10",
    "se_non_finite": "10",
    "se_non_positive": "10",
    "se_from_clipped_pval": "10",
    "beta_invalid": "10",
    "beta_zero": "10",
    "beta_se_z_discordant": "10",
    "z_pval_discordant": "10",
    "info_out_of_range": "11",
    "info_missing": "11",
    "final_missing_chr": "13",
    "final_missing_pos": "13",
    "final_missing_snp": "13",
    "final_missing_ea": "13",
    "final_missing_oa": "13",
    "final_missing_eaf": "13",
    "final_missing_beta": "13",
    "final_missing_se": "13",
    "final_missing_zscore": "13",
    "final_missing_pval": "13",
    "final_missing_info": "13",
    "final_missing_n": "13",
}

REASONS = MappingProxyType(dict((k, " ".join(v.split())) for k, v in _REASONS.items()))
REASON_STEPS = MappingProxyType(dict(_REASON_STEPS))

#: The three appended columns, in the order they are written.  Reason is last.
REJECT_COLUMNS = ("reject_step", "reject_detail", "reject_reason")

SOURCE_INPUT_ROW_COLUMN = "__postgwas_source_input_row"
_MASK_COLUMN = "__postgwas_reject_mask__"


class RejectReasonError(KeyError):
    """A reject reason that is not in the registry."""


class ReconciliationError(AssertionError):
    """Rows in did not equal rows out plus rows rejected.

    Derived from AssertionError because the plan requires reconciliation to be a
    hard assertion: if it does not balance, a variant leaked or was counted
    twice, and the run must stop.
    """


class RejectOutputError(RuntimeError):
    """Required rejected-variant provenance could not be written or read."""


# =============================================================================
# 2. The collector
# =============================================================================


class RejectCollector(object):
    """Collects every rejected variant for one chromosome (or one dataset).

    ``source_snapshot`` is either the immutable input-frame snapshot or the
    chromosome-specific Parquet snapshot written by the dataset stage. Only
    source-row IDs and rejection metadata are accumulated during scientific
    processing, keeping transformed values out of the provenance output.
    """

    def __init__(
        self,
        source_snapshot,
        logger,
        out_path,
        delimiter,
        compress=True,
        source_row_column=SOURCE_INPUT_ROW_COLUMN,
    ):
        self.source_row_column = str(source_row_column)
        self._source_frame = None
        self._source_path = None
        if isinstance(source_snapshot, pl.DataFrame):
            self._source_frame = source_snapshot.clone()
            source_schema = dict(source_snapshot.schema)
        elif isinstance(source_snapshot, (str, os.PathLike)):
            self._source_path = os.fspath(source_snapshot)
            if not os.path.isfile(self._source_path):
                raise RejectOutputError(
                    "Immutable source-row snapshot is missing: %s"
                    % self._source_path
                )
            try:
                source_schema = dict(
                    pl.scan_parquet(self._source_path).collect_schema()
                )
            except Exception as exc:
                raise RejectOutputError(
                    "Immutable source-row snapshot could not be read from %s "
                    "(%s: %s)."
                    % (self._source_path, type(exc).__name__, exc)
                ) from exc
        else:
            raise TypeError(
                "RejectCollector source_snapshot must be a Polars DataFrame or "
                "Parquet path, got %r" % type(source_snapshot).__name__
            )

        if self.source_row_column not in source_schema:
            raise RejectOutputError(
                "Immutable source-row snapshot has no stable row-ID column %r."
                % self.source_row_column
            )
        self.original_columns = [
            str(column)
            for column in source_schema
            if column != self.source_row_column
        ]
        collisions = sorted(set(self.original_columns) & set(REJECT_COLUMNS))
        if collisions:
            raise RejectOutputError(
                "Input column name(s) conflict with rejected-variant metadata: %s. "
                "Rename them before harmonisation."
                % ", ".join(collisions)
            )
        used_columns = set(source_schema) | set(REJECT_COLUMNS)
        self._event_order_column = "__postgwas_reject_event_order__"
        while self._event_order_column in used_columns:
            self._event_order_column += "_"
        used_columns.add(self._event_order_column)
        self._source_found_column = "__postgwas_reject_source_found__"
        while self._source_found_column in used_columns:
            self._source_found_column += "_"
        self.logger = logger
        self.compress = bool(compress)
        self.delimiter = str(delimiter)
        if len(self.delimiter) != 1:
            raise ValueError("Reject-file delimiter must be exactly one character")

        path = str(out_path)
        if self.compress and not path.endswith(".gz"):
            path += ".gz"
        self.out_path = path

        self.columns = self.original_columns + list(REJECT_COLUMNS)
        self._frames = []          # type: List[pl.DataFrame]
        self._counts = {}          # type: Dict[str, int]
        self._by_step = {}         # type: Dict[str, int]
        self._total = 0
        self.flushed_at = None     # type: Optional[str]
        self.flush_count = 0

    # -- the one way to remove a variant -----------------------------------
    def reject(self, df, mask, step, reason, detail=None):
        # type: (pl.DataFrame, Union[pl.Expr, pl.Series], str, str, Optional[str]) -> pl.DataFrame
        """Split `df` on `mask`, record the rejected rows, return the survivors.

        A row whose mask value is null is REJECTED, never dropped from both
        sides.  The mask is materialised once, so the keep half and the reject
        half are exact complements of each other by construction.
        """
        if reason not in REASONS:
            raise RejectReasonError(
                "'%s' is not a registered reject reason. Add it to REASONS in "
                "harmonisation/rejects.py with a plain-English description; a reason "
                "with no description must never reach the output. Registered reasons: %s"
                % (reason, ", ".join(sorted(REASONS)))
            )
        if not isinstance(df, pl.DataFrame):
            raise TypeError("reject() needs a polars DataFrame, got %r" % type(df).__name__)
        if not step:
            raise ValueError("reject() needs the step it was called from, for the record")

        before = df.height

        if isinstance(mask, pl.Series):
            if mask.len() != before:
                raise ValueError(
                    "the mask has %d values but the frame has %d rows"
                    % (mask.len(), before)
                )
            mask_expr = pl.lit(mask.fill_null(True)).alias(_MASK_COLUMN)
        elif isinstance(mask, pl.Expr):
            mask_expr = mask.fill_null(True).alias(_MASK_COLUMN)
        else:
            raise TypeError(
                "reject() needs a polars expression or a boolean Series as its mask, "
                "got %r" % type(mask).__name__
            )

        working = df.drop(_MASK_COLUMN) if _MASK_COLUMN in df.columns else df
        tagged = working.with_columns(mask_expr)
        if tagged.schema[_MASK_COLUMN] != pl.Boolean:
            tagged = tagged.with_columns(
                pl.col(_MASK_COLUMN).cast(pl.Boolean, strict=False).fill_null(True)
            )

        partitions = tagged.partition_by(
            _MASK_COLUMN, as_dict=True, maintain_order=True,
        )
        empty = tagged.head(0)
        rejected = partitions.get((True,), empty).drop(_MASK_COLUMN)
        survivors = partitions.get((False,), empty).drop(_MASK_COLUMN)
        after = survivors.height

        if rejected.height:
            self._record(rejected, step, reason, detail)

        if self.logger is not None:
            plain = REASONS[reason]
            if detail:
                plain = "%s (%s)" % (plain, detail)
            self.logger.qc(
                "%s" % reason.replace("_", " "),
                plain,
                before,
                after,
                reason=reason,
                step=step,
            )
        return survivors

    def _record(self, rejected, step, reason, detail):
        # type: (pl.DataFrame, str, str, Optional[str]) -> None
        if self.source_row_column not in rejected.columns:
            raise RejectOutputError(
                "A rejected working row has no stable source-row ID column %r, "
                "so its original study values cannot be recovered."
                % self.source_row_column
            )
        count = rejected.height
        events = (
            rejected.select(
                pl.col(self.source_row_column)
                .cast(pl.UInt64, strict=True)
                .alias(self.source_row_column)
            )
            .with_row_index(self._event_order_column, offset=self._total)
            .with_columns([
                pl.lit(str(step), dtype=pl.Utf8).alias("reject_step"),
                pl.lit(
                    None if detail is None else str(detail), dtype=pl.Utf8,
                ).alias("reject_detail"),
                pl.lit(str(reason), dtype=pl.Utf8).alias("reject_reason"),
            ])
        )
        self._frames.append(events)
        self._counts[reason] = self._counts.get(reason, 0) + count
        self._by_step[str(step)] = self._by_step.get(str(step), 0) + count
        self._total += count

    # -- accounting ---------------------------------------------------------
    def counts(self):
        # type: () -> Dict[str, int]
        """Rejected variants by reason, in the order the reasons first fired."""
        return dict(self._counts)

    def counts_by_step(self):
        # type: () -> Dict[str, int]
        return dict(self._by_step)

    def total(self):
        # type: () -> int
        """Every variant this collector has taken out of the frame."""
        return int(self._total)

    def summary(self):
        # type: () -> Dict[str, Any]
        """For the run manifest and the parent's per-chromosome footer."""
        return {
            "path": self.out_path,
            "total": self.total(),
            "by_reason": self.counts(),
            "by_step": self.counts_by_step(),
            "flushed": self.flush_count > 0,
        }

    def frame(self):
        # type: () -> pl.DataFrame
        """Materialise original values plus rejection metadata in event order."""
        if not self._frames:
            return pl.DataFrame(
                dict((name, pl.Series(name, [], dtype=pl.Utf8)) for name in self.columns)
            ).select(self.columns)
        events = (
            self._frames[0]
            if len(self._frames) == 1
            else pl.concat(self._frames, how="vertical")
        )
        if events.get_column(self.source_row_column).n_unique() != events.height:
            raise RejectOutputError(
                "A stable source-row ID was rejected more than once; rejected "
                "provenance cannot be written without duplicating a study row."
            )

        wanted_ids = events.select(self.source_row_column)
        source = (
            self._source_frame.lazy()
            if self._source_frame is not None
            else pl.scan_parquet(self._source_path)
        )
        try:
            source_rows = (
                source.select([
                    pl.col(self.source_row_column)
                    .cast(pl.UInt64, strict=True)
                    .alias(self.source_row_column),
                    *[
                        pl.col(name).cast(pl.Utf8, strict=False).alias(name)
                        for name in self.original_columns
                    ],
                ])
                .join(wanted_ids.lazy(), on=self.source_row_column, how="semi")
                .collect()
            )
        except Exception as exc:
            raise RejectOutputError(
                "Original values for rejected rows could not be loaded from the "
                "immutable source snapshot (%s: %s)."
                % (type(exc).__name__, exc)
            ) from exc

        if (
            source_rows.height != events.height
            or source_rows.get_column(self.source_row_column).n_unique()
            != source_rows.height
        ):
            found = source_rows.get_column(self.source_row_column).n_unique()
            raise RejectOutputError(
                "Immutable source-row snapshot matched %d unique row ID(s) for "
                "%d rejected variant(s). Every rejected row must match exactly once."
                % (found, events.height)
            )

        source_rows = source_rows.with_columns(
            pl.lit(True).alias(self._source_found_column)
        )
        joined = events.join(
            source_rows,
            on=self.source_row_column,
            how="left",
        )
        if joined.get_column(self._source_found_column).null_count():
            raise RejectOutputError(
                "At least one rejected source-row ID is absent from the immutable "
                "source snapshot."
            )
        return (
            joined.sort(self._event_order_column)
            .select(self.original_columns + list(REJECT_COLUMNS))
        )

    # -- output --------------------------------------------------------------
    def flush(self):
        # type: () -> str
        """Write the file.  Header-only when nothing was rejected.

        It is safe to call twice and rewrites the whole file each time.  A
        failed write raises because rejection provenance is a required output
        whenever collection is enabled; callers that are already handling an
        analysis failure may record this as a secondary error.
        """
        frame = self.frame()
        directory = os.path.dirname(os.path.abspath(self.out_path))
        try:
            if directory:
                os.makedirs(directory, exist_ok=True)
            _write_table(frame, self.out_path, self.compress, self.delimiter)
            self.flush_count += 1
            self.flushed_at = self.out_path
            if self.logger is not None:
                if frame.height:
                    self.logger.info(
                        "Wrote %s rejected variants, with all their original columns, to %s"
                        % ("{:,}".format(frame.height), self.out_path)
                    )
                else:
                    self.logger.info(
                        "No variants were rejected. Wrote the header-only reject file to %s, "
                        "so that an absent file means the step never ran." % self.out_path
                    )
        except Exception as exc:
            if self.logger is not None:
                self.logger.error(
                    "The required rejected-variants file could not be written to %s "
                    "(%s: %s)."
                    % (self.out_path, type(exc).__name__, exc)
                )
            raise RejectOutputError(
                "Required rejected-variants output could not be written to %s "
                "(%s: %s)." % (self.out_path, type(exc).__name__, exc)
            ) from exc
        return self.out_path

    def __repr__(self):
        return "<RejectCollector %s rejected, %d reasons, -> %s>" % (
            self.total(), len(self._counts), self.out_path
        )


# =============================================================================
# 3. Module-level helpers
# =============================================================================


def _write_table(frame, path, compress, delimiter):
    if compress:
        with gzip.open(path, "wb") as handle:
            frame.write_csv(handle, separator=delimiter)
    else:
        frame.write_csv(path, separator=delimiter)


def _read_table(path, delimiter):
    """Read a reject file with every column as text, so the schemas always match."""
    return pl.read_csv(
        path,
        separator=delimiter,
        infer_schema_length=0,
        has_header=True,
        truncate_ragged_lines=False,
    )


def _validate_reject_output(path, delimiter, expected_rows, compressed):
    """Validate one written reject table without loading it back into memory."""
    opener = gzip.open if compressed else open
    try:
        with opener(
            path, "rt", encoding="utf-8", newline="",
        ) as handle:
            reader = csv.reader(handle, delimiter=delimiter)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise RejectOutputError(
                    "Combined rejected-variants output is empty and has no header: %s"
                    % path
                ) from exc
            missing = [column for column in REJECT_COLUMNS if column not in header]
            if missing:
                raise RejectOutputError(
                    "Combined rejected-variants output %s is missing required "
                    "column(s): %s." % (path, ", ".join(missing))
                )
            observed_rows = sum(1 for _row in reader)
    except RejectOutputError:
        raise
    except Exception as exc:
        raise RejectOutputError(
            "Combined rejected-variants output could not be validated at %s "
            "(%s: %s)." % (path, type(exc).__name__, exc)
        ) from exc
    if observed_rows != int(expected_rows):
        raise RejectOutputError(
            "Combined rejected-variants output %s contains %d row(s), but %d "
            "were written from the source rejection files. Source files were "
            "not removed."
            % (path, observed_rows, int(expected_rows))
        )
    return observed_rows


def concat_reject_files(
    paths, out_path, delimiter, logger=None, compress=None, remove_sources=False,
):
    """Merge the per-chromosome reject files into the one dataset file.

    Missing or unreadable inputs fail finalisation: an absent file means that
    an expected stage left no auditable record. Header-only files contribute
    nothing but confirm that the stage ran and rejected no variants. When
    ``remove_sources`` is true, source shards are removed only after an atomic
    write has passed required-column and exact-row-count validation.
    """
    delimiter = str(delimiter)
    if len(delimiter) != 1:
        raise RejectOutputError(
            "Reject-file delimiter must be exactly one character."
        )
    if compress is None:
        compress = str(out_path).endswith(".gz")
    if compress and not str(out_path).endswith(".gz"):
        out_path = str(out_path) + ".gz"
    out_path = str(out_path)

    frames = []      # type: List[pl.DataFrame]
    used = []        # type: List[str]
    missing = []     # type: List[str]
    for path in paths:
        path = str(path)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            missing.append(path)
            continue
        try:
            frame = _read_table(path, delimiter)
        except Exception as exc:
            missing.append(path)
            if logger is not None:
                logger.error(
                    "The required reject file %s could not be read (%s: %s)."
                    % (path, type(exc).__name__, exc)
                )
            continue
        used.append(path)
        if frame.height:
            frames.append(frame)

    if missing:
        raise RejectOutputError(
            "%d required rejected-variant file%s %s missing or unreadable: %s"
            % (
                len(missing),
                "" if len(missing) == 1 else "s",
                "is" if len(missing) == 1 else "are",
                ", ".join(missing),
            )
        )

    if frames:
        merged = pl.concat(frames, how="diagonal")
    elif used:
        merged = _read_table(used[0], delimiter).head(0)
    else:
        merged = pl.DataFrame(
            dict((name, pl.Series(name, [], dtype=pl.Utf8)) for name in REJECT_COLUMNS)
        )

    ordered = [c for c in merged.columns if c not in REJECT_COLUMNS]
    ordered += [c for c in REJECT_COLUMNS if c in merged.columns]
    merged = merged.select(ordered)

    directory = os.path.dirname(os.path.abspath(out_path))
    temporary = out_path + ".part%d" % os.getpid()
    try:
        if directory:
            os.makedirs(directory, exist_ok=True)
        _write_table(merged, temporary, compress, delimiter)
        _validate_reject_output(
            temporary, delimiter, merged.height, compressed=compress,
        )
        os.replace(temporary, out_path)
    except RejectOutputError:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    except Exception as exc:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        if logger is not None:
            logger.error(
                "The combined rejected-variants file could not be written to %s "
                "(%s: %s)." % (out_path, type(exc).__name__, exc)
            )
        raise RejectOutputError(
            "Combined rejected-variants output could not be written to %s "
            "(%s: %s)." % (out_path, type(exc).__name__, exc)
        ) from exc

    removed_sources = []
    retained_sources = []
    if remove_sources:
        destination = os.path.abspath(out_path)
        for path in used:
            if os.path.abspath(path) == destination:
                continue
            try:
                os.unlink(path)
                removed_sources.append(path)
            except OSError as exc:
                retained_sources.append(path)
                if logger is not None:
                    warning = (
                        "The consolidated rejected-variants file is valid, but "
                        "source shard %s could not be removed (%s: %s)."
                        % (path, type(exc).__name__, exc)
                    )
                    method = getattr(logger, "warn", None) or logger.error
                    method(warning)

    result = {
        "path": out_path,
        "rows": merged.height,
        "files_used": used,
        "files_missing": missing,
        "source_files_removed": removed_sources,
        "source_files_retained": retained_sources,
    }
    if logger is not None:
        logger.info(
            "Combined %d reject file%s into %s: %s variants in total."
            % (len(used), "" if len(used) == 1 else "s", out_path,
               "{:,}".format(merged.height))
        )
        if remove_sources:
            logger.info(
                "Removed %d validated source rejection file%s; %d remain."
                % (
                    len(removed_sources),
                    "" if len(removed_sources) == 1 else "s",
                    len(retained_sources),
                )
            )
    return result


def reconcile(rows_in, rows_out, collector, logger=None, label=None):
    # type: (int, int, Any, Any, Optional[str]) -> Dict[str, Any]
    """Hard assertion: rows in, minus everything rejected, must equal rows out.

    If it does not balance a variant leaked or was counted twice, so this raises
    ReconciliationError (an AssertionError) naming the discrepancy.  Accepts a
    RejectCollector or a plain total.
    """
    if hasattr(collector, "total"):
        rejected = int(collector.total())
        by_reason = collector.counts() if hasattr(collector, "counts") else {}
    else:
        rejected = int(collector)
        by_reason = {}

    rows_in = int(rows_in)
    rows_out = int(rows_out)
    expected = rows_in - rejected
    difference = rows_out - expected
    where = label or "this chromosome"

    table = [
        "rows read ......................... %s" % "{:,}".format(rows_in),
        "  - rejected (all reasons) ....... %s" % "{:,}".format(rejected),
        "  = rows exported ................ %s" % "{:,}".format(rows_out),
    ]

    if difference != 0:
        if difference < 0:
            what = (
                "%s variants have vanished: they are neither in the exported data nor in "
                "the rejected file." % "{:,}".format(-difference)
            )
        else:
            what = (
                "%s variants have been counted twice: the exported data holds more rows "
                "than were read minus those rejected." % "{:,}".format(difference)
            )
        lines = [
            "The variant counts for %s do not balance." % where,
            "",
        ] + table + [
            "",
            what,
            "Expected %s rows exported, found %s."
            % ("{:,}".format(expected), "{:,}".format(rows_out)),
        ]
        if by_reason:
            lines.append("")
            lines.append("Rejected by reason:")
            for reason, count in sorted(by_reason.items(), key=lambda kv: (-kv[1], kv[0])):
                lines.append("  %-28s %s" % (reason, "{:,}".format(count)))
        message = "\n".join(lines)
        if logger is not None:
            logger.error(message)
        raise ReconciliationError(message)

    if logger is not None:
        logger.info("\n".join(table + ["balances"]))

    return {
        "rows_in": rows_in,
        "rows_out": rows_out,
        "rejected": rejected,
        "reject_counts": by_reason,
        "balanced": True,
    }
