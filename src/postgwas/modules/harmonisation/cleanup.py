"""Merging and clean-up of the per-chromosome intermediate files.

Every ``os.system`` call is gone.  ``os.system`` never raises, so the
``try/except`` around each one was dead code and each return status was thrown
away: a failed move was indistinguishable from a successful one, and the
delete that followed ran regardless.  Everything here now uses
``pathlib``/``shutil`` or ``subprocess.run`` on an argument list, so no value
out of the user's config CSV is ever handed to a shell.

Three ordering rules are load-bearing:

1. Each merged side output is written to a temporary file and renamed into
   place.  All required current-run sources must exist before any merge begins;
   otherwise finalisation fails instead of accepting a file left by an earlier
   run.  Shell ``>`` truncated the destination *before* ``cat`` ran, so a
   re-run, a second config row sharing an output directory, or a clean run with
   no errors destroyed the previously merged file.
2. Files are deleted only after the copies at their destination have been
   verified.  The gwas2vcf input files matched both the move glob and the
   delete glob, so a failed move was followed by an unrecoverable delete.
3. Per-chromosome sources are deleted only when their content actually reached
   the merged file.
"""

import csv
import gzip
import json
import os
import shutil
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from postgwas.core.paths import (
    configured_output_matches,
    configured_output_path,
    resolve_executable,
)
from postgwas.core.processes import run_checked_command

from .gwas2vcf_export import (
    GWAS2VCF_SUMMARY_CHROMOSOME_COLUMN,
    GWAS2VCF_SUMMARY_COLUMN_COUNT_COLUMN,
    GWAS2VCF_SUMMARY_COLUMNS,
    GWAS2VCF_SUMMARY_KEY_COLUMN,
    GWAS2VCF_SUMMARY_MISSING_COLUMN,
    GWAS2VCF_SUMMARY_ROW_COUNT_COLUMN,
    GWAS2VCF_SUMMARY_STATUS_COLUMN,
    GWAS2VCF_SUMMARY_SUCCESS_STATUS,
)
from .shared.runtime import emit_message


__all__ = [
    "finalise_harmonisation_outputs",
    "remove_merged_gwas2vcf_intermediate",
    "remove_partial_chromosome_outputs",
]


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
_emit = partial(emit_message, screen_when_unlogged=True)


def _safe_name(value, what):
    """A config value used to build a filename must not carry a path in it."""
    text = str(value).strip()
    if not text:
        raise ValueError("%s is empty; it is used to build output filenames." % what)
    if os.sep in text or (os.altsep and os.altsep in text) or "\0" in text:
        raise ValueError(
            "%s = %r contains a path separator, which would write outside the "
            "output directory." % (what, text)
        )
    return text


def _summary_integer(value, field, source, row_number, minimum):
    """Parse one generated adapter-summary integer without numeric coercion."""
    text = "" if value is None else str(value).strip()
    try:
        number = int(text)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            "GWAS-to-VCF summary %s row %d has invalid %s=%r; expected an "
            "integer." % (source, row_number, field, value)
        ) from exc
    if number < minimum:
        raise RuntimeError(
            "GWAS-to-VCF summary %s row %d has invalid %s=%d; expected a "
            "value of at least %d."
            % (source, row_number, field, number, minimum)
        )
    return number


def _merge_adapter_summaries(
    sources,
    destination,
    delimiter,
    expected_rows_by_chromosome,
    logger=None,
):
    """Atomically merge validated chromosome audits with exactly one header."""
    if not isinstance(delimiter, str) or len(delimiter) != 1:
        raise RuntimeError(
            "GWAS-to-VCF summary delimiter must be exactly one character."
        )
    expected_columns = list(GWAS2VCF_SUMMARY_COLUMNS)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part%d" % os.getpid())
    rows_by_chromosome = {}
    try:
        with temporary.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(
                output,
                fieldnames=expected_columns,
                delimiter=delimiter,
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for expected_chromosome, source in sources:
                source = Path(source)
                with source.open("r", encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(
                        handle, delimiter=delimiter, strict=True,
                    )
                    if reader.fieldnames != expected_columns:
                        raise RuntimeError(
                            "GWAS-to-VCF summary %s has header %s; expected %s. "
                            "The previous merged audit was left unchanged."
                            % (
                                source,
                                reader.fieldnames,
                                expected_columns,
                            )
                        )
                    source_rows = 0
                    observed_row_counts = set()
                    observed_column_counts = set()
                    observed_keys = set()
                    for row_number, row in enumerate(reader, start=2):
                        if None in row or any(
                            row.get(column) is None for column in expected_columns
                        ):
                            raise RuntimeError(
                                "GWAS-to-VCF summary %s row %d does not match its "
                                "validated header." % (source, row_number)
                            )
                        chromosome = str(
                            row[GWAS2VCF_SUMMARY_CHROMOSOME_COLUMN]
                        ).strip()
                        if chromosome != expected_chromosome:
                            raise RuntimeError(
                                "GWAS-to-VCF summary %s row %d reports chromosome "
                                "%r; expected chromosome %s."
                                % (
                                    source,
                                    row_number,
                                    chromosome,
                                    expected_chromosome,
                                )
                            )
                        status = str(
                            row[GWAS2VCF_SUMMARY_STATUS_COLUMN]
                        ).strip()
                        if status != GWAS2VCF_SUMMARY_SUCCESS_STATUS:
                            raise RuntimeError(
                                "GWAS-to-VCF summary %s row %d has status %r; "
                                "only successfully exported chromosome audits may "
                                "be merged."
                                % (source, row_number, status)
                            )
                        key = str(row[GWAS2VCF_SUMMARY_KEY_COLUMN]).strip()
                        if not key or key in observed_keys:
                            raise RuntimeError(
                                "GWAS-to-VCF summary %s row %d has an empty or "
                                "duplicate mapping key %r."
                                % (source, row_number, key)
                            )
                        observed_keys.add(key)
                        num_rows = _summary_integer(
                            row[GWAS2VCF_SUMMARY_ROW_COUNT_COLUMN],
                            GWAS2VCF_SUMMARY_ROW_COUNT_COLUMN,
                            source,
                            row_number,
                            minimum=0,
                        )
                        num_cols = _summary_integer(
                            row[GWAS2VCF_SUMMARY_COLUMN_COUNT_COLUMN],
                            GWAS2VCF_SUMMARY_COLUMN_COUNT_COLUMN,
                            source,
                            row_number,
                            minimum=1,
                        )
                        n_missing = _summary_integer(
                            row[GWAS2VCF_SUMMARY_MISSING_COLUMN],
                            GWAS2VCF_SUMMARY_MISSING_COLUMN,
                            source,
                            row_number,
                            minimum=0,
                        )
                        if n_missing > num_rows:
                            raise RuntimeError(
                                "GWAS-to-VCF summary %s row %d reports n_missing=%d "
                                "for only %d exported row(s)."
                                % (source, row_number, n_missing, num_rows)
                            )
                        observed_row_counts.add(num_rows)
                        observed_column_counts.add(num_cols)
                        writer.writerow(
                            dict((column, row[column]) for column in expected_columns)
                        )
                        source_rows += 1
                    if source_rows == 0:
                        raise RuntimeError(
                            "GWAS-to-VCF summary %s has a header but no audit rows."
                            % source
                        )
                    if len(observed_row_counts) != 1:
                        raise RuntimeError(
                            "GWAS-to-VCF summary %s contains inconsistent num_rows "
                            "values: %s."
                            % (source, sorted(observed_row_counts))
                        )
                    if len(observed_column_counts) != 1:
                        raise RuntimeError(
                            "GWAS-to-VCF summary %s contains inconsistent num_cols "
                            "values: %s."
                            % (source, sorted(observed_column_counts))
                        )
                    observed_rows = next(iter(observed_row_counts))
                    expected_rows = expected_rows_by_chromosome.get(
                        expected_chromosome
                    )
                    if expected_rows is not None and observed_rows != expected_rows:
                        raise RuntimeError(
                            "GWAS-to-VCF summary %s reports %d exported row(s) for "
                            "chromosome %s, but chromosome reconciliation reports %d. "
                            "The audit and scientific row accounting disagree."
                            % (
                                source,
                                observed_rows,
                                expected_chromosome,
                                expected_rows,
                            )
                        )
                    rows_by_chromosome[expected_chromosome] = observed_rows
        os.replace(str(temporary), str(destination))
    except RuntimeError:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise RuntimeError(
            "Cannot write the validated GWAS-to-VCF summary audit to %s "
            "(%s: %s). The previous merged audit and chromosome source audits "
            "were left unchanged."
            % (destination, type(exc).__name__, exc)
        ) from exc
    total_rows = sum(rows_by_chromosome.values())
    _emit(
        logger,
        "Validated and merged %d GWAS-to-VCF summary file(s) into %s: %d "
        "chromosome(s), %d adapter-input row(s), one header."
        % (
            len(sources),
            destination,
            len(rows_by_chromosome),
            total_rows,
        ),
    )
    return {
        "path": destination,
        "rows_by_chromosome": rows_by_chromosome,
        "total_rows": total_rows,
    }


def _merge_identical_json_mappings(sources, destination, logger=None):
    # type: (Sequence[Path], Path, Any) -> Optional[Path]
    """Validate chromosome mappings and write one atomic JSON document."""
    sources = [Path(source) for source in sources if Path(source).is_file()]
    if not sources:
        _emit(
            logger,
            "Nothing matched for %s, so the existing file was left untouched."
            % destination.name,
        )
        return None

    mappings = []
    for source in sources:
        try:
            with source.open("r", encoding="utf-8") as handle:
                mapping = json.load(handle)
        except (OSError, ValueError) as exc:
            raise RuntimeError(
                "Cannot read GWAS-to-VCF column mapping %s as JSON: %s."
                % (source, exc)
            ) from exc
        if not isinstance(mapping, dict):
            raise RuntimeError(
                "GWAS-to-VCF column mapping %s must contain one JSON object."
                % source
            )
        mappings.append(mapping)

    first = mappings[0]
    differing = [
        str(source)
        for source, mapping in zip(sources[1:], mappings[1:])
        if mapping != first
    ]
    if differing:
        raise RuntimeError(
            "Per-chromosome GWAS-to-VCF column mappings are inconsistent. "
            "The first mapping is %s; differing mapping file(s): %s. Column "
            "positions must be identical before one dataset-level mapping can "
            "be reported."
            % (sources[0], ", ".join(differing))
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part%d" % os.getpid())
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(first, handle, indent=2)
            handle.write("\n")
        os.replace(str(temporary), str(destination))
    except OSError as exc:
        try:
            temporary.unlink()
        except OSError:
            pass
        _emit(
            logger,
            "Could not write %s (%s); the previous file is unchanged."
            % (destination, exc),
            warn=True,
        )
        return None
    _emit(
        logger,
        "Validated %d chromosome mapping file(s) and wrote %s."
        % (len(sources), destination),
    )
    return destination


def _unlink(path, logger=None):
    """Remove one file or directory. Returns True when something went away."""
    path = Path(path)
    try:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(str(path))
        else:
            path.unlink()
    except FileNotFoundError:
        return False
    except OSError as exc:
        _emit(logger, "Could not remove %s: %s" % (path, exc), warn=True)
        return False
    return True


def remove_merged_gwas2vcf_intermediate(
    vcf_path, *, output_directory, logger=None,
):
    """Delete one validated raw adapter VCF and its exact tabix/CSI indexes.

    The caller supplies the path resolved from ``output_layout.merged_raw_vcf``.
    Wildcards are deliberately forbidden so cleanup cannot affect another study
    or either final build-specific VCF. A cleanup error is fatal: a successful
    run must not claim that the configured intermediate was removed when it was
    left behind.
    """
    output_root = Path(output_directory).expanduser().resolve()
    raw_vcf = Path(vcf_path).expanduser()
    try:
        raw_vcf.resolve(strict=False).relative_to(output_root)
    except ValueError as exc:
        raise RuntimeError(
            "Refusing to remove GWAS-to-VCF intermediate outside the dataset "
            "output directory %s: %s." % (output_root, raw_vcf)
        ) from exc
    if raw_vcf.is_symlink():
        raise RuntimeError(
            "Refusing to remove the GWAS-to-VCF intermediate because its "
            "configured path is a symbolic link: %s." % raw_vcf
        )
    if not raw_vcf.is_file() or raw_vcf.stat().st_size == 0:
        raise RuntimeError(
            "Cannot remove the validated GWAS-to-VCF intermediate because it "
            "is missing or empty: %s." % raw_vcf
        )

    candidates = [
        raw_vcf,
        Path(str(raw_vcf) + ".tbi"),
        Path(str(raw_vcf) + ".csi"),
    ]
    for path in candidates[1:]:
        if not path.exists():
            continue
        if path.is_symlink() or not path.is_file():
            raise RuntimeError(
                "Refusing to remove an invalid index path for the GWAS-to-VCF "
                "intermediate: %s." % path
            )
    removed = []
    for path in candidates:
        if not path.exists():
            continue
        try:
            path.unlink()
        except OSError as exc:
            raise RuntimeError(
                "Could not remove the GWAS-to-VCF intermediate artifact %s: %s."
                % (path, exc)
            ) from exc
        removed.append(path)

    _emit(
        logger,
        "Removed raw GWAS-to-VCF intermediate after successful final "
        "validation: %s%s."
        % (
            raw_vcf,
            " (and %d index file(s))" % (len(removed) - 1)
            if len(removed) > 1 else "",
        ),
    )
    return [str(path) for path in removed]


def _compress(paths, threads, executable=None, logger=None):
    """Compress each file with pigz, falling back to Python's gzip."""
    paths = [Path(p) for p in paths if Path(p).is_file()]
    if not paths:
        return []
    done = []
    pigz = None
    if executable:
        try:
            pigz = resolve_executable(executable, "pigz executable")
        except RuntimeError as exc:
            _emit(
                logger,
                "%s; compressing with Python gzip instead." % exc,
                warn=True,
            )
    if pigz:
        try:
            run_checked_command(
                [pigz, "-p", str(threads), *[str(path) for path in paths]],
                "Compressing archived GWAS-to-VCF inputs",
                logger=logger,
                error_type=RuntimeError,
            )
            return [p.with_name(p.name + ".gz") for p in paths]
        except RuntimeError as exc:
            _emit(
                logger,
                "pigz failed (%s); compressing with Python instead. The files "
                "themselves are safe either way." % exc,
                warn=True,
            )
    for path in paths:
        target = path.with_name(path.name + ".gz")
        try:
            with path.open("rb") as raw, gzip.open(str(target), "wb") as packed:
                shutil.copyfileobj(raw, packed)
        except OSError as exc:
            _emit(logger, "Could not compress %s: %s" % (path, exc), warn=True)
            continue
        path.unlink()
        done.append(target)
    return done


# ======================================================================
# 1. Merge and clean up after a whole dataset
# ======================================================================
def finalise_harmonisation_outputs(
    output_dir: str,
    gwas_outputname,
    output_layout: Dict[str, str],
    summary_delimiter: str,
    expected_chromosomes: Sequence[str],
    expected_rows_by_chromosome: Dict[str, Optional[int]],
    threads: int,
    compression_executable: Optional[str],
    policies=None,
    logger=None,
):
    """
    Merge the per-chromosome gwas2vcf side files, park the gwas2vcf inputs, and
    remove the per-chromosome intermediates that have been merged.

    Parameters
    ----------
    output_dir : str
        Directory holding the per-chromosome files.
    gwas_outputname : str
        Prefix used in filenames.
    summary_delimiter : str
        Resolved delimiter used by the chromosome adapter-summary files.
    expected_chromosomes : sequence of str
        Completed chromosomes whose side outputs belong in the final audit.
    expected_rows_by_chromosome : mapping
        Independently reconciled exported rows, or ``None`` only when rejection
        collection was explicitly disabled and reconciliation was unavailable.
    policies : Policies, optional
    logger : PipelineLogger, optional
        When absent the messages are printed, exactly as before.

    Returns
    -------
    dict
        The merged mapping and summary paths, plus what was moved and removed.
    """
    outdir = Path(output_dir)
    name = _safe_name(gwas_outputname, "gwas_outputname")
    chromosomes = [
        _safe_name(chromosome, "expected chromosome")
        for chromosome in expected_chromosomes
    ]
    if not chromosomes:
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: no "
            "completed chromosomes were provided." % name
        )
    if len(chromosomes) != len(set(chromosomes)):
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: the "
            "completed chromosome list contains duplicates: %s."
            % (name, ", ".join(chromosomes))
        )
    supplied_expected = dict(expected_rows_by_chromosome or {})
    if set(supplied_expected) != set(chromosomes):
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: expected "
            "row accounting must name exactly the completed chromosomes %s; "
            "received %s."
            % (
                name,
                ", ".join(chromosomes),
                ", ".join(sorted(str(key) for key in supplied_expected)),
            )
        )
    expected_rows = {}
    for chromosome in chromosomes:
        value = supplied_expected[chromosome]
        if value is None:
            expected_rows[chromosome] = None
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise RuntimeError(
                "Cannot finalise GWAS-to-VCF side outputs for dataset %s: "
                "independent row accounting for chromosome %s must be a "
                "non-negative integer or unavailable, found %r."
                % (name, chromosome, value)
            )
        expected_rows[chromosome] = value

    qc_dir = configured_output_path(outdir, output_layout["qc_directory"])
    qc_dir.mkdir(parents=True, exist_ok=True)

    values = {"dataset_id": name, "chromosome": "*"}
    mapping_path = configured_output_path(
        outdir, output_layout["adapter_merged_mapping"], **values,
    )
    summary_path = configured_output_path(
        outdir, output_layout["gwas2vcf_summary"], **values,
    )

    result = {
        "dict": None,
        "summary": None,
        "moved": [],
        "removed": [],
        "not_moved": [],
    }

    # ------------------------------------------------------------------
    # 1. Validate and merge only completed chromosomes' side outputs.
    #    The summary merge writes one header, validates its scientific row
    #    accounting, and replaces the delivered audit atomically.
    # ------------------------------------------------------------------
    mapping_sources = [
        configured_output_path(
            outdir,
            output_layout["adapter_mapping"],
            dataset_id=name,
            chromosome=chromosome,
        )
        for chromosome in chromosomes
    ]
    summary_sources = [
        configured_output_path(
            outdir,
            output_layout["adapter_summary"],
            dataset_id=name,
            chromosome=chromosome,
        )
        for chromosome in chromosomes
    ]
    missing = [
        str(path)
        for path in mapping_sources + summary_sources
        if not path.is_file() or path.stat().st_size == 0
    ]
    if missing:
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: required "
            "current per-chromosome side output(s) are missing or empty: %s. "
            "Existing merged files were not accepted because they may belong "
            "to an earlier run. Re-run the missing chromosome work or use a "
            "clean output directory."
            % (name, ", ".join(missing))
        )

    mapping_written = _merge_identical_json_mappings(
        mapping_sources, mapping_path, logger=logger,
    )
    if mapping_written is None:
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: the "
            "current mapping files could not be merged into %s. Any existing "
            "destination was not accepted as current-run output."
            % (name, mapping_path)
        )
    summary_result = _merge_adapter_summaries(
        list(zip(chromosomes, summary_sources)),
        summary_path,
        summary_delimiter,
        expected_rows,
        logger=logger,
    )
    result["dict"] = str(mapping_written)
    result["summary"] = str(summary_result["path"])
    result["adapter_input_rows_by_chromosome"] = dict(
        summary_result["rows_by_chromosome"]
    )
    result["total_adapter_input_rows"] = int(summary_result["total_rows"])
    merged_sources = {
        "dict": mapping_sources,
        "summary": summary_sources,
    }

    # ------------------------------------------------------------------
    # 2. Park the gwas2vcf input files, then compress them in place.
    # ------------------------------------------------------------------
    target_dir = configured_output_path(
        outdir, output_layout["adapter_input_archive"], **values,
    )
    vcf_inputs = sorted(configured_output_matches(
        outdir, output_layout["adapter_input"], **values,
    ))
    moved = []
    not_moved = []
    try:
        target_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        # Without a destination there is nowhere to verify a move, so the
        # inputs stay exactly where they are and nothing below deletes them.
        _emit(
            logger,
            "Could not create %s (%s), so the gwas2vcf input files were left in "
            "place and nothing was deleted." % (target_dir, exc),
            warn=True,
        )
        not_moved = list(vcf_inputs) or [outdir]
        vcf_inputs = []
    for source in vcf_inputs:
        destination = target_dir / source.name
        try:
            shutil.move(str(source), str(destination))
        except (OSError, shutil.Error) as exc:
            not_moved.append(source)
            _emit(logger, "Could not move %s: %s" % (source, exc), warn=True)
            continue
        # Verified at the destination BEFORE anything is deleted.
        if destination.is_file() and destination.stat().st_size > 0:
            moved.append(destination)
        else:
            not_moved.append(source)
            _emit(
                logger,
                "%s reported success but nothing arrived at %s."
                % (source, destination),
                warn=True,
            )

    if moved:
        _compress(
            moved, threads, executable=compression_executable, logger=logger,
        )
    result["moved"] = [str(p) for p in moved]
    result["not_moved"] = [str(p) for p in not_moved]

    # ------------------------------------------------------------------
    # 3. Remove the per-chromosome tables.
    #    '{name}*_chr*.tsv' also matches '{name}_chr7_vcf_input.tsv', so this
    #    used to delete the gwas2vcf inputs whenever the move above had failed
    #    - unrecoverably, and without an error.
    # ------------------------------------------------------------------
    removed = []
    move_complete = not not_moved
    if not move_complete:
        _emit(
            logger,
            "%d gwas2vcf input file(s) are not at their destination, so nothing "
            "was deleted. Fix the move and re-run the clean-up." % len(not_moved),
            warn=True,
        )
    else:
        for pattern_name in (
            "chromosome_table", "chromosome_source_snapshot",
            "external_eaf_partition", "external_info_partition",
        ):
            for candidate in sorted(configured_output_matches(
                outdir, output_layout[pattern_name], **values,
            )):
                if not candidate.is_file():
                    continue
                if _unlink(candidate, logger=logger):
                    removed.append(candidate)

        # ------------------------------------------------------------------
        # 4. Remove the per-chromosome side files, but only the ones whose
        #    content actually reached the merged file.
        # ------------------------------------------------------------------
        for key, pattern_name in (
            ("dict", "adapter_mapping"),
            ("summary", "adapter_summary"),
        ):
            merged = set(str(p) for p in merged_sources.get(key, []))
            if not merged:
                continue
            for candidate in sorted(configured_output_matches(
                outdir, output_layout[pattern_name], **values,
            )):
                if str(candidate) not in merged:
                    continue
                if _unlink(candidate, logger=logger):
                    removed.append(candidate)

    result["removed"] = [str(p) for p in removed]
    _emit(
        logger,
        "Clean-up finished: %d file(s) parked, %d removed."
        % (len(moved), len(removed)),
    )
    return result


# ======================================================================
# 2. Clean one chromosome before a retry  (plan Part 7)
# ======================================================================
def remove_partial_chromosome_outputs(
    outdir, sample_id, chromosome, output_layout, logger=None,
):
    # type: (Any, Any, Any, Any) -> int
    """Remove one chromosome's partial outputs before it is retried.

    Without this, round one's truncated ``chr7_GRCh37_CSQ.vcf.gz`` survives a
    round-two failure and reaches the merge, because ``validate_and_fix_vcf``
    only rejects files under 100 bytes.

    The validated chromosome input table and immutable source snapshot are
    deliberately retained because every retry reads them again. The log is
    also retained: it is opened in append mode, so every attempt for this
    chromosome stays in the one file.

    Returns the number of files and directories removed.
    """
    base = Path(outdir)
    sample = _safe_name(sample_id, "sample_id")
    chrom = _safe_name(chromosome, "chromosome")
    if not base.is_dir():
        _emit(logger, "Nothing to clean: %s does not exist." % base)
        return 0

    victims = []
    seen = set()
    values = {
        "dataset_id": sample,
        "chromosome": chrom,
        "build": "*",
        "target_build": "*",
    }
    for pattern_name in (
        "chromosome_reject", "adapter_input", "adapter_mapping",
        "adapter_summary", "adapter_output_vcf",
        "chromosome_raw_vcf", "chromosome_normalized_vcf",
        "chromosome_id_vcf",
        "chromosome_frequency_vcf", "chromosome_annotated_vcf",
        "chromosome_lifted_vcf", "chromosome_not_lifted_vcf",
        "missing_eaf", "out_of_range_eaf",
        "post_orientation_duplicates",
    ):
        for path in configured_output_matches(
            base, output_layout[pattern_name] + "*", **values,
        ):
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            victims.append(path)

    temp_dir = configured_output_path(
        base, output_layout["sort_temp_directory"], **values,
    )
    if temp_dir.is_dir():
        victims.append(temp_dir)

    removed = 0
    for path in victims:
        # Never remove a log: they are append-mode so that every retry attempt
        # is preserved in a single file.
        if path.suffix == ".log" or "logs" in path.parts:
            continue
        if _unlink(path, logger=logger):
            removed += 1

    _emit(
        logger,
        "Cleaned %d partial output(s) for chromosome %s before retrying; the log "
        "was kept." % (removed, chromosome),
    )
    return removed
