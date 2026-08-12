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

from .shared.runtime import emit_message


__all__ = [
    "finalise_harmonisation_outputs",
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


def _concatenate(sources, destination, logger=None):
    # type: (Sequence[Path], Path, Any) -> Optional[Path]
    """Concatenate `sources` into `destination`, or do nothing if there are none.

    Written to a temporary file in the destination directory and renamed into
    place, so the destination is never left truncated or half written.
    """
    sources = [Path(s) for s in sources if Path(s).is_file()]
    if not sources:
        _emit(
            logger,
            "Nothing matched for %s, so the existing file was left untouched."
            % destination.name,
        )
        return None
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part%d" % os.getpid())
    try:
        with temporary.open("wb") as out:
            for source in sources:
                with source.open("rb") as handle:
                    shutil.copyfileobj(handle, out)
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
        "Merged %d file(s) into %s." % (len(sources), destination),
    )
    return destination


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
    # 1. Validate and merge the per-chromosome side outputs.
    #    Shell '>' truncated the destination before cat ran, so an empty
    #    glob left a 0-byte file where the previous merge used to be.
    # ------------------------------------------------------------------
    side_outputs = (
        (
            "dict", "mapping", "adapter_mapping", mapping_path,
            _merge_identical_json_mappings,
        ),
        (
            "summary", "summary", "adapter_summary", summary_path,
            _concatenate,
        ),
    )
    source_groups = {
        key: sorted(configured_output_matches(
            outdir, output_layout[pattern_name], **values,
        ))
        for key, _label, pattern_name, _destination, _merge in side_outputs
    }
    missing = [
        label
        for key, label, _pattern_name, _destination, _merge in side_outputs
        if not source_groups[key]
    ]
    if missing:
        raise RuntimeError(
            "Cannot finalise GWAS-to-VCF side outputs for dataset %s: no "
            "current per-chromosome %s file(s) were found. Existing merged "
            "files were not accepted because they may belong to an earlier "
            "run. Re-run the missing chromosome work or use a clean output "
            "directory."
            % (name, " or ".join(missing))
        )

    merged_sources = {}
    for key, label, _pattern_name, destination, merge in side_outputs:
        sources = source_groups[key]
        written = merge(sources, destination, logger=logger)
        if written is None:
            raise RuntimeError(
                "Cannot finalise GWAS-to-VCF side outputs for dataset %s: "
                "the current %s files could not be merged into %s. Any "
                "existing destination was not accepted as current-run output."
                % (name, label, destination)
            )
        result[key] = str(written)
        merged_sources[key] = sources

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
        "chromosome_raw_vcf", "chromosome_original_vcf",
        "chromosome_normalized_vcf", "chromosome_id_vcf",
        "chromosome_frequency_vcf", "chromosome_annotated_vcf",
        "chromosome_lifted_vcf", "chromosome_not_lifted_vcf",
        "missing_eaf", "out_of_range_eaf",
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
