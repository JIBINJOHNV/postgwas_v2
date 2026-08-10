"""PLINK reference inspection shared by MAGMA mapping strategies."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from postgwas.core.paths import require_nonempty_file
from postgwas.modules.magma.errors import MagmaError


def read_reference_bim_matches(
    ld_reference_prefix: str | Path,
    input_config,
    identifiers: set[str],
) -> tuple[pd.DataFrame, int]:
    """Stream one BIM and retain exact records requested by the caller."""
    bim = require_nonempty_file(
        "%s%s" % (ld_reference_prefix, input_config.bim_extension),
        "LD-reference BIM file",
        error_type=MagmaError,
    )
    roles = list(input_config.bim_columns)
    try:
        delimiter = re.compile(input_config.table_delimiter_pattern)
    except re.error as exc:
        raise MagmaError("Invalid configured BIM delimiter pattern: %s" % exc) from exc
    indexes = {role: roles.index(role) for role in roles}
    matched: dict[str, tuple[str, int, str]] = {}
    variant_count = 0
    try:
        with bim.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                fields = delimiter.split(raw.strip())
                if len(fields) != len(roles):
                    raise MagmaError(
                        "LD-reference BIM line %d has %d fields; expected %d."
                        % (line_number, len(fields), len(roles))
                    )
                identifier = fields[indexes["variant_id"]].strip()
                if not identifier:
                    raise MagmaError(
                        "LD-reference BIM contains an empty variant ID at line %d."
                        % line_number
                    )
                try:
                    position = int(fields[indexes["position"]])
                except ValueError as exc:
                    raise MagmaError(
                        "LD-reference BIM contains a non-integer position at line %d."
                        % line_number
                    ) from exc
                chromosome = fields[indexes["chromosome"]].strip()
                first = fields[indexes["allele1"]].strip().upper()
                second = fields[indexes["allele2"]].strip().upper()
                if (
                    position < 1
                    or chromosome in input_config.invalid_chromosome_labels
                    or first in input_config.invalid_allele_labels
                    or second in input_config.invalid_allele_labels
                ):
                    raise MagmaError(
                        "LD-reference BIM contains an invalid coordinate or allele "
                        "at line %d." % line_number
                    )
                variant_count += 1
                if identifier not in identifiers:
                    continue
                if identifier in matched:
                    raise MagmaError(
                        "LD-reference BIM repeats matched variant ID %s at line %d."
                        % (identifier, line_number)
                    )
                normalized = re.sub(
                    input_config.chromosome_prefix_pattern, "", chromosome,
                ).upper()
                normalized = input_config.chromosome_aliases.get(
                    normalized, normalized,
                )
                matched[identifier] = (
                    normalized,
                    position,
                    "|".join(sorted((first, second))),
                )
    except OSError as exc:
        raise MagmaError("Cannot read LD-reference BIM %s: %s" % (bim, exc)) from exc
    if variant_count == 0:
        raise MagmaError("LD-reference BIM contains no variants: %s" % bim)
    frame = pd.DataFrame.from_dict(
        matched,
        orient="index",
        columns=["CHR_NORM", "BP_NORM", "ALLELE_KEY"],
    )
    frame.index.name = "variant_id"
    return frame, variant_count


__all__ = ["read_reference_bim_matches"]
