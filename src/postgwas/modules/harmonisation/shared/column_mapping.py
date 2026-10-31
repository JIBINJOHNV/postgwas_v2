"""Scientific column roles shared by input checks and resolved transformations."""

from collections.abc import Iterable, Mapping

from postgwas.core.values import optional_text


# These keys define the sample-sheet/engine protocol, not guessed input headers.
COLUMN_CONFIG_KEYS = (
    "chr_col", "pos_col", "chr_pos_col", "snp_id_col",
    "ea_col", "oa_col", "eaf_col",
    "beta_or_col", "se_col", "imp_z_col", "pval_col",
    "ncontrol_col", "ncase_col", "imp_info_col",
)


def column_mapping_groups(
    mapping: Mapping, *, keys: Iterable[str] = COLUMN_CONFIG_KEYS,
) -> dict[str, list[str]]:
    """Group explicitly selected mappings by physical name, preserving order."""
    groups: dict[str, list[str]] = {}
    for key in keys:
        column = optional_text(mapping.get(key))
        if column is not None:
            groups.setdefault(column, []).append(key)
    return groups


def validate_scientific_column_mappings(
    mapping: Mapping, *, stage: str, error_type=ValueError,
) -> None:
    """Reject unrelated scientific roles sharing one resolved working column.

    ``beta_or_col`` and ``beta_col`` are the input/resolved names for the same
    effect role. Their alias is intentional; sharing an effect with coordinates,
    identifiers, uncertainty, frequency, or annotation values is never valid.
    This representation invariant cannot be waived by a header-warning policy.
    """
    conflicts = []
    for column, keys in column_mapping_groups(
        mapping, keys=(*COLUMN_CONFIG_KEYS, "beta_col"),
    ).items():
        roles = {"beta_or_col" if key == "beta_col" else key for key in keys}
        if len(roles) > 1:
            conflicts.append("%s -> %r" % (" and ".join(keys), column))
    if conflicts:
        raise error_type(
            "%s has incompatible scientific column mappings: %s. "
            "Assign each scientific role its own column; selected values "
            "must not be overwritten by another role."
            % (stage, "; ".join(conflicts))
        )
