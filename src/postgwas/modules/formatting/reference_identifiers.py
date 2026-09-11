"""Reference-driven formatter variant-ID selection shared by pipeline consumers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from postgwas.core.variant_identifiers import inspect_bim_identifier_type


@dataclass(frozen=True)
class BimIdentifierRequirement:
    """One downstream consumer's PLINK BIM identifier requirement."""

    consumer: str
    formatter_target: str
    bim_file: str | Path
    column_roles: Sequence[str]
    delimiter_pattern: str


def configure_reference_variant_identifiers(
    args,
    formatting_config,
    requirements: Sequence[BimIdentifierRequirement],
) -> dict[str, dict]:
    """Inspect required BIMs and configure each formatter target exactly once.

    Consumers that share a formatter artifact must agree on the observed BIM
    identifier convention. An explicit global formatter identifier type is
    also rejected when it conflicts with a required reference.
    """
    identifier_policy = formatting_config.variant_identifiers
    observations = dict(getattr(args, "variant_id_observations", None) or {})
    target_types = dict(identifier_policy.target_types)
    target_types.update(getattr(args, "variant_id_types", None) or {})
    detected_by_target = {
        target: (str(observation["variant_id_type"]), ", ".join(observation.get("consumers", (target,))))
        for target, observation in observations.items()
    }
    detected_by_target.update(getattr(args, "_variant_id_requirements", None) or {})

    for requirement in requirements:
        resolved_bim = Path(requirement.bim_file).expanduser().resolve()
        cached = observations.get(requirement.formatter_target, {})
        cached_consumers = set(cached.get("consumers", ()))
        # Reuse belongs to the shared file-identity/contract cache, never to a
        # formatter target or a path-only observation from an earlier call.
        summary = inspect_bim_identifier_type(
            resolved_bim,
            column_roles=list(requirement.column_roles),
            delimiter_pattern=requirement.delimiter_pattern,
            rsid_pattern=identifier_policy.rsid_pattern,
            unique_id_template=identifier_policy.unique_id_template,
            chromosome_prefix_pattern=(
                formatting_config.chromosome_labels.prefix_pattern
            ),
            chromosome_aliases=formatting_config.chromosome_labels.aliases,
        )
        detected_type = summary.identifier_type
        observation = {
            "bim_file": str(summary.bim_file),
            "variant_id_type": detected_type,
            "variants": summary.variants,
            "rsids": summary.rsids,
            "unique_ids": summary.unique_ids,
            "other_ids": summary.other_ids,
            "missing_ids": summary.missing_ids,
            "duplicate_ids": summary.duplicate_ids,
            "duplicate_rows": summary.duplicate_rows,
        }

        previous = detected_by_target.get(requirement.formatter_target)
        if previous is not None and previous[0] != detected_type:
            raise ValueError(
                "%s and %s share formatter target %s but their BIM files use "
                "different identifier conventions (%s versus %s). Use compatible "
                "references or run the modules separately."
                % (
                    previous[1], requirement.consumer,
                    requirement.formatter_target, previous[0], detected_type,
                )
            )
        detected_by_target[requirement.formatter_target] = (
            detected_type, requirement.consumer,
        )
        requested_type = getattr(args, "variant_id_type", None)
        if requested_type is not None and requested_type != detected_type:
            raise ValueError(
                "--variant-id-type %s conflicts with the %s BIM identifier type %s."
                % (requested_type, requirement.consumer, detected_type)
            )
        target_types[requirement.formatter_target] = detected_type
        observation["consumers"] = sorted(
            cached_consumers | {requirement.consumer}
        )
        references = dict(cached.get("references", {}))
        references[str(resolved_bim)] = dict(observation)
        observation["references"] = references
        observations[requirement.formatter_target] = observation

    args.variant_id_types = target_types
    args.variant_id_observations = observations
    args._variant_id_requirements = detected_by_target
    return observations


def configure_required_variant_identifier_type(
    args,
    formatting_config,
    *,
    consumer: str,
    formatter_target: str,
    required_type: str,
) -> None:
    """Set one protocol-required target identifier type or reject conflicts."""
    requested_type = getattr(args, "variant_id_type", None)
    if requested_type is not None and requested_type != required_type:
        raise ValueError(
            "--variant-id-type %s conflicts with the %s requirement for %s."
            % (requested_type, consumer, required_type)
        )
    target_types = dict(
        getattr(args, "variant_id_types", None)
        or formatting_config.variant_identifiers.target_types
    )
    configured_type = target_types.get(formatter_target)
    if configured_type is not None and configured_type != required_type:
        raise ValueError(
            "%s requires formatter target %s to use %s, but configuration "
            "requests %s."
            % (
                consumer, formatter_target, required_type, configured_type,
            )
        )
    target_types[formatter_target] = required_type
    args.variant_id_types = target_types
    requirements = dict(getattr(args, "_variant_id_requirements", None) or {})
    requirements[formatter_target] = (required_type, consumer)
    args._variant_id_requirements = requirements


__all__ = [
    "BimIdentifierRequirement", "configure_reference_variant_identifiers",
    "configure_required_variant_identifier_type",
]
