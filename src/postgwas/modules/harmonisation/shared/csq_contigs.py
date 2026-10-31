"""Explicit native CSQ chromosome lookup, without rewriting resource files.

bcftools removes a literal source prefix and prepends the destination prefix.
The protocol order is VCF, GFF, FASTA; output VCF labels stay unchanged.
https://samtools.github.io/bcftools/bcftools.html#csq
"""

from collections.abc import Iterable, Mapping

from postgwas.config.models.modules.harmonisation import HarmonisationCsqPrefixes


CSQ_CONTIG_ROLES = ("vcf", "gff", "fasta")
CSQ_PREFIX_CONFIG_KEY = "modules.harmonisation.vcf_processing.csq_prefixes_by_build"


def resolve_csq_prefixes(vcf_config: Mapping, build: str) -> dict[str, str | None]:
    """Use the canonical schema for both public and direct-library callers."""
    by_build = vcf_config.get("csq_prefixes_by_build")
    if not isinstance(by_build, Mapping) or build not in by_build:
        raise ValueError(
            "Required CSQ chromosome naming is not configured for %s. Set %s.%s "
            "with explicit vcf, gff and fasta prefixes (null for no prefix)."
            % (build, CSQ_PREFIX_CONFIG_KEY, build)
        )
    return HarmonisationCsqPrefixes.model_validate(by_build[build]).model_dump()


def csq_prefix_arguments(prefixes: Mapping[str, str | None]) -> list[str]:
    """Emit the native option only when the validated mapping needs it."""
    values = [prefixes[role] for role in CSQ_CONTIG_ROLES]
    if all(value is None for value in values):
        return []
    return ["--unify-chr-names", ",".join("-" if value is None else value for value in values)]


def csq_contig_lookup(labels: Iterable[str], prefix: str | None, role: str) -> dict[str, str]:
    """Reject distinct raw names that would collapse under a prefix mapping.

    Unrelated alternative contigs need not share a primary-chromosome prefix.
    Requested contigs are checked separately using native exact lookup rules.
    """
    result = {}
    for label in sorted(labels):
        logical = label[len(prefix):] if prefix and label.startswith(prefix) else label
        if not logical:
            raise ValueError("%s chromosome label %r becomes empty after prefix removal" % (role, label))
        if logical in result and result[logical] != label:
            raise ValueError(
                "Ambiguous %s chromosome names %r and %r collapse to %r under "
                "the configured prefix %r; do not silently select one."
                % (role, result[logical], label, logical, prefix)
            )
        result[logical] = label
    return result


def csq_target_contig(vcf_contig: str, prefixes: Mapping[str, str | None], role: str) -> str:
    """Mirror native source-prefix removal and destination-prefix addition."""
    source_prefix = prefixes["vcf"]
    if source_prefix and not vcf_contig.startswith(source_prefix):
        raise ValueError(
            "VCF chromosome %r does not start with the configured vcf prefix %r"
            % (vcf_contig, source_prefix)
        )
    logical = vcf_contig[len(source_prefix):] if source_prefix else vcf_contig
    if not logical:
        raise ValueError("VCF chromosome becomes empty after the configured prefix removal")
    return (prefixes[role] or "") + logical
