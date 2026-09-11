"""Public execution contract for the MAGMAcovar-only pipeline."""

from __future__ import annotations


# These MAGMA arguments control competitive gene-set analysis.  They are not
# part of a MAGMAcovar-only pipeline, whose MAGMA dependency produces only the
# gene-association result consumed by the gene-property analysis.
MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS = (
    ("gene_set_file", "--gene-set-file"),
    ("minimum_gene_id_overlap", "--minimum-gene-id-overlap"),
    ("gene_set_identifier_mismatch", "--gene-set-identifier-mismatch"),
    (
        "alternate_gene_id_duplicate_policy",
        "--alternate-gene-id-duplicate-policy",
    ),
    ("gene_location_alternate_id_type", "--gene-location-alternate-id-type"),
)


__all__ = ["MAGMACOVAR_EXCLUDED_PATHWAY_OPTIONS"]
