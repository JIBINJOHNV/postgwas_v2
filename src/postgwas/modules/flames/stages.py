"""Ordered scientific progress contract for the PostGWAS FLAMES wrapper."""

FLAMES_STAGES = (
    "Validate the fine-mapping, MAGMA, MAGMAcovar, and PoPS inputs",
    "Validate the FLAMES model, annotation resources, and runtime",
    "Run and validate FLAMES locus annotation",
    "Run and validate FLAMES gene scoring",
    "Validate and publish all FLAMES outputs",
)


__all__ = ["FLAMES_STAGES"]
