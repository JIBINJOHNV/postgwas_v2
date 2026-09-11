"""Small reusable pipeline-preflight evidence for isolated unit tests."""

from pathlib import Path

from postgwas.core.vcf import IndexedVcfValidation


def pipeline_input_vcf_evidence() -> dict[str, object]:
    vcf = Path("/validated/test.vcf.gz")
    return {
        "input_vcf": {
            "indexed": IndexedVcfValidation(
                vcf=vcf,
                bcftools="/validated/bcftools",
                header="##fileformat=VCFv4.2\n",
                genome_build="GRCh37",
                contigs=("1",),
                dataset_id="study",
                sample="study",
                variant_count=1,
                required_fields=(),
                size_bytes=1,
                modified_ns=1,
                index_identities=((Path(str(vcf) + ".tbi"), 1, 1),),
            ),
            "harmonised": {"genome_build": "GRCh37"},
        }
    }


__all__ = ["pipeline_input_vcf_evidence"]
