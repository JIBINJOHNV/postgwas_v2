# dbSNP Build 157 GRCh37 resource preparation

This workflow retains the official NCBI files unchanged and produces three
analysis-ready derived resources:

1. A GRCh37.p13 dbSNP VCF restricted to primary assembled molecules and renamed
   from RefSeq accessions to `1`–`22`, `X`, `Y`, and `MT`.
2. A gzip-compressed previous-to-merged-into rsID table derived from
   `refsnp-merged.json.bz2`, including older `dbsnp1_merges` members.
3. A gzip-compressed withdrawn-rsID status table, also carrying historical IDs
   that had merged into a subsequently withdrawn RefSNP.

The merged table columns are `PREVIOUS_RSID`, `MERGED_INTO_RSID`, `RELATION`,
`MERGE_BUILD`, and `MERGE_DATE`. The withdrawn table columns are `RSID`,
`WITHDRAWN_ROOT_RSID`, `RELATION`, `WITHDRAWN_DATE`, and `LAST_UPDATE_BUILD`.
Most merged records provide one or more current rsIDs in `MERGED_INTO_RSID`.
NCBI also distributes rare merged records whose official `merged_into` array is
empty. Those old IDs are retained with a blank `MERGED_INTO_RSID` and an
explicit `*_without_current_target` relation. Observation-movement rsIDs are
not substituted for an absent official merge target because they describe
allele remapping rather than an asserted RefSNP merge relationship.

The original VCF, Tabix index, merged JSON, withdrawn JSON, upstream MD5 files,
JSON README, and GRCh37.p13 assembly report remain under `downloads/`.

## Chromosome policy

The mapping is derived from NCBI's assembly report. Only records whose
`Sequence-Role` is explicitly configured as `assembled-molecule` are retained.
Non-primary alt, patch, unlocalized, and unplaced records are excluded from the
derived VCF, with their record count logged and stored in the manifest. They
remain available in the unchanged source VCF.

Coordinates, REF, ALT, FILTER, and current rsIDs in the `ID` column are retained
without normalization. All INFO annotations and their header definitions are
removed from the derived VCF, so its INFO column is `.`. The unchanged source
VCF retains every original INFO annotation. The source VCF header must identify
dbSNP Build 157 and GRCh37.p13 before conversion begins.

## Run

```bash
python tools/resource_preparation/dbsnp_build157_grch37/prepare.py \
  --config tools/resource_preparation/dbsnp_build157_grch37/config.yaml \
  --stage all
```

Valid stages are `download`, `transform`, `validate`, and `all`. Downloads are
resumable, atomic, size-checked, and verified against NCBI's published MD5s.
Derived outputs are written atomically and indexed with CSI.

## Official sources

- https://ftp.ncbi.nih.gov/snp/latest_release/VCF/
- https://ftp.ncbi.nih.gov/snp/latest_release/JSON/
- https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/328
- https://ftp.ncbi.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt
