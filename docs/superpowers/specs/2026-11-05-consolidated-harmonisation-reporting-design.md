# Consolidated harmonisation reporting design

Status: implemented and locally verified; publication to main explicitly authorized. Exact server execution remains unverified.

## Evidence and scope

Reviewed implementation: `/Users/JJOHN41/Downloads/gwas_sumstat/.implementation/postgwas_v2`, commit `cc2bd0718ac7c4e53452041c6408bb63fa564157`.

User transcript: `/Users/JJOHN41/.codex/attachments/57ece1cd-a6d6-420c-9613-950bb14ec5fa/Pasted text.txt`.

The supplied excerpt contains two complete study-decision reports, two selected-build fields, and 44 chromosome build rows: all 22 chromosomes printed twice. Effect and P-value evidence first appears in automatic-inference warnings and then in both reports. Input completeness and the removal ledger are separate blocks with overlapping presentation. Current native-AF reporting additionally repeats screen counts/cutoff inside one subsection.

Verified causes: `service.py:3401` prints the full report before external AF approval; `_report_study_failure()` at `service.py:506` prints it again on the subsequent AF failure. `study_reporting.py:106` expands every chromosome on each render. Routine inference warning cards originate at `study_properties.py:451,714`. The actual production call to `infer_genome_build()` occurs once at `service.py:3127`; chromosome evidence is aggregated from existing matches.

## Required readable output

1. **Input validation:** one section combining dataset/source identity, effective column mappings, missing values, conversion results and the removal ledger. Refer to established file numbers rather than repeating long source inventories. If a template/source has not yet been registered, show its exact declaration once; never invent a file number or read a file only for presentation. Keep original input and post-filter conversion denominators explicit.
2. **Study-wide decisions and evidence:** one detailed report after dataset preparation, AF approval and INFO resolution have reached a terminal preparation outcome, immediately before worker launch. Include build, effect, SE, P-value, BETA/Z, strand, AF A/B and INFO decisions together. Preparation approval is not final analysis/output success.
3. **Genome build:** one dataset decision plus pooled supporting evidence and action. No routine chromosome-by-chromosome build rows in readable terminal, screen or dataset logs. A blocking conflict may name affected chromosomes concisely in its single dataset failure reason.
4. **Failures:** if preparation fails before the report, show one report containing only measured/known evidence, the failed check, its counts/thresholds and the stop reason. Never show unperformed checks as passed. Later failures report the failing stage/new findings without replaying unchanged study decisions.
5. **Warnings:** routine automatic effect/P-value notices appear once inside the consolidated report in orchestrated runs. Preserve their warning severity, reasons, recommendations and stored evidence. Declaration conflicts and fatal warnings remain immediate. Direct library calls retain current visible notices by default.
6. **Progress and final results:** retain mandatory stage start/completion/failure messages, chromosome processing progress, and final output/completeness summaries. These describe different stages; they must not be suppressed as duplicate scientific reports.

## Scientific and audit boundaries

- Change reporting ownership/timing, not scientific calculations, source selection, filtering, thresholds, AF approvals, stage order, chromosome partitioning or checkpoint identity.
- Preserve internal chromosome build evidence, distinct loci, normal/fallback gates and existing fallback-only chromosome veto behavior. Full evidence remains in existing structured stage extras/manifests, including failed-run manifests; routine readable per-chromosome observations are removed.
- Preserve all supplied external AF identities and absent-native-column status. Show A, conditional diagnostic and B together; do not publish a pending selection as an approved final source.
- Never sum unrelated field-missing counts into removed-row counts. In the supplied case, missingness refers to 11,749,765 input rows; conversion/inference refers to 11,434,316 retained rows; 315,449 is a removal count, not an additional loss every time it is mentioned.
- Reuse existing renderers, `_announce()`, `_report_study_failure()`, shared UI/progress, screen recorder and file-validation references. Do not add text matching, global suppression flags, a second logger or a parallel reporting framework. Presentation state must not enter scientific decisions or approval digests.
- No new scientific configuration/defaults or dependencies; protected `harmonisation/adapters/` remains untouched.

## Acceptance

On successful or failed preparation: one detailed study report, one selected-build field if known, no expanded chromosome build listing, no duplicate inference cards, one combined input-QC section, and no duplicate native AF screen values. Visible and hidden transcripts and the canonical dataset log retain the same known results and warnings. Scientific frames, counts, source contracts and build decisions remain unchanged; rendering performs no file rereads, joins or inference. Multiple datasets receive independent reports.

This design does not establish the exact server executable version or claim the pasted server run has been rerun.
