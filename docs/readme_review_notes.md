# README review notes

Findings from reviewing the PostGWAS README against the implementation
(2026-08-10). The README was rewritten from the code; this file records what was
**deliberately left out or worded cautiously**, and the code/doc conflicts that
need a decision. No code was changed.

Line references come from static reading of the current working tree, not from
executing the pipeline.

## 1. Behaviour that would break a documented workflow

These are not documentation problems, but they determine what the README can
honestly promise.

| # | Issue | Location |
|---|---|---|
| 1 | `imputation` reads `args.corr_method`, but the flag is `--correlation-method` with no `dest=`, so argparse produces `args.correlation_method`. No parser defines `corr_method`. Post-processing would raise `AttributeError` in both direct and pipeline mode. | `modules/imputation/service.py:48`, `cli/common.py:1244` |
| 2 | `postgwas qc` in direct mode reads `args.bcftools`, which its own parser never defines; the registry adds that parser only for pipeline mode. | `modules/qc_summary/service.py:15`, `modules/qc_summary/cli.py:103-108`, `pipeline/registry.py:366` |
| 3 | The imputation runner unpacks the harmonisation result one level too high (`"GRCh37" in outputs`), but harmonisation returns `{dataset_id: {"GRCh37": …}}`. | `pipeline/runners.py:334-335`, `modules/harmonisation/cli.py:1157` |
| 4 | `run_ldsc_direct` and `run_qc_summary_direct` never return their result dict, so `ctx["heritability"]` and `ctx["qc_summary"]` are always `None`. | `modules/ldsc/service.py:12`, `modules/qc_summary/service.py:8-16` |
| 5 | `ld_clump` wraps both clumping methods in bare `except Exception` that only prints a warning; a step can "succeed" having produced nothing. | `modules/ld_clumping/service.py:32-33, 51-52` |

Because of 1–3 the README documents `--apply-imputation` and `postgwas qc` as
existing interfaces without asserting they are validated end-to-end.

## 2. Documentation that conflicts with the code

| # | Issue | Suggested fix |
|---|---|---|
| 6 | `docs/wiki/reference/output-structure.md` never mentions the `NN_<module>` step numbering that `pipeline/runners.py:49-69` actually creates. | Add the numbering to the wiki page (the README now documents it). |
| 7 | Resolved: harmonisation now writes each dataset screen report, records its path, and uses the shared `--show-screen` / `--hide-screen` controls. Every analysis command also writes the configured shared screen transcript. | Keep the README and logging documentation aligned with the shared controls. |
| 8 | `--comparison-af-source` help renders `resource_examples` such as `GRCh37_1000G_freq_chr[1..22,X,Y].vcf.gz`, omitting the `{build}/external_af/vcf_files/` prefix that `resource_layout` actually requires. | Align the examples with `resource_layout`. |
| 9 | `examples/configs/harmonisation/run_config.yaml` sets `run.dataset_id`, `run.overwrite`, `run.resume`; harmonisation reads only `run.output_directory`. | Trim the example. |
| 10 | `--fixed-info` is CLI-only — setting `fixed_info.value` in a run-config is a hard error, though the key exists in the config model and packaged defaults. | Documented in the README as a flag only. |
| 11 | The sample-sheet `delimiter` column does not control parsing; the reader resolves the separator from the `input.delimiter` policy and uses the sheet value only for a cross-check whose message claims the declared delimiter "remains authoritative". | Decide which wins, then fix the message. |
| 12 | Compressed-input support is inconsistent: `open_text` and the sample-sheet generator accept `.bz2` / `.xz` / `.zip`, but the bulk read path is `polars.read_csv` on the raw path. | Verify before advertising those formats. The README says only that the sample sheet declares the columns. |

## 3. Settings that look configurable but are inert

Documenting these as user knobs would be wrong, so the README does not mention
them.

- `pipeline.modules` and `pipeline.stop_after` are validated but never drive
  execution; module selection comes only from `--modules`
  (`config/models/application.py:70-76` vs `pipeline/cli.py:80,129-145`). The
  README states this explicitly because `config export --pipeline` writes a
  `pipeline:` block that a user would reasonably expect to work.
- Module YAML is never loaded at runtime for `imputation`, `ld_annotation`,
  `ld_clumping`, `qc_summary` and `manhattan`; those modules read argparse
  attributes instead. Affected keys include `ld_clumping.{lead_pvalue, clump_r2,
  lead_r2, window_kb, merge_distance_bp, remove_mhc}` and
  `manhattan.{significance_threshold, suggestive_threshold, file_format}`.
- Filtering defaults disagree between layers: `filtering.yaml` has
  `maf_min: 0.01`, `info_min: 0.7`, `remove_mhc: true`; the shared argparse
  parser has `0.002`, `0.3`, `false`. Pipeline mode uses the argparse values,
  `postgwas sumstat_filter` uses the YAML values. **No filtering thresholds are
  quoted in the README for this reason.**
- `--threads`, `--memory-gb` and `--seed` are resolved by
  `resolve_compute_args` from `load_configuration()` *without* the user's
  `--run-config`, then passed downstream as explicit CLI values — so
  `execution.threads` / `memory_gb` / `random_seed` set in a run-config are
  overridden by hardware detection (`cli/compute.py:74-85`, `pipeline/cli.py:352`).
  The README describes precedence as defaults → YAML → CLI, which is the
  documented intent; this path contradicts it.
- `modules.enrichment` (providers, thresholds) is never read — all nine
  providers always run. `modules.ldsc` is read only by the single-cell
  `ldsc_celltype` runner, not by `heritability`.
- `modules.allele_orientation` is a required configuration block that no runtime
  code reads.
- `resources.executables.{finemap, ldstore, bgenix, ldsc, munge_sumstats}` are
  declared but ignored; those tools are resolved by `shutil.which` or hardcoded
  names. Only `magma`, `gcta`, `plink`, `rscript`, `bcftools` and `mixer` are
  honoured.
- The packaged profiles `minimal`, `standard`, `comprehensive` are unreachable
  from any CLI (`profile=` is a Python-API keyword; `tests/test_configuration.py:259-267`
  asserts `--profile` is absent). Not documented.
- `logging.capture_external_tools`, `logging.include_timestamps`,
  `execution.fail_fast`, `pipeline.fail_fast` have no read sites.

## 4. Method and scope points worth an explicit decision

- `ld_clump` region-based clumping hardcodes significance at `LP >= 7.3`
  (`ld_prune_region.py:198`) instead of using `--lead-p`, so the two methods in
  the same run use different significance definitions.
- The module registers a PLINK binary parser but never invokes PLINK; clumping
  is bcftools + tabix over precomputed `.ld.gz` tables (`pipeline/registry.py:166`).
- Two MHC windows exist in the codebase (`finemap` defaults 6:25–35 Mb versus
  `pops.py` 6:20–40 Mb).
- `single_cell`'s static registry dependency is `("magma",)`, which is wrong for
  `ldsc_celltype` (needs `formatter`). Correct only because `pipeline/cli.py`
  passes `dependency_overrides`.
- scDRS produces no normalised `results/*.tsv`, unlike `magma_celltype` and
  `ldsc_celltype`, so multiple-testing normalisation is not uniform across the
  three single-cell tools.
- `gcta_cojo`, `gcta_gene`, `mixer`, `heritability` and `enrichment` are
  terminal: nothing consumes their results. Stated in the README under
  limitations.
- FINEMAP `prior_k` is exposed in the model and defaults but unconditionally
  rejected at runtime (`fine_mapping/service.py:175-182`).
- `pathway_enrichment` hardcodes a conda env name (`enricher`), an R zip binary
  path (`/opt/conda/envs/enricher/bin/zip`), `tax_id=9606`, and roughly 200
  Enrichr libraries with a 0.5 s sleep each. It also ignores `--reference-set`
  and reads two options (`string_score`, `fdr_thr`) that no parser defines.

## 5. Documentation coverage gaps

- There is **no user-guide page for `gcta_cojo`**, although it is a first-class
  command and pipeline target. `docs/wiki/modules/` has `gcta-gene.md` only, so
  the README module table lists a command the documentation index cannot link.
- `kpops` and `caldera` pages live at `docs/modules/*.md` rather than
  `docs/wiki/modules/`, unlike every other module page, and therefore do not
  follow the enforced module-page heading structure
  (`tests/test_wiki_docs.py:67-71` only checks `docs/wiki/modules/*.md`).
- No page documents the `allele_orientation` placeholder or the terminal-module
  boundary; both are now summarised in the README.

## 6. Repository hygiene (blocks a clean public release)

- `LICENSE` contains a `[tool.poetry.scripts]` snippet, not a licence. **The
  README therefore has no licence section** — add the intended licence text
  first.
- `pyproject.toml` still carries placeholder metadata: description "PostGWAS
  analysis toolkit"  with a `# Customize` comment, author `Jibin
  <johnjibinv@email.com>`, and a stale comment on `[project.scripts]`.
- Untracked run artefacts sit in the repository root: `ldsc.log`, `magma.log`,
  `.DS_Store`, `output/`, `results/`.
- Dockerfile reproducibility gaps: `snpEff_latest_core.zip` and
  `remotes::install_github('jianyanglab/gsmr2')` are unpinned; GCTA is installed
  twice (pinned 1.95.0 at `/usr/local/bin/gcta`, plus conda `gcta`), and with the
  conda env first on `PATH` the default `executables.gcta: gcta64` resolves to
  the conda build, so the pinned binary is unused.
- Python version drift: `requires-python >= 3.10`, conda env pins 3.10, docs CI
  uses 3.11.
- Default resource path contains a typo — `.../software_resources/resourses/...`
  — and it is echoed in `postgwas resources --help`
  (`resources/magma_functional_mapping/config.yaml:2`).

## 7. Constraints the README must keep satisfying

`tests/test_wiki_docs.py` enforces the following on the root README, and CI runs
`tools/docs/validate_wiki_cli.py` against it:

- a `## Documentation` heading and the phrase "GitHub Wiki access is not enabled";
- a link to **every** page listed in `docs/wiki.yml`, and the **first**
  occurrence of those links must follow manifest order — so in-body links to
  wiki pages will break the test unless they form a prefix of that order;
- the literal string `--sample-sheet studies.csv`;
- absence of `JIBINJOHNV/postgwas.git`, `JIBINJOHNV/postgwas_v2/wiki`,
  `jibinjv/postgwas:1.3`, `sumstat_file`;
- every `postgwas …` command in a fenced block must use real commands, real long
  options and real choice values for the installed CLI.

Verified after the rewrite: `pytest tests/test_wiki_docs.py` → 14 passed;
`tools/docs/validate_wiki_cli.py` → 140 documented commands validated.
