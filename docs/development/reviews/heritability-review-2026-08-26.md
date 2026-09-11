# PostGWAS `heritability` (LDSC) — module review

**Date:** 2026-08-26 · **Scope:** `src/postgwas/modules/ldsc/` (`service.py` 465, `ldsc_runner.py` 419, `cli.py` 162), `config/models/modules/ldsc.py`, `config/defaults/modules/ldsc.yaml`, `docs/wiki/modules/ldsc.md`, `tests/test_ldsc_runtime.py`, and the prevalence handoff from `formatter`.
**Method:** read-only on the tree — no file was modified. Every constructed flag was checked against the CBIIT `ldsc.py` / `munge_sumstats.py` command-line contract, and every packaged default against LDSC's own default. Line numbers refer to the working tree as of 2026-08-26.

**Follow-up (2026-08-27):** I1 and I2 have been implemented. PostGWAS now
compares explicit pipeline sample prevalence with the formatter-returned value,
and it validates both the intercept and attenuation ratio across observed and
liability logs before publishing outputs. The screen and canonical result log
report those shared regression findings once and identify liability h² as a
conversion. I3 and the optional findings remain open.

**Summary.** This is a thin, disciplined wrapper around a pinned external tool, and it is written the way such a wrapper should be: every flag maps to a real LDSC option, `null` means *omit the flag* so LDSC keeps its own data-dependent behaviour, results are parsed and validated rather than trusted, outputs are staged and published only after verification, and a test asserts the packaged defaults still match the pinned upstream. There are no Critical findings.

Three things are worth attention: the LD-score regression is **run twice** to obtain a value that is a closed-form scalar multiple of one already computed; a user-supplied sample prevalence **overrides the data-derived one with no comparison**, which is the one path to a quietly wrong liability heritability; and the module with the longest external steps is the one with **no resume**.

---

## IMPORTANT

### I1 — The LD-score regression runs twice for a result that differs by a scalar

**Where:** `ldsc_runner.py:345-364` (observed) and `:381-402` (liability).

When a population prevalence is configured, `run_ldsc` issues a second `ldsc.py --h2` over the *same* munged sumstats, the *same* `--ref-ld-chr` and `--w-ld-chr`, the *same* `--n-blocks`, `--chisq-max`, `--two-step` and `--intercept-h2` — differing only by `--samp-prev` and `--pop-prev`. LDSC applies the liability conversion as a post-hoc scalar to h² and its standard error; the intercept, the attenuation ratio and the underlying jackknife regression are unchanged. The second run therefore re-reads 22 chromosomes of LD scores and repeats a 200-block jackknife in order to multiply by

```
K²(1 − K)² / [ P(1 − P) · φ(z)² ]
```

**The counter-argument is legitimate and should be weighed.** Letting LDSC own the conversion avoids reimplementing a formula whose φ term and threshold are easy to get subtly wrong, and the cost is minutes rather than hours — `munge_sumstats` over a genome-wide file is often the slower step anyway. Reasonable people ship this deliberately.

**What is not a judgement call** is that nothing in the output says the two are one analysis. `_record_findings` (`service.py:153-172`) logs `h2`, `intercept` and `ratio` from *both* logs, and the terminal summary prints both blocks. A reader sees two intercepts and two ratios and may reasonably take them as two estimates that happen to agree.

The cheap improvement keeps LDSC as the authority and costs one comparison: **assert that the two logs' intercepts match**, and present the liability block as a *conversion of* the observed estimate rather than as a parallel result. That assertion also converts the claim in this finding into something the project tests rather than assumes — and it is the kind of check that would catch a future LDSC version changing the conversion's scope.

### I2 — A user-supplied sample prevalence overrides the data-derived one with no comparison

**Where:** `service.py:101-113` (`_resolve_prevalence`); the value it declines to consult is produced at `formatting/exporters/ldsc.py:57-97`.

```python
source = "configuration_or_cli" if module.sample_prevalence is not None else "unset"
if module.population_prevalence is None:
    return source
if module.sample_prevalence is None:
    module.sample_prevalence = _formatter_sample_prevalence(ctx)
    source = "formatter_return_value"
```

The formatter's value is read **only** when the user supplied none. This precedence is documented (`ldsc.md`: "An explicit `--samp-prev` overrides the formatter return"), the source is recorded, and `_prevalence_source_label` renders it in plain language on screen — the provenance handling here is genuinely good.

What is missing is the check. The formatter has already computed the case fraction from the actual per-variant counts and returns it *together with* its minimum, maximum, and the case- and control-count ranges. A user who types `--samp-prev 0.5` out of habit — a very common habit, because many people think of it as the ascertainment design rather than as a property of the file — for a study whose real case fraction is 0.11 gets a liability h² wrong by roughly a factor of two, and nothing in the run mentions the disagreement.

Reading the formatter's value even when an override is present, and warning when the two differ materially, uses data the module already receives and costs one comparison. This is the single highest-value change in the module: it is the only path here to a quietly wrong scientific number.

### I3 — No completion manifest and no resume, in the module with the longest external steps

`run_ldsc_direct` contains none of the `core.completion` machinery that `qc_summary` and `formatter` use — no configuration digest, no input/output fingerprints, no validated resume. It detects existing owned outputs and refuses (`service.py:392-401`), which is safe, but a re-run then repeats `munge_sumstats` over the whole genome-wide file and both `--h2` regressions from scratch.

The foundation is already there: `_publish_staged_result` (`:117-152`) runs everything into a temporary directory, verifies each expected file is present and non-empty, refuses to publish any path outside the staged prefix, and only then moves them into place. That is exactly the boundary a completion manifest hangs on.

---

## OPTIONAL

**O1 — An all-missing INFO column becomes `NA` and may empty the munge step.** The formatter's LDSC export lists `EAF: FRQ` and `INFO: INFO` as `trailing_columns`, not `required_columns` (`formatting.yaml:258, 262`), so a study whose harmonised VCF carries no `FORMAT/SI` still gets an `INFO` column filled with the configured null token. This module then always passes `--info-min 0.9`. If the pinned `munge_sumstats` treats a non-numeric INFO as failing that threshold, every variant is dropped and the failure surfaces as a munge error that does not name the cause. Many published GWAS do not report imputation quality, so this is not an exotic input — worth confirming against the pinned commit, and if it holds, omitting the column (or the flag) when the column is entirely null turns a confusing failure into a clear one.

**O2 — HapMap3 restriction is applied twice.** The formatter already filters to `--merge-alleles` rsIDs with compatible allele pairs, deliberately, so that its duplicate policy operates on the reference-matched set (`formatting.yaml:118-123`). `munge_sumstats` then applies `--merge-alleles` again over the same file. This is harmless — the second pass is LDSC's own contract and cannot be dropped — but it deserves a line in the docs so nobody removes the formatter pass believing it to be the redundant one.

**O3 — `ldsc_runner.py` hosts a builder this module never calls.** `build_h2_cts_command` (`:132-150`) is used only by `single_cell/methods/ldsc_celltype/runner.py:573`. Sharing the command builder is right; hosting it in a module that does not use it makes the dependency invisible from both sides. A shared `ldsc_commands` home would make the coupling explicit.

**O4 — the liability skip is logged only when a `.record`-capable logger is present** (`ldsc_runner.py:407-411`). Trivial, but this module is otherwise scrupulous about recording every decision, and "we did not run the liability analysis, and why" is one worth keeping unconditional.

---

## Verified correct — do not re-flag

**Every constructed flag is a real LDSC option, used as documented.** `--h2`, `--ref-ld-chr`, `--w-ld-chr`, `--out`, `--n-blocks`, `--intercept-h2`, `--two-step`, `--chisq-max`, `--not-M-5-50`, `--print-cov`, `--print-delete-vals`, `--samp-prev`/`--pop-prev`; and for munging, `--sumstats`, `--out`, `--merge-alleles`, `--info-min`, `--maf-min`, `--chunksize`, `--n-min`, `--keep-maf`.

**The packaged defaults are LDSC's own defaults**, not invented ones: `--info-min 0.9`, `--maf-min 0.01`, `--chunksize 5000000`, `--n-blocks 200`, and `.l2.M_5_50` unless `--not-M-5-50`. `tests/test_ldsc_runtime.py:122` asserts the expected packaged values, the command-builder regression test asserts that data-dependent options remain omitted, and the Docker regression test protects the exact upstream commit. Together, these checks guard the configuration-to-tool boundary; the defaults test does not dynamically inspect upstream source.

**`null` means "omit the flag", so LDSC keeps its data-dependent behaviour.** The YAML says so at the top, and it matters for exactly the settings where LDSC's default is computed rather than fixed. For PostGWAS input with per-variant `N`, the pinned `munge_sumstats.py` default for `--n-min` is the pandas-default 90th percentile of `N` divided by 1.5, not 0.67 times maximum `N`. In single-annotation h² with an unconstrained intercept, pinned LDSC supplies a two-step cutoff of 30. In multi-annotation h², it instead supplies `--chisq-max` behaviour equivalent to `max(0.001 × max(N), 80)` when no explicit maximum was provided. Hard-coding any of those in PostGWAS would duplicate conditional logic owned by the pinned tool. See the pinned [`process_n`](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py#L308-L324) and [`estimate_h2`](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/sumstats.py#L310-L354) implementations.

**`intercept` and `two_step` are mutually exclusive**, enforced twice — in the pydantic model (`ldsc.py:138-143`) and again in `build_h2_command` (`:100-103`) — mirroring LDSC's own constraint. `--samp-prev` and `--pop-prev` are all-or-nothing (`:96-99`), and a sample prevalence supplied alone correctly does *not* trigger a liability run.

**Reference validation matches what LDSC actually reads.** `_chromosome_prefix` appends the path separator so LDSC's `1..22` reader resolves `<dir>/1.l2.ldscore.gz`, and `validate_reference_files` checks precisely those files for all 22 autosomes — including the subtlety that the `.M`/`.M_5_50` file is required for the **reference** LD scores but *not* for the **weights**, which is exactly LDSC's requirement. The M-suffix checked follows the `use_m_5_50` setting, so switching to `--not-M-5-50` validates `.l2.M` instead. `chromosomes` is pinned to `Literal[22]` with a docstring explaining that this is an upstream protocol invariant rather than a tunable.

**Results are parsed, not trusted.** `extract_ldsc_metrics` requires both an h² and an intercept, validates the `estimate (se)` form with a strict regex, rejects non-finite values and negative standard errors, accepts `constrained to X` when the intercept was fixed, and preserves LDSC's `Ratio < 0` message rather than coercing it into a number. `run_checked_command` receives `expected_outputs`, so a zero exit that produced no log is rejected. Six tests cover these paths individually.

**Publication is atomic and guarded.** Everything runs into a `TemporaryDirectory` under a staging root; `_publish_staged_result` verifies each staged file is present and non-empty, refuses any path that does not sit under the staged prefix, and only then moves them to their final names. Existing owned outputs are detected before any work starts, and the owned-path list conservatively includes the optional `.cov`/`.delete`/`.part_delete` files even when they were not requested.

**The prevalence handoff itself is well built.** `_formatter_sample_prevalence` reads the formatter's returned value and — per its own docstring — never derives one; it raises with a specific message when the formatter returned none, naming the trait type. The resolved source is logged and rendered. I2 is about adding a comparison, not about the mechanism.

**The documentation is honest about what cannot be checked.** It states plainly that observed-scale h² is not liability-scale h², that population prevalence is not the sample case fraction, and that standard LDSC LD-score files carry no machine-readable build or ancestry — so `population` and `genome_build` are recorded as declarations and the user must verify reference compatibility. Recording an unverifiable declaration *as* a declaration is the right call.

---

## Documentation

`docs/wiki/modules/ldsc.md` is accurate and its Limitations section is unusually candid. Two additions:

- Say that the observed and liability runs are the same regression, and that their intercept and ratio are therefore identical by construction (I1). Without it, the two output blocks read as independent estimates.
- Note that a user-supplied `--samp-prev` is used as given and is not checked against the case fraction the formatter computed from the data (I2), and that the formatter's value — plus its range — is available in the formatter's own report for comparison.

---

## Suggested order of work

1. **I2** — compare a user-supplied sample prevalence against the formatter's data-derived value and warn on material disagreement. Smallest change here, and the only one that protects a scientific number.
2. **I1** — assert the two logs' intercepts match and present the liability result as a conversion rather than a second estimate; then decide, with that assertion in place, whether the second regression is worth its cost.
3. **O1** — confirm `munge_sumstats`' behaviour on an all-`NA` INFO column against the pinned commit, and guard it if it drops everything.
4. **I3** — a completion manifest, hung on the existing staging boundary.
5. **O2–O4**, and the two documentation additions.

Nothing here changes a number the module reports today. Step 1 stops a user from reporting a wrong one.
