# Unified GWAS frequency screening and replacement acceptance

Status: approved specification; implemented locally on codex/unified-af-screening, with final verification recorded separately.

## Agreed behavior

The selected frequency, A, is either the mapped GWAS frequency column or the frequency merged from the user's initially selected external file. Preserve A and its source provenance. Existing source-selection precedence remains unchanged. Screen frequencies in the direction of the original listed GWAS effect allele, before later canonical ALT/effect transformations.

Apply one common MAF-like screen to native A, externally merged A, and a subsequently selected replacement B. Screen the GWAS rows, not the full external panel. A row contributes once to the screen; do not average chromosome percentages or filter the screen to the later diagnostic subset.

Usable values are finite numeric frequencies in [0,1]. Calculate p = count(usable AF <= 0.5) / count(usable AF). The result is MAF-like only when p is strictly greater than the configured eaf.maf_decision_cutoff (current default 0.95). Equality does not trigger the diagnostic. Zero usable values stop the dataset; they must not be treated as an acceptable non-MAF-like result.

## Initial frequency A

1. Preserve A. If its screen is not MAF-like, retain it and proceed through remaining normal QC. Skip only the additional semantic default-reference comparison; reference use for strand and other QC remains unchanged.
2. If A looks MAF-like, merge the configured default-panel frequency onto the same original GWAS effect allele. Let that effect-aligned value be f.
3. Use the requested diagnostic convention: other_AF = 1-f and folded_AF = min(f,1-f). Do not change A during comparison. This is effect-allele-versus-all-other-alleles folding at multiallelic sites, not a claim about a universally defined multiallelic site MAF.
4. Do not exclude diagnostic evidence solely because the site is multiallelic. Preserve the other existing protections: compatible build/allele representation, unique and independently resolved orientation, palindrome safeguards, valid mapped frequencies, and duplicate/conflict handling.
5. Preserve the existing absolute diagnostic count, fixed-N coverage, correlation, error and profile thresholds. Keep N as all study rows entering the initial AF preflight; diagnostic exclusions do not reduce N. No automatic change from 0.50 to 0.45 is included.

Let x=A, r=f, and m=min(r,1-r). The current defaults remain:

| Requirement | Rule 1: retain A as EAF | Rule 2: MAF support |
|---|---|---|
| Shared usable pairs | >= 1,000 | >= 1,000 |
| Shared fixed-N coverage | >= 0.50 | >= 0.50 |
| Correlation | corr(x,r) >= 0.95 | corr(x,m) >= 0.95 |
| Mean absolute difference | mean(abs(x-r)) <= 0.10 | mean(abs(x-m)) <= 0.10 |
| Reference profile | proportion(r <= 0.5) > 0.90 | proportion(r > 0.5) > 0.20 |

Use configured values rather than independent hardcoded defaults. Numeric equality passes count, coverage, correlation and error; profile equality fails. Pearson remains the default and configured Spearman remains supported. Both rules are evaluated after the shared gates. Only rule 1 passing retains A. Only rule 2 passing selects the explicit replacement. Both passing, neither passing, or failed shared evidence stops the dataset under the existing failure policy. Exact numerical agreement is not required.

Preserve raw effect-aligned r and B alongside their folded views. Do not overwrite r or B with min(f,1-f): the unchanged minor/major-profile rules must still see the raw aligned values. In particular, globally folding r would force its major fraction to zero and prevent rule 2 from passing.

The existing requirement for independent default evidence when classifying an initially selected external source remains. The same physical source or its aliases cannot be treated as independent confirmation. This requirement is distinct from existing labelled same-source reuse in the replacement consistency check.

## Replacement B

When the initial diagnostic supports MAF, require both --eaf-replacement-file and --eaf-replacement-column. The source may now follow either native A or externally merged A. This extends the current code, which stops rather than selecting replacement after initial external MAF support.

Merge B through the same existing external-frequency machinery, aligned to the original listed GWAS effect allele. Freeze A rather than overwriting it. Run the common screen on the resulting B column across the GWAS rows before adoption or further replacement validation.

- No usable B: stop the affected dataset.
- B looks MAF-like: reject B and stop the affected dataset. Do not run another default semantic comparison, select another replacement, retry, or fall back to A/default frequencies. This is the user's explicit replacement acceptance rule, even if another concordance statistic could pass.
- B does not look MAF-like: preserve B and apply the remaining existing QC before final adoption.

Because all other rules are retained, the existing replacement consistency gate remains after B passes its screen: compare unchanged A with min(B,1-B), using the same original-EA direction. Preserve its current count/coverage, correlation >=0.95, error <=0.10, and proportion(B>0.5)>0.20 defaults and configured overrides. Its diagnostic eligibility must follow the same newly requested convention, without reintroducing a multiallelic-only exclusion. A failed consistency gate still stops the dataset.

## Boundaries

- The change concerns screening, diagnostic value alignment and replacement acceptance. Preserve the existing complete-site multiallelic attachment logic, final strand/effect transformations, output interfaces and protected adapters.
- Do not globally remove structural guards from strand orientation or export. Do not add ALT/ALT export support as part of this change.
- Implementation review exposed an unsupported native multiallelic EAF swap: complementing a study REF frequency would combine other ALT frequencies. Preserve diagnostic acceptance, but reject/fail that unsupported native output swap under the existing strand ambiguity policy and retain original values in rejection provenance. Do not substitute population AF. External/replacement complete-site assignments are unaffected.
- Preserve missing/invalid AF policies, duplicate policies, per-chromosome usable attachment coverage (default 0.80), and final output checks. A global diagnostic pass does not override these checks.
- Frequency source and storage location remain explicit. User-supplied external AF is a supplied frequency input; reference-population provenance must not be relabelled as study-measured frequency.
- This screening rule is a configured decision heuristic, not proof of source-column semantics.
- Preserve current thresholds, CLI options, builds, resources and input payloads. No cloud mutation, input rewrite, threshold relaxation or source switching without the defined transition is included.

## Required records

Record source file/column/mapping/identity; original GWAS EA direction; A and B roles; screen scope; total GWAS rows, usable and unusable values, count <=0.5, fraction, cutoff and result; fixed diagnostic N and n; exclusion reasons; default-panel comparison evidence; replacement screen/consistency results; exact stop reason; and whether final source adoption occurred. A failure must never be presented as approved EAF.

For externally supplied A, retain immutable private row-bound evidence and its descriptor through replacement validation. Do not fabricate a native eaf_col or change the fact that the original GWAS frequency column was absent. Downstream source-selection, concordance and VCF provenance consumers must follow the selected external/replacement source correctly.

## Approaches considered

The user selected the shared merge plus 1-effect_AF diagnostic convention. A full per-allele multiallelic map and an automatic lowering of diagnostic coverage were discussed earlier and are not part of this specification. Existing complete-site attachment remains reusable; diagnostic and attachment policies must remain separately scoped.

## Verification baseline and limits

The owning repository inspected is /Users/JJOHN41/Documents/developing_software/postgwas_v2 at commit 2cd5a98. Existing tests and local fetal-WLM chromosome 22 diagnostics were executed earlier in this chat, but they do not test this unimplemented flow. The available real extract is full fetal WLM, not the exact current server wlm_sex_gest_plus_p_fetal_wlm_pair.gz payload. Full server validation remains required.
