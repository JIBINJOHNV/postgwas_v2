# Fine-mapping scientific corrections

The public function names and existing return keys are preserved. Both engines now use a fixed 95% credible-set target.

| Area | Implemented correction | Scientific reason |
|---|---|---|
| SuSiE credible sets | Every `susie_rss` call, including recovery calls, explicitly uses `coverage = 0.95`; exported component coverage is checked. | Recovery must not silently change the inferential target. |
| SuSiE PIP export | FLAMES `prob1` is the SuSiE global PIP for each member of the component-defined 95% credible set. | Component `alpha` defines membership; global PIP is the SNP-level inclusion probability expected by FLAMES. |
| SuSiE mapping | Model indices are restored before mapping PIP, alpha, effect moments, and log Bayes factors to variants. | Sorting or merging without index restoration can attach posterior values to the wrong SNP. |
| Alleles and LD | Summary alleles are aligned to the allele counted by the LD matrix; z-score/beta signs are flipped on swaps. Strand-ambiguous A/T and C/G SNPs are excluded without frequency evidence. | Fine-mapping is invalid when association signs and LD signs use different effect alleles. |
| Sample size | Invalid or missing effective sample size causes a locus failure; the old `10,000` fallback was removed. | An invented sample size changes Bayes factors and posterior probabilities. |
| LD validation | Variant order, coordinates, alleles, dimensions, finite values, diagonal, symmetry, correlation bounds, eigenvalues, and LD/z consistency are checked. | FINEMAP and SuSiE require the z vector and LD matrix to describe the same ordered, allele-coded variants. |
| LD repair | SuSiE recovery records the repair and rejects a repair whose maximum correlation change exceeds 0.05. FINEMAP does not silently repair LD. | Large repairs can manufacture a materially different likelihood. |
| FINEMAP credible sets | `--prob-cred-set` is always `0.95`; the selected causal-count model is deterministic, and each component is exported separately. | The model posterior and SNP PIP are different quantities and must not be interchanged. |
| Defaults and audit trail | Scientific/resource controls are named in `defaults.py` or `susie/defaults.r`; each run writes effective configuration, stage logs, software versions, percentage completed, and work remaining. | Results must be reproducible and failures must not be hidden behind console-only messages or unexplained literals. |
| SuSiE probability semantics | Credible-set membership is the 95% single-effect component set (`alpha`), while FLAMES `prob1` is the model-wide SNP PIP required by FLAMES. Both definitions and both coverage checks are recorded. | Component posterior weights and model-wide SNP PIPs answer different questions; silently treating them as interchangeable is scientifically misleading. |
| FINEMAP-to-FLAMES IDs | Credible-set rsIDs are mapped through the harmonized `.z` file to `CHR:BP:A1_A2`; IDs are no longer guessed by replacing underscores. | An rsID does not intrinsically encode chromosome, position, or allele orientation. |
| FLAMES files | Each file contains only `index`, `cred1`, and `prob1`; `indexfile.txt` contains `Filename`, `GenomicLocus`, and `Annotfiles`; a manifest records coverage/build/model metadata. | This follows the FLAMES example and the minimal formatter used by FUMA, without attaching misleading model-level comments to SNP-level columns. |
| Failure handling | Zero-task/all-failed runs raise errors, recovery placeholders are replaced in place, and only current successful loci enter FLAMES indexing. | A pipeline must not report success or merge stale output when no valid locus was produced. |
| Reproducibility | Software versions, input harmonization QC, per-locus status, LD QC, failure reasons, and genome build are written. | These fields are needed to reproduce and audit fine-mapping results. |

Reference implementations reviewed:

- FLAMES example credible set: <https://github.com/Marijn-Schipper/FLAMES/blob/master/example_data/locus_files/locus_1.cred1>
- FLAMES example index: <https://github.com/Marijn-Schipper/FLAMES/blob/master/example_data/indexfile.txt>
- FLAMES input documentation: <https://github.com/Marijn-Schipper/FLAMES/blob/master/README.md>
- FUMA FLAMES formatter: <https://github.com/vufuma/FUMA-webapp/blob/master/scripts/flames/format_finemapping_out.py>
- FUMA FLAMES workflow: <https://github.com/vufuma/FUMA-webapp/blob/master/scripts/flames/run_flames.py>

The included unit tests verify exact FLAMES columns and IDs, reject a cumulative FINEMAP credible set below 95%, and protect the existing return-key contracts.
