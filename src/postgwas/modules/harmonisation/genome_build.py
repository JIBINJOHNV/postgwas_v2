"""Dataset step 05 - infer the genome build of the input.

Implements the build rows of plan v3 section 6.4: an explicit ``build.mode``
instead of the ``--genome-build`` flag the current error message advertises but
that does not exist, absolute-count, chromosome-scoped reference-coverage and
testable-row evidence floors (``build.min_match_count``,
``build.min_reference_match_fraction`` and ``build.min_match_fraction``), and
optional deduplication of the references. Input-wide match fractions remain
reported as audit evidence but are not a decision threshold because the build
references are representative marker panels rather than variant inventories.

Build evidence follows the GWAS Catalog summary-statistics convention that a
variant occupies one row while either allele order is valid. Consequently, a
study row contributes at most one match to each candidate build, even if the
reference repeats a variant or contains both allele orientations. Sources:
https://www.ebi.ac.uk/gwas/docs/summary-statistics-format and
https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics.
"""

import polars as pl

from .shared.runtime import resolve_policies, step_context
from .shared.variant_columns import (
    canonicalize_variant_frame,
    has_canonical_variant_columns,
    palindromic_snp_expression,
    read_reference_variant_table,
    reverse_complement_expression,
)

__all__ = ["infer_genome_build", "GenomeBuildError", "POLICY_KEYS"]


POLICY_KEYS = [
    "build.mode",
    "build.confidence_ratio",
    "build.min_match_count",
    "build.min_reference_match_fraction",
    "build.min_match_fraction",
    "build.deduplicate_reference",
    "build.low_coverage_fallback_enabled",
    "build.low_coverage_min_match_fraction",
    "build.low_coverage_max_competitor_fraction",
    "build.low_coverage_min_unique_loci",
    "build.low_coverage_chromosome_min_competitor_loci",
    "build.low_coverage_chromosome_min_competitor_fraction",
    "chromosome.strip_chr_prefix",
    "chromosome.strip_leading_zero",
    "chromosome.rename_map",
]


class GenomeBuildError(RuntimeError):
    """The genome build cannot be inferred from the data supplied.

    A typed exception rather than sys.exit(), so the dataset-level retry can
    catch it.
    """


def _load_reference(
    path,
    label,
    step,
    deduplicate,
    column_mapping,
    policies,
):
    """Read one build reference and normalise its coordinates."""
    step.input(
        "Genome-build reference",
        build=label,
        file=str(path),
        delimiter=column_mapping.get("delimiter"),
        chromosome_column=column_mapping.get("chr"),
        position_column=column_mapping.get("pos"),
        reference_allele_column=column_mapping.get("a1"),
        alternate_allele_column=column_mapping.get("a2"),
    )
    raw, _detected = read_reference_variant_table(
        path,
        column_mapping,
        policies,
        error_type=GenomeBuildError,
        description="%s genome-build reference" % label,
    )
    configured_columns = [
        column_mapping["chr"],
        column_mapping["pos"],
        column_mapping["a1"],
        column_mapping["a2"],
    ]
    raw = raw.rename(dict(zip(
        configured_columns,
        ["__build_chr", "__build_pos", "__build_ref", "__build_alt"],
    )))
    reference, _normalization = canonicalize_variant_frame(
        raw,
        {
            "chr": "__build_chr",
            "pos": "__build_pos",
            "ea": "__build_alt",
            "oa": "__build_ref",
        },
        policies,
        error_type=GenomeBuildError,
        warn=step.warn,
        label="%s reference" % label,
    )

    if deduplicate:
        before = reference.height
        reference = reference.unique(
            subset=[
                "__build_chr",
                "__build_pos",
                "__build_ref",
                "__build_alt",
            ],
            keep="first",
        )
        step.info(
            "{} reference deduplicated: {:,} rows before, {:,} after ({:,} "
            "repeats removed).".format(
                label, before, reference.height, before - reference.height
            )
        )
    else:
        step.info("{} reference rows: {:,}.".format(label, reference.height))
    return reference


def infer_genome_build(
    t_df: pl.DataFrame,
    reference_files: dict[str, str],
    reference_column_mapping: dict[str, str],
    sample_column_dict: dict,
    logger=None,
    policies=None,
    ctx=None,
    step_number: int = 5,
    step_total: int = 8,
):
    """
    Infer the genome build by matching variant positions and alleles.

    Parameters
    ----------
    t_df : pl.DataFrame
        Subset of GWAS summary statistics.
    reference_files : dict
        Configured genome-build name to reference-table path.
    reference_column_mapping : dict
        Configured chromosome, position, allele and delimiter fields.
    sample_column_dict : dict
        Contains mapping for 'chr_col', 'pos_col', 'ea_col', 'oa_col'.
    logger : PipelineLogger, optional
        When given, and ``ctx`` is not, this function opens its own step on it.
    policies : Policies, optional
        When None the registry defaults are used, so behaviour is unchanged.
    ctx : StepContext, optional
        Pass this when the caller opened the step itself.

    Returns
    -------
    dict
        genome_build_info with counts, percentages, and inferred build.
        'total_variants' is the denominator behind the percentages: the number
        of study rows whose exact chromosome and position occurs in at least
        one configured reference.
        'input_variants' is the cleaned study-row count supplied to this step.
        'reference_percentages' measures unique matched reference allele records
        against the unique records on chromosomes present in the study.
    """
    policies = resolve_policies(policies)
    mode = policies.get("build.mode")
    confidence_ratio = float(policies.get("build.confidence_ratio"))
    min_match_count = int(policies.get("build.min_match_count"))
    min_reference_match_fraction = float(
        policies.get("build.min_reference_match_fraction")
    )
    min_match_fraction = float(policies.get("build.min_match_fraction"))
    deduplicate = bool(policies.get("build.deduplicate_reference"))
    build_names = list(reference_files)
    fallback_enabled = policies.get("build.low_coverage_fallback_enabled")
    if mode == "auto" and fallback_enabled and len(build_names) != 2:
        raise GenomeBuildError(
            "Automatic fallback requires exactly two configured references when "
            "'build.low_coverage_fallback_enabled' is true; supply two references "
            "or disable that policy."
        )
    if len(build_names) < 2:
        raise GenomeBuildError(
            "Genome-build inference requires at least two configured references."
        )

    with step_context(
        logger, ctx,
        number=step_number, total=step_total,
        title="Genome build", operation="genome_build.infer_genome_build",
        rows_in=t_df.height, policy_keys=POLICY_KEYS,
    ) as step:
        chr_col = sample_column_dict["chr_col"]
        pos_col = sample_column_dict["pos_col"]
        ea_col = sample_column_dict["ea_col"]
        oa_col = sample_column_dict["oa_col"]

        # --- Prepare GWAS dataframe ---------------------------------------
        if not has_canonical_variant_columns(t_df, sample_column_dict):
            # Keep direct library calls safe; the normal pipeline reaches this
            # step with dataset-level canonical columns and avoids this scan.
            t_df, _normalization = canonicalize_variant_frame(
                t_df,
                {"chr": chr_col, "pos": pos_col, "ea": ea_col, "oa": oa_col},
                policies,
                error_type=GenomeBuildError,
                warn=step.warn,
                label="study",
            )

        references = {
            name: _load_reference(
                path,
                name,
                step,
                deduplicate,
                reference_column_mapping,
                policies,
            )
            for name, path in reference_files.items()
            if mode == "auto" or name == mode
        }

        # One coordinate join supplies both build evidence and strand evidence.
        # Palindromic SNPs are excluded from strand counts because forward and
        # reverse-complement representations are indistinguishable.
        study_row = "__build_study_row"
        study = t_df.with_row_index(study_row)
        study_chromosomes = (
            study.select(pl.col(chr_col).alias("__build_chr"))
            .drop_nulls()
            .unique(maintain_order=True)
        )
        reference_marker_columns = [
            "__build_marker_chr",
            "__build_marker_pos",
            "__build_ref",
            "__build_alt",
        ]

        def match_evidence(ref_df: pl.DataFrame):
            relevant_reference = ref_df.join(
                study_chromosomes,
                on="__build_chr",
                how="semi",
            ).with_columns([
                pl.col("__build_chr").alias("__build_marker_chr"),
                pl.col("__build_pos").alias("__build_marker_pos"),
            ])
            reference_markers = relevant_reference.select(
                reference_marker_columns
            ).unique()
            joined = study.join(
                relevant_reference,
                left_on=[chr_col, pos_col],
                right_on=["__build_chr", "__build_pos"],
                how="inner",
            ).with_columns([
                reverse_complement_expression(
                    pl.col(ea_col)
                ).alias("__build_ea_rc"),
                reverse_complement_expression(
                    pl.col(oa_col)
                ).alias("__build_oa_rc"),
            ])
            ea = pl.col(ea_col)
            oa = pl.col(oa_col)
            ref = pl.col("__build_ref")
            alt = pl.col("__build_alt")
            forward = ((ea == alt) & (oa == ref)) | ((ea == ref) & (oa == alt))
            reverse = (
                ((pl.col("__build_ea_rc") == alt) & (pl.col("__build_oa_rc") == ref))
                | ((pl.col("__build_ea_rc") == ref) & (pl.col("__build_oa_rc") == alt))
            )
            snv = (
                (ea.str.len_chars() == 1)
                & (oa.str.len_chars() == 1)
                & (ref.str.len_chars() == 1)
                & (alt.str.len_chars() == 1)
            )
            non_palindromic = snv & ~palindromic_snp_expression(ea, oa)
            classified = joined.with_columns([
                forward.alias("__build_forward"),
                reverse.alias("__build_reverse"),
                non_palindromic.alias("__build_informative"),
            ])
            coordinate_rows = classified.select(study_row, chr_col).unique()
            matched_rows = classified.filter(
                pl.col("__build_forward") | pl.col("__build_reverse")
            )
            matched = matched_rows.select(study_row).unique().height
            matched_reference_markers = matched_rows.select(
                reference_marker_columns
            ).unique().height
            informative = classified.filter(pl.col("__build_informative"))
            forward_rows = informative.filter(
                pl.col("__build_forward") & ~pl.col("__build_reverse")
            ).select(study_row).unique().height
            reverse_rows = informative.filter(
                pl.col("__build_reverse") & ~pl.col("__build_forward")
            ).select(study_row).unique().height
            ambiguous_rows = informative.filter(
                pl.col("__build_forward") & pl.col("__build_reverse")
            ).select(study_row).unique().height
            # DNA SNP substitutions are an algorithmic evidence invariant:
            # indels, missing/ambiguous bases and identical alleles cannot
            # supply the fallback's distinct SNP-locus floor.
            snp_loci = matched_rows.filter(
                # Genomic positions are one-based; direct library calls may
                # bypass upstream coordinate QC, so invalid loci cannot count.
                (pl.col(pos_col) > 0)
                & ea.str.contains(r"^[ACGT]$")
                & oa.str.contains(r"^[ACGT]$")
                & ref.str.contains(r"^[ACGT]$")
                & alt.str.contains(r"^[ACGT]$")
                & (ea != oa) & (ref != alt)
            ).select(chr_col, pos_col).unique()
            return {
                "summary": {
                    "matches": matched,
                    "coordinate_hits": coordinate_rows.height,
                    "reference_markers": reference_markers.height,
                    "matched_reference_markers": matched_reference_markers,
                    "forward": forward_rows,
                    "reverse": reverse_rows,
                    "ambiguous": ambiguous_rows,
                },
                "coordinate_rows": coordinate_rows,
                "matched_rows": matched_rows.select(study_row, chr_col).unique(),
                "matched_loci": matched_rows.select(chr_col, pos_col).unique(),
                "snp_loci": snp_loci,
                "reference_chromosomes": reference_markers.get_column(
                    "__build_marker_chr"
                ).unique().to_list(),
            }

        # --- Compute match statistics -------------------------------------
        evidence_by_build = {
            name: match_evidence(reference)
            for name, reference in references.items()
        }
        orientation_evidence = {
            name: evidence["summary"] for name, evidence in evidence_by_build.items()
        }
        coordinate_rows = [
            evidence["coordinate_rows"] for evidence in evidence_by_build.values()
        ]
        testable_rows = (
            pl.concat(coordinate_rows, how="vertical").unique()
            if coordinate_rows else study.select(study_row, chr_col).head(0)
        )
        testable = testable_rows.height
        untestable = t_df.height - testable
        step.info(
            "{:,} of {:,} variants have an exact coordinate in at least one "
            "configured build reference; {:,} are not build-testable.".format(
                testable, t_df.height, untestable
            )
        )
        matches = {
            name: orientation_evidence.get(name, {}).get("matches", 0)
            for name in build_names
        }
        coordinate_hits = {
            name: orientation_evidence.get(name, {}).get("coordinate_hits", 0)
            for name in build_names
        }
        reference_marker_counts = {
            name: orientation_evidence.get(name, {}).get("reference_markers", 0)
            for name in build_names
        }
        matched_reference_marker_counts = {
            name: orientation_evidence.get(name, {}).get(
                "matched_reference_markers", 0
            )
            for name in build_names
        }
        fractions = {
            name: (count / testable if testable else 0.0)
            for name, count in matches.items()
        }
        input_fractions = {
            name: (count / t_df.height if t_df.height else 0.0)
            for name, count in matches.items()
        }
        reference_fractions = {
            name: (
                matched_reference_marker_counts[name] / count
                if count else 0.0
            )
            for name, count in reference_marker_counts.items()
        }

        # A locus is exclusive only if NO allele-compatible study evidence at
        # that canonical coordinate supports another reference, even when it
        # comes from a different study row/allele pair at the same position.
        exclusive_loci = {}
        for name, item in evidence_by_build.items():
            other_loci = [other["matched_loci"] for other_name, other in evidence_by_build.items()
                          if other_name != name]
            exclusive_loci[name] = (
                item["snp_loci"].join(pl.concat(other_loci).unique(),
                             on=[chr_col, pos_col], how="anti")
                if other_loci else item["snp_loci"].head(0)
            )
        distinct_loci = {
            "assessed_builds": list(evidence_by_build),
            "compatible_snp_loci": {
                name: item["snp_loci"].height for name, item in evidence_by_build.items()
            },
            "exclusive_snp_loci": {
                name: item.height for name, item in exclusive_loci.items()
            } if mode == "auto" else {},
        }
        # Aggregate compact evidence frames rather than scanning the study or
        # reference joins separately for every chromosome.
        def chromosome_counts(frame):
            return dict(frame.group_by(chr_col).len().iter_rows())

        testable_by_chromosome = chromosome_counts(testable_rows)
        matches_by_chromosome = {
            name: chromosome_counts(item["matched_rows"]) for name, item in evidence_by_build.items()
        }
        loci_by_chromosome = {
            name: chromosome_counts(item) for name, item in exclusive_loci.items()
        }
        per_chromosome = {}
        for chromosome in study_chromosomes.get_column("__build_chr").to_list():
            denominator = testable_by_chromosome.get(chromosome, 0)
            counts = {name: item.get(chromosome, 0)
                      for name, item in matches_by_chromosome.items()}
            per_chromosome[chromosome] = {
                "testable_variants": denominator,
                "matches": counts,
                "match_fractions": {
                    name: count / denominator if denominator else 0.0
                    for name, count in counts.items()
                },
                "exclusive_snp_loci": {
                    name: item.get(chromosome, 0) for name, item in loci_by_chromosome.items()
                } if mode == "auto" else {},
                "competitor_veto": False,
                "assessment": "not_assessed" if mode != "auto" else "unassessed",
            }

        total_matches = sum(matches.values())
        match_fraction = max(fractions.values())
        input_match_fraction = max(input_fractions.values())
        reference_match_fraction = max(reference_fractions.values())
        ambiguous_reason = None
        confidence = None
        winner = None
        normal_failed_gates = []
        if mode != "auto":
            inferred_build = mode
            confidence = None
            match_fraction = fractions.get(mode, 0.0)
            input_match_fraction = input_fractions.get(mode, 0.0)
            reference_match_fraction = reference_fractions.get(mode, 0.0)
        elif testable == 0:
            inferred_build = "Ambiguous"
            ambiguous_reason = (
                "no study coordinate occurs in any configured build reference"
            )
        elif total_matches == 0:
            inferred_build = "Ambiguous"
            ambiguous_reason = (
                "none of the {:,} coordinate-testable variants had compatible "
                "alleles in any configured reference".format(testable)
            )
        else:
            ordered = sorted(matches.items(), key=lambda item: item[1], reverse=True)
            winner, winning_matches = ordered[0]
            tied = len(ordered) > 1 and winning_matches == ordered[1][1]
            confidence = winning_matches / total_matches
            reference_match_fraction = reference_fractions[winner]
            normal_failed_gates = [key for key, passed in (
                ("unique_winner", not tied),
                ("build.min_match_count", winning_matches >= min_match_count),
                ("build.min_reference_match_fraction", reference_match_fraction >= min_reference_match_fraction),
                ("build.min_match_fraction", match_fraction >= min_match_fraction),
                ("build.confidence_ratio", confidence >= confidence_ratio),
            ) if not passed]
            if tied:
                inferred_build = "Ambiguous"
                ambiguous_reason = "the leading genome-build references tied"
            elif winning_matches < min_match_count:
                inferred_build = "Ambiguous"
                ambiguous_reason = (
                    "%s matched only %s allele-compatible variants, below the "
                    "%s required by 'build.min_match_count'"
                    % (
                        winner,
                        "{:,}".format(winning_matches),
                        "{:,}".format(min_match_count),
                    )
                )
            elif reference_match_fraction < min_reference_match_fraction:
                inferred_build = "Ambiguous"
                ambiguous_reason = (
                    "only %s of %s unique %s reference markers on the study's "
                    "chromosomes matched (%.4f%%), below the %.2f%% required "
                    "by 'build.min_reference_match_fraction'"
                    % (
                        "{:,}".format(matched_reference_marker_counts[winner]),
                        "{:,}".format(reference_marker_counts[winner]),
                        winner,
                        reference_match_fraction * 100,
                        min_reference_match_fraction * 100,
                    )
                )
            elif match_fraction < min_match_fraction:
                inferred_build = "Ambiguous"
                ambiguous_reason = (
                    "only %.4f%% of the %s coordinate-testable variants "
                    "matched %s, below the %.2f%% required by "
                    "'build.min_match_fraction'"
                    % (
                        match_fraction * 100,
                        "{:,}".format(testable),
                        winner,
                        min_match_fraction * 100,
                    )
                )
            elif confidence >= confidence_ratio:
                inferred_build = winner
            else:
                inferred_build = "Ambiguous"
                ambiguous_reason = (
                    "no build reached the %.2f share of allele matches required "
                    "by 'build.confidence_ratio' (best: %s %.3f)"
                    % (confidence_ratio, winner, confidence)
                )

        fallback = {
            "enabled": fallback_enabled, "attempted": False, "accepted": False,
            "status": "disabled" if not fallback_enabled else "not_applicable",
            "winner": winner, "competitor": None, "gates": {}, "failed_gates": [],
            "veto_chromosomes": [],
            "missing_reference_chromosomes": {},
            "thresholds": {
                key: policies.get("build." + key) for key in (
                    "low_coverage_min_match_fraction", "low_coverage_max_competitor_fraction",
                    "low_coverage_min_unique_loci", "low_coverage_chromosome_min_competitor_loci",
                    "low_coverage_chromosome_min_competitor_fraction",
                )
            },
        }
        acceptance_path = "declared" if mode != "auto" else (
            "normal" if inferred_build != "Ambiguous" else "ambiguous"
        )
        if mode == "auto" and winner is not None and len(build_names) == 2:
            competitor = next(name for name in build_names if name != winner)
            fallback["competitor"] = competitor
            # A zero competing count is interpretable only where both panels
            # represent the chromosome. This comparability invariant applies
            # to fallback alone; it does not alter normal/declared decisions.
            fallback["missing_reference_chromosomes"] = {
                name: sorted(set(testable_by_chromosome) - set(item["reference_chromosomes"]))
                for name, item in evidence_by_build.items()
                if set(testable_by_chromosome) - set(item["reference_chromosomes"])
            }
            for chromosome, item in per_chromosome.items():
                competitor_loci = item["exclusive_snp_loci"][competitor]
                veto = (
                    competitor_loci >= policies.get("build.low_coverage_chromosome_min_competitor_loci")
                    and item["match_fractions"][competitor] >= policies.get("build.low_coverage_chromosome_min_competitor_fraction")
                )
                item["competitor_veto"] = veto
                item["assessment"] = (
                    "competitor_dominated" if veto else "sparse" if
                    max(item["exclusive_snp_loci"].values()) < policies.get("build.low_coverage_chromosome_min_competitor_loci")
                    else "no_veto_detected"
                ) if item["testable_variants"] else "unassessed"
                if veto:
                    fallback["veto_chromosomes"].append(chromosome)
            fallback["gates"] = {
                "build.low_coverage_min_match_fraction": fractions[winner] >= policies.get("build.low_coverage_min_match_fraction"),
                "build.low_coverage_max_competitor_fraction": fractions[competitor] <= policies.get("build.low_coverage_max_competitor_fraction"),
                "build.low_coverage_min_unique_loci": distinct_loci["exclusive_snp_loci"][winner] >= policies.get("build.low_coverage_min_unique_loci"),
                "chromosome_veto": not fallback["veto_chromosomes"],
                "comparable_reference_scope": not fallback["missing_reference_chromosomes"],
            }
            fallback["failed_gates"] = [key for key, passed in fallback["gates"].items() if not passed]
        if fallback_enabled and mode == "auto":
            if acceptance_path == "normal":
                fallback["status"] = "not_needed"
            elif normal_failed_gates == ["build.min_reference_match_fraction"]:
                fallback["attempted"] = True
                fallback["accepted"] = not fallback["failed_gates"]
                fallback["status"] = "accepted" if fallback["accepted"] else "rejected"
                if fallback["accepted"]:
                    inferred_build = winner
                    ambiguous_reason = None
                    acceptance_path = "low_coverage_fallback"
                    step.warn(
                        "Experimental dominant-build inference accepted via low_coverage_fallback; "
                        "reference coverage alone was waived. Match fractions are not a posterior "
                        "probability and small mixed components can pass."
                    )
                else:
                    ambiguous_reason += "; low-coverage fallback failed: " + ", ".join(fallback["failed_gates"])
                    if fallback["veto_chromosomes"]:
                        ambiguous_reason += " (competitor-dominated chromosomes: %s)" % ", ".join(fallback["veto_chromosomes"])
                    if fallback["missing_reference_chromosomes"]:
                        ambiguous_reason += (
                            " (references missing testable chromosomes: %s; supply comparable "
                            "panels covering these chromosomes)" % fallback["missing_reference_chromosomes"]
                        )
            else:
                ambiguous_reason += (
                    "; low-coverage fallback is inapplicable: reference coverage must be the sole failed "
                    "normal gate; failed gates: %s" % ", ".join(normal_failed_gates or ["no_compatible_evidence"])
                )

        evidence = {
            "input variants": "{:,}".format(t_df.height),
            "coordinate-testable variants": "{:,}".format(testable),
            "not build-testable": "{:,}".format(untestable),
        }
        for name in build_names:
            evidence["%s coordinate hits" % name] = "{:,}".format(
                coordinate_hits[name]
            )
            evidence["%s allele matches" % name] = "{:,} ({:.2f}%)".format(
                matches[name], fractions[name] * 100,
            )
            evidence["%s input-wide match fraction" % name] = "{:.2f}%".format(
                input_fractions[name] * 100,
            )
            evidence["%s relevant reference markers" % name] = "{:,}".format(
                reference_marker_counts[name]
            )
            evidence["%s matched reference markers" % name] = (
                "{:,} ({:.2f}%)".format(
                    matched_reference_marker_counts[name],
                    reference_fractions[name] * 100,
                )
            )
        evidence["minimum allele matches"] = "{:,}".format(min_match_count)
        evidence["minimum relevant-reference coverage"] = "{:.2f}%".format(
            min_reference_match_fraction * 100
        )
        evidence["minimum match fraction"] = "{:.2f}%".format(
            min_match_fraction * 100
        )
        evidence["winning share of matches"] = (
            "n/a" if confidence is None else "{:.3f}".format(confidence)
        )
        evidence["acceptance path"] = acceptance_path
        evidence["low_coverage_fallback"] = fallback["status"]
        for name, count in distinct_loci["exclusive_snp_loci"].items():
            evidence["%s exclusive SNP loci" % name] = "{:,}".format(count)
        step.observed("Genome-build distinct loci", **distinct_loci)
        step.observed("Genome-build low-coverage fallback", **fallback)
        decision = inferred_build
        if ambiguous_reason:
            decision = "%s - %s" % (inferred_build, ambiguous_reason)
        step.decide("Genome build", evidence, decision)
        if inferred_build == "Ambiguous":
            step.warn(
                "The genome build could not be inferred. Set policy "
                "'build.mode' to one of the configured builds to state it outright."
            )

        # --- Summary output -----------------------------------------------
        genome_build_info = {
            "inferred_build": inferred_build,
            "matches": matches,
            "percentages": {
                name: round(fraction * 100, 2)
                for name, fraction in fractions.items()
            },
            "input_percentages": {
                name: round(fraction * 100, 2)
                for name, fraction in input_fractions.items()
            },
            "reference_marker_counts": reference_marker_counts,
            "matched_reference_marker_counts": matched_reference_marker_counts,
            "reference_percentages": {
                name: round(fraction * 100, 2)
                for name, fraction in reference_fractions.items()
            },
            "total_variants": testable,
            "input_variants": t_df.height,
            "testable_variants": testable,
            "untestable_variants": untestable,
            "coordinate_hits": coordinate_hits,
            "match_fraction": match_fraction,
            "input_match_fraction": input_match_fraction,
            "reference_match_fraction": reference_match_fraction,
            "confidence": confidence,
            "mode": mode,
            "forced": mode != "auto",
            "ambiguous_reason": ambiguous_reason,
            "strand_evidence": orientation_evidence,
            "acceptance_path": acceptance_path,
            "normal_failed_gates": normal_failed_gates,
            "distinct_loci": distinct_loci,
            "per_chromosome": per_chromosome,
            "low_coverage_fallback": fallback,
        }
        step.extra["genome_build"] = genome_build_info
        step.set_rows(t_df.height, removed=0)
    return genome_build_info
