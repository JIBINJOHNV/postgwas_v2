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
            coordinate_rows = classified.select(study_row).unique()
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
            return (
                {
                    "matches": matched,
                    "coordinate_hits": coordinate_rows.height,
                    "reference_markers": reference_markers.height,
                    "matched_reference_markers": matched_reference_markers,
                    "forward": forward_rows,
                    "reverse": reverse_rows,
                    "ambiguous": ambiguous_rows,
                },
                coordinate_rows,
            )

        # --- Compute match statistics -------------------------------------
        evidence_by_build = {
            name: match_evidence(reference)
            for name, reference in references.items()
        }
        orientation_evidence = {
            name: evidence[0] for name, evidence in evidence_by_build.items()
        }
        coordinate_rows = [
            evidence[1] for evidence in evidence_by_build.values()
        ]
        testable = (
            pl.concat(coordinate_rows, how="vertical").unique().height
            if coordinate_rows else 0
        )
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

        total_matches = sum(matches.values())
        match_fraction = max(fractions.values())
        input_match_fraction = max(input_fractions.values())
        reference_match_fraction = max(reference_fractions.values())
        ambiguous_reason = None
        confidence = None
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
        }
        step.set_rows(t_df.height, removed=0)
    return genome_build_info
