"""Dataset step 05 - infer the genome build of the input.

Implements the build rows of plan v3 section 6.4: an explicit ``build.mode``
instead of the ``--genome-build`` flag the current error message advertises but
that does not exist, absolute and relative evidence floors
(``build.min_match_count`` and ``build.min_match_fraction``), optional
deduplication of the references, and a denominator that counts only variants
whose exact coordinate occurs in at least one configured build reference.

Build evidence follows the GWAS Catalog summary-statistics convention that a
variant occupies one row while either allele order is valid. Consequently, a
study row contributes at most one match to each candidate build, even if the
reference repeats a variant or contains both allele orientations. Sources:
https://www.ebi.ac.uk/gwas/docs/summary-statistics-format and
https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics.
"""

import polars as pl

from postgwas.core.dataframes import (
    chromosome_expression,
    count_non_null,
    position_expression,
    validate_cast_retention,
)
from postgwas.core.io.tables import read_delimited_table

from .shared.runtime import resolve_policies, step_context
from .shared.variant_columns import has_canonical_variant_columns

__all__ = ["infer_genome_build", "GenomeBuildError", "POLICY_KEYS"]


POLICY_KEYS = [
    "build.mode",
    "build.confidence_ratio",
    "build.min_match_count",
    "build.min_match_fraction",
    "build.deduplicate_reference",
    "chromosome.strip_chr_prefix",
]


def _reverse_complement(column: str) -> pl.Expr:
    """Reverse-complement a sequence allele using IUPAC DNA bases A/C/G/T."""
    return (
        pl.col(column)
        .cast(pl.Utf8)
        .str.to_uppercase()
        .str.replace_many(["A", "C", "G", "T"], ["T", "G", "C", "A"])
        .str.reverse()
    )


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
    strip_chr_prefix,
    column_mapping,
    policies,
):
    """Read one build reference and normalise its coordinates."""
    step.input(
        "Genome-build reference",
        build=label,
        file=str(path),
        delimiter=column_mapping["delimiter"],
        chromosome_column=column_mapping["chr"],
        position_column=column_mapping["pos"],
        reference_allele_column=column_mapping["a1"],
        alternate_allele_column=column_mapping["a2"],
    )
    raw, _detected = read_delimited_table(
        path,
        column_mapping["delimiter"],
        candidates=list(policies.get("input.delimiter_candidates")),
        minimum_columns=int(policies.get("input.delimiter_min_columns")),
        maximum_columns=int(policies.get("input.delimiter_max_columns")),
        sample_lines=int(policies.get("input.delimiter_sample_rows")),
        null_values=list(policies.get("input.null_values")),
        infer_schema_length=int(policies.get("input.schema_inference_rows")),
        error_type=GenomeBuildError,
        description="%s genome-build reference" % label,
    )
    configured_columns = [
        column_mapping["chr"],
        column_mapping["pos"],
        column_mapping["a1"],
        column_mapping["a2"],
    ]
    missing = [column for column in configured_columns if column not in raw.columns]
    if missing:
        raise GenomeBuildError(
            "%s genome-build reference is missing configured columns: %s"
            % (label, ", ".join(missing))
        )
    raw = raw.select(configured_columns).rename(dict(zip(
        configured_columns,
        ["__build_chr", "__build_pos", "__build_ref", "__build_alt"],
    )))
    schema = dict(raw.schema)
    position_before = count_non_null(raw, "__build_pos")
    reference = raw.with_columns([
        chromosome_expression(
            "__build_chr",
            schema["__build_chr"],
            strip_chr_prefix=strip_chr_prefix,
            strip_leading_zero=True,
        ),
        position_expression("__build_pos", schema["__build_pos"]),
        pl.col("__build_ref").cast(pl.Utf8).str.to_uppercase(),
        pl.col("__build_alt").cast(pl.Utf8).str.to_uppercase(),
    ])
    validate_cast_retention(
        position_before,
        count_non_null(reference, "__build_pos"),
        "%s reference position" % label,
        error_type=GenomeBuildError,
        warn=step.warn,
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
        'input_variants' is the unfiltered row count.
    """
    policies = resolve_policies(policies)
    mode = policies.get("build.mode")
    confidence_ratio = float(policies.get("build.confidence_ratio"))
    min_match_count = int(policies.get("build.min_match_count"))
    min_match_fraction = float(policies.get("build.min_match_fraction"))
    deduplicate = bool(policies.get("build.deduplicate_reference"))
    strip_chr_prefix = bool(policies.get("chromosome.strip_chr_prefix"))
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
            schema = dict(t_df.schema)
            position_before = count_non_null(t_df, pos_col)
            t_df = t_df.with_columns([
                chromosome_expression(
                    chr_col, schema[chr_col], strip_chr_prefix=strip_chr_prefix,
                    strip_leading_zero=True,
                ),
                position_expression(pos_col, schema[pos_col]),
                pl.col(ea_col).cast(pl.Utf8).str.to_uppercase(),
                pl.col(oa_col).cast(pl.Utf8).str.to_uppercase(),
            ])
            validate_cast_retention(
                position_before,
                count_non_null(t_df, pos_col),
                "study position",
                error_type=GenomeBuildError,
                warn=step.warn,
            )

        references = {
            name: _load_reference(
                path,
                name,
                step,
                deduplicate,
                strip_chr_prefix,
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

        def match_evidence(ref_df: pl.DataFrame):
            joined = study.join(
                ref_df,
                left_on=[chr_col, pos_col],
                right_on=["__build_chr", "__build_pos"],
                how="inner",
            ).with_columns([
                _reverse_complement(ea_col).alias("__build_ea_rc"),
                _reverse_complement(oa_col).alias("__build_oa_rc"),
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
            pair = pl.concat_str([ea, oa])
            non_palindromic = snv & ~pair.is_in(["AT", "TA", "CG", "GC"])
            classified = joined.with_columns([
                forward.alias("__build_forward"),
                reverse.alias("__build_reverse"),
                non_palindromic.alias("__build_informative"),
            ])
            coordinate_rows = classified.select(study_row).unique()
            matched = classified.filter(
                pl.col("__build_forward") | pl.col("__build_reverse")
            ).select(study_row).unique().height
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
        fractions = {
            name: (count / testable if testable else 0.0)
            for name, count in matches.items()
        }

        total_matches = sum(matches.values())
        match_fraction = max(fractions.values())
        ambiguous_reason = None
        confidence = None
        if mode != "auto":
            inferred_build = mode
            confidence = None
            match_fraction = fractions.get(mode, 0.0)
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
        evidence["minimum allele matches"] = "{:,}".format(min_match_count)
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
            "total_variants": testable,
            "input_variants": t_df.height,
            "testable_variants": testable,
            "untestable_variants": untestable,
            "coordinate_hits": coordinate_hits,
            "match_fraction": match_fraction,
            "confidence": confidence,
            "mode": mode,
            "forced": mode != "auto",
            "ambiguous_reason": ambiguous_reason,
            "strand_evidence": orientation_evidence,
        }
        step.set_rows(t_df.height, removed=0)
    return genome_build_info
