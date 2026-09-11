# Manhattan plots

`postgwas pipeline --modules manhattan` and `postgwas manhattan` consume an
indexed, single-sample GWAS-VCF. The declared genome build is validated from
the header, not guessed from a default. Direct `--genome-build` is an optional
assertion that must agree with that declaration. `--pheno`, when supplied,
must match the validated single sample.

Settings resolve as packaged YAML, then run YAML, then explicit CLI options.
`modules.manhattan.file_format` selects PDF or PNG unless `--pdf` or `--png`
explicitly selects an output; those two flags are mutually exclusive.
`output_file`, `log_file`, and `plot_data_file` are relative patterns rendered
inside the output directory. `single_chromosome_height`, `genomewide_height`,
and `png_dpi` configure the automatic figure dimensions and raster resolution.
The original explicit width, height, spacing, font and threshold flags remain
available. `--max-height` is a ceiling in untransformed negative-log10-P units.
`axis_label_rows` staggers crowded chromosome tick labels without moving
points. `caption_template` and `caption_font_size` control the caption, which
reports the plotted subset, the display transformation and empty chromosomes.

The default fields are FORMAT/LP and FORMAT/AF. Configurable `pvalue_field`
and `frequency_field` must declare one numeric value per association (Number=1
or Number=A for biallelic records). Multiallelic plotted records fail rather
than silently selecting one alternate allele.
Optional allelic-shift mode instead requires the two numeric counts in the
configured INFO `allelic_shift_field`. It uses a symmetric exact two-sided
binomial test, capped at probability one. Allelic-shift mode cannot be combined
with the allele-frequency filter. Coding highlighting uses `consequence_info_field`,
`consequence_field`, and `coding_terms`, and requires a functional bcftools
split-vep plugin. Missing consequences are retained as unhighlighted records.

Preflight loads the actual R package stack before analysis. The native log
records R, optparse, data.table and ggplot2 versions. The pilot environment
required an explicit `R_LIBS_USER` pointing to the conda R library because a
user-level R library shadowed it with incompatible binaries; this selection
is retained in the recorded pilot commands and does not alter user settings.

The headerless bcftools query is read with `header=FALSE`, preserving its first
record. Chromosome offsets are named by assembly chromosome, so chromosomes
absent from the input cannot shift later points into the wrong chromosome.
The bundled GRCh37/GRCh38 chromosome lengths are assembly reference invariants;
unsupported chromosome labels or out-of-assembly coordinates fail explicitly.
Chromosomes 1–22 and X are supported. The configured `chromosome_label_policy`
defaults to `strip_prefix` with `chromosome_prefix: chr`, preserving accepted
chr-prefixed labels without modifying the VCF. The audit retains `source_chrom`
and the log counts changed labels. Set the policy to `exact` to disallow this
normalization. M/MT are unsupported by the bundled chromosome lengths and fail
explicitly instead of creating undefined coordinates.

For LP values above threshold t, the display is t log10(LP) / log10(t).
Values at or below t are unchanged; t must exceed one to maintain a monotone,
continuous scale. The same transformation applies to the configured
significance and suggestive reference lines and plot ceiling. The point audit
retains the original LP, plotted LP, chromosome, original position and cumulative
position. The log records the resolved configuration, command, input count,
retained count and exclusion count. Empty results or a failed query cannot be
reported as a completed plot. This module produces Manhattan plots, not QQ plots.

Method and interface references:

- [bcftools query](https://samtools.github.io/bcftools/bcftools.html#query)
- [data.table fread header interpretation](https://rdatatable.gitlab.io/data.table/reference/fread.html)
- [bcftools split-vep and keep-sites](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html)
- [Original MIT-licensed plotting adapter](https://github.com/freeseek/score/blob/master/assoc_plot.R)
- [R binomial distribution functions](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html)
- [UCSC hg19 chromosome lengths](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes)
- [UCSC hg38 chromosome lengths](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)

The upstream GRCh37 adapter overstated lengths of chromosomes 1, 20 and 22
by 1,000 bp each. They are corrected to 249250621, 63025520 and 51304566 bp,
respectively, matching UCSC and the pilot VCF contig declarations. This changes
cumulative display coordinates, not the source positions or association values.

The native regression tests cover first-record retention, nonconsecutive
chromosomes, fractional plotting settings, reference thresholds, YAML/CLI
precedence, exact binomial boundary cases, missing LP, empty plots, and invalid
inputs. Pilot validation does not establish the visual appearance of every
possible configuration or annotation dataset.
