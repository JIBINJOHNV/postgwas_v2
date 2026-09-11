#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2021-2025 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

options(error = function() {traceback(3); q("no", 1)})

assoc_plot_version <- '2025-08-19-postgwas-validated'
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
if (capabilities()[['cairo']]) options(bitmapType = 'cairo')

parser <- OptionParser('usage: assoc_plot.R [options] --genome <GRCh37|GRCh38>|--cytoband <cytoband.txt.gz> --vcf|--tbx <file>')
parser <- add_option(parser, c('--genome'), type = 'character', help = 'genome assembly (e.g. GRCh38)', metavar = '<assembly>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--nauto'), type = 'integer', help = 'number of autosomes', metavar = '<integer>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--pheno'), type = 'character', help = 'phenotype to select from GWAS-VCF file', metavar = '<string>')
parser <- add_option(parser, c('--as'), action = 'store_true', default = FALSE, help = 'input VCF file has allelic shift information')
parser <- add_option(parser, c('--csq'), action = 'store_true', default = FALSE, help = 'whether coding variant should be flagged as red')
parser <- add_option(parser, c('--tbx'), type = 'character', help = 'input REGENIE/PLINK summary statistics', metavar = '<file.gz>')
parser <- add_option(parser, c('--region'), type = 'character', help = 'region to plot', metavar = '<region>')
parser <- add_option(parser, c('--min-af'), type = 'double', help = 'minimum minor allele frequency [0]', metavar = '<float>')
parser <- add_option(parser, c('--min-lp'), type = 'double', help = 'minimum -log10 p-val [2]', metavar = '<integer>')
parser <- add_option(parser, c('--loglog-pval'), type = 'double', help = '-log10 p-val threshold for using log-log scale in manhattan plot', metavar = '<integer>')
parser <- add_option(parser, c('--cyto-ratio'), type = 'double', help = 'plot to cytoband ratio [25]', metavar = '<integer>')
parser <- add_option(parser, c('--max-height'), type = 'double', help = 'raw negative-log10-P ceiling (before display transformation)', metavar = '<float>')
parser <- add_option(parser, c('--spacing'), type = 'integer', help = 'spacing between chromosomes [10]', metavar = '<integer>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'double', help = 'inches width of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--height'), type = 'double', help = 'inches height of the output file [3.5/7.0]', metavar = '<float>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', help = 'font size [12]', metavar = '<integer>')

parser <- add_option(parser, c('--bcftools'), type = 'character')
parser <- add_option(parser, c('--plot-data'), type = 'character')
parser <- add_option(parser, c('--threads'), type = 'integer')
parser <- add_option(parser, c('--significance'), type = 'double')
parser <- add_option(parser, c('--suggestive'), type = 'double')
parser <- add_option(parser, c('--single-height'), type = 'double')
parser <- add_option(parser, c('--genome-height'), type = 'double')
parser <- add_option(parser, c('--dpi'), type = 'integer')
parser <- add_option(parser, c('--lp-field'), type = 'character')
parser <- add_option(parser, c('--af-field'), type = 'character')
parser <- add_option(parser, c('--as-field'), type = 'character')
parser <- add_option(parser, c('--consequence-field'), type = 'character')
parser <- add_option(parser, c('--consequence-info-field'), type = 'character')
parser <- add_option(parser, c('--coding-terms'), type = 'character')
parser <- add_option(parser, c('--chromosome-label-policy'), type = 'character')
parser <- add_option(parser, c('--chromosome-prefix'), type = 'character')
parser <- add_option(parser, c('--axis-label-rows'), type = 'integer')
parser <- add_option(parser, c('--caption-fontsize'), type = 'double')
parser <- add_option(parser, c('--caption-template'), type = 'character')
parser <- add_option(parser, c('--input-records'), type = 'integer')
parser <- add_option(parser, c('--check-runtime'), action = 'store_true', default = FALSE)

args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)
runtime_packages <- c('optparse', 'data.table', 'ggplot2')
runtime_versions <- vapply(runtime_packages, function(x) as.character(packageVersion(x)), character(1))
runtime_version <- paste(c(R.version.string, paste(runtime_packages, runtime_versions)), collapse = '; ')
if (args$check_runtime) {
  cat(runtime_version, '\n')
  quit(status = 0)
}
write(runtime_version, stderr())

# The Python boundary supplies schema-validated canonical values.
required_options <- c('nauto', 'min_af', 'min_lp', 'loglog_pval', 'cyto_ratio', 'spacing', 'width', 'fontsize', 'bcftools', 'plot_data', 'threads', 'significance', 'suggestive', 'single_height', 'genome_height', 'dpi', 'lp_field', 'af_field', 'as_field', 'consequence_field', 'consequence_info_field', 'coding_terms')
missing_options <- required_options[vapply(required_options, function(x) is.null(args[[x]]), logical(1))]
if (is.null(args$chromosome_label_policy) || is.null(args$chromosome_prefix)) stop('Missing resolved chromosome label policy or prefix')
if (any(vapply(c('axis_label_rows', 'caption_fontsize', 'caption_template', 'input_records'), function(x) is.null(args[[x]]), logical(1)))) stop('Missing resolved axis or caption settings')
if (length(missing_options)) stop(paste('Missing resolved adapter options:', paste(missing_options, collapse = ', ')))
if (args$loglog_pval <= 1) stop('loglog-pval must exceed one for a monotone log-log scale')
setDTthreads(args$threads)

write(paste('assoc_plot.R', assoc_plot_version, 'http://github.com/freeseek/score'), stderr())

if (is.null(args$genome) && is.null(args$cytoband)) {print_help(parser); stop('either --genome or --cytoband is required\nTo download the cytoband file run:\nwget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz')}
if (!is.null(args$genome) && !is.null(args$cytoband)) {print_help(parser); stop('cannot use --genome and --cytoband at the same time')}
if (!is.null(args$genome) && args$genome != 'GRCh37' && args$genome != 'GRCh38') {print_help(parser); stop('--genome accepts only GRCh37 or GRCh38')}
if (is.null(args$vcf) && is.null(args$tbx)) {print_help(parser); stop('either --vcf or --tbx is required')}
if (!is.null(args$vcf) && !is.null(args$tbx)) {print_help(parser); stop('cannot use --vcf and --tbx at the same time')}
if (is.null(args$vcf) && !is.null(args$pheno))  {print_help(parser); stop('--pheno requires --vcf')}
if (is.null(args$vcf) && args$as) {print_help(parser); stop('--as requires --vcf')}
if (args$as && args$min_af > 0) {print_help(parser); stop('--as cannot be used with --min-af')}
if (is.null(args$vcf) && is.null(args$as_vcf) && args$csq)  {print_help(parser); stop('--csq requires --vcf')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option')}

if (!is.null(args$cytoband)) {
  df_cyto <- setNames(read.table(args$cytoband, sep = '\t', header = FALSE), c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
  df_cyto$chrom <- gsub('chr', '', df_cyto$chrom)
  chrlen <- tapply(df_cyto$chromEnd, df_cyto$chrom, max)
  chrs <- unique(df_cyto$chrom)
  modified_chrs <- gsub('^M[T]?$', args$nauto + 4, gsub('^Y$', args$nauto + 2, gsub('^X$', args$nauto + 1, chrs)))
  ord <- order(suppressWarnings(as.numeric(modified_chrs)))
  chrs <- chrs[ord]

  idx_p <- df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11'
  idx_q <- df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11'
  if (sum(idx_p) > 0 && sum(idx_q)) {
    df_cen <- rbind(cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1),
                    cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1/2),
                    cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = 0),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1/2),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = 0))
  }

  chrlen <- chrlen[c(1:args$nauto, 'X')]
} else if ( args$genome == 'GRCh37' ) {
  # Assembly invariants: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
  chrlen <- setNames(c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560), c(1:22,'X'))
} else if ( args$genome == 'GRCh38' ) {
  chrlen <- setNames(c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895), c(1:22,'X'))
}


# ---------------------------------
# Inflate widths of small chromosomes , it will increase the chromosome width 
# ---------------------------------

if ( !is.null(args$vcf) ) {
  fmt <- '%CHROM\\t%POS\\t%ALT'
  names <- c('chrom', 'pos', 'alt')
  if (args$as) {
    fmt <- paste0(fmt, '\\t%', args$as_field, '{0}\\t%', args$as_field, '{1}')
    names <- c(names, 'as0', 'as1')
  } else {
    fmt <- paste0(fmt, '[\\t%', args$lp_field, '{0}]')
    names <- c(names, 'lp')
  }
  if (!args$as && (args$min_af>0 || args$min_lp>0)) {
    opt_filter <- ' --include \''
    if (args$min_af>0) opt_filter <- paste0(opt_filter, 'FORMAT/', args$af_field, '>', args$min_af, ' & FORMAT/', args$af_field, '<1-', args$min_af)
    if (args$min_af>0 && args$min_lp>0) opt_filter <- paste0(opt_filter, ' & ')
    if (args$min_lp>0) opt_filter <- paste0(opt_filter, 'FORMAT/', args$lp_field, '>=', args$min_lp)
    opt_filter <- paste0(opt_filter, '\'')
  } else if (args$as) {
    opt_filter <- paste0(" --include ", shQuote(paste0("sum(", args$as_field, ")>0")))
  } else opt_filter <- ''
  if (!is.null(args$region)) opt_regions <- paste0(' --regions ', shQuote(args$region)) else opt_regions <- ''
  if (!is.null(args$pheno)) opt_samples <- paste0(' --samples ', shQuote(args$pheno)) else opt_samples <- ''
  if (args$csq) {
    fmt <- paste0(fmt, '\\t%', args$consequence_field)
    names <- c(names, 'consequence')
    # --keep-sites preserves records without annotations; only colouring changes.
    cmd <- paste0(shQuote(args$bcftools), ' +split-vep --annotation ', shQuote(args$consequence_info_field), ' --keep-sites --format ', shQuote(paste0(fmt, '\\n')), opt_filter, opt_regions, ' ', shQuote(args$vcf))
  } else {
    cmd <- paste0(shQuote(args$bcftools), ' query --format ', shQuote(paste0(fmt, '\\n')), opt_filter, opt_regions, opt_samples, ' ', shQuote(args$vcf))
  }
} else {
  if ( !is.null(args$region) ) {
    cmd <- past0e('tabix --print-header "', args$tbx, '" ', strsplit(args$region,','))
  } else {
    cmd <- paste0('zcat "', args$tbx, '"')
  }
  if (!is.null(args$min_af)) {
    filter <- paste0('$af>', args$min_af, ' && $af<', 1-args$min_af, ' && $lp!="NA" && $lp>=', args$min_lp)
  } else {
    filter <- paste0('$lp!="NA" && $lp>=', args$min_lp)
  }
  cmd <- paste(cmd, '| awk \'NR==1 {for (i=1; i<=NF; i++) f[$i] = i; if ("CHROM" in f) chrom=f["CHROM"]; else chrom=f["#CHROM"]; if ("GENPOS" in f) pos=f["GENPOS"]; else pos=f["POS"]; if ("A1FREQ" in f) af=f["A1FREQ"]; else af=f["A1_FREQ"]; if ("NEG_LOG10_P" in f) lp=f["NEG_LOG10_P"]; else if ("LOG10P" in f) lp=f["LOG10P"]; else lp=f["LOG10_P"]} NR==1 || NR>1 &&', filter , '{print $chrom"\\t"$pos"\\t"$lp}\'')
  names <- c('chrom', 'pos', 'lp')
}

write(paste('Command:', cmd), stderr())
query_file <- tempfile()
query_status <- system(paste(cmd, '>', shQuote(query_file)))
if (query_status != 0) {
  unlink(query_file)
  stop(paste('Association query failed with exit status', query_status))
}
df <- tryCatch(
  setNames(fread(file = query_file, sep = '\t', header = is.null(args$vcf), na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names),
  finally = unlink(query_file)
)

input_rows <- nrow(df)
if (!is.null(args$vcf) && any(grepl(',', df$alt, fixed = TRUE))) stop('Manhattan plots require biallelic associations; split multiallelic records before plotting')
if (args$as) {
  counts <- c(df$as0, df$as1)
  if (any(!is.finite(counts) | counts < 0 | counts != floor(counts))) stop('Allelic-shift counts must be finite nonnegative integers')
  # A symmetric exact two-sided binomial probability cannot exceed one.
  df$lp = -pmin(0, log(2) + pbinom(pmin(df$as0, df$as1), df$as0 + df$as1, .5, log.p = TRUE)) / log(10)
  df <- df[df$lp >= args$min_lp & !is.na(df$lp),]
} else {
  df <- df[!is.na(df$lp),]
}

if (nrow(df) == 0) {
  stop('Nothing to be plotted after the configured filters')
}

write(paste('Query rows:', input_rows, '; retained:', nrow(df), '; missing/filtered LP:', input_rows - nrow(df)), stderr())
if (any(!is.finite(df$lp) | df$lp < 0)) stop('Association LP values must be finite and nonnegative')
df$source_chrom <- as.character(df$chrom)
df$chrom <- df$source_chrom
if (args$chromosome_label_policy == 'strip_prefix') {
  selected <- startsWith(df$chrom, args$chromosome_prefix)
  df$chrom[selected] <- substring(df$chrom[selected], nchar(args$chromosome_prefix) + 1)
} else if (args$chromosome_label_policy != 'exact') stop('Unsupported chromosome label policy')
write(paste('Chromosome label policy:', args$chromosome_label_policy, '; labels changed:', sum(df$chrom != df$source_chrom)), stderr())
unknown <- setdiff(unique(df$chrom), names(chrlen))
if (length(unknown)) stop(paste('Chromosome labels absent from the selected assembly:', paste(unknown, collapse = ', ')))
if (any(!is.finite(df$pos) | df$pos < 1 | df$pos != floor(df$pos) | df$pos > chrlen[df$chrom])) stop('Positions must be valid one-based coordinates in the selected assembly')
df$chrom <- factor(df$chrom, levels = names(chrlen))
# Index offsets by chromosome name, never by subset factor codes.
offsets <- setNames(c(0, head(cumsum(args$spacing * 1e6 + chrlen), -1)), names(chrlen))
df$chrompos <- unname(offsets[as.character(df$chrom)]) + df$pos
df$raw_lp <- df$lp

if (is.null(args$height)) {
  if (length(unique(df$chrom)) == 1) {
    args$height <- args$single_height
  } else {
    args$height <- args$genome_height
  }
}

if (!is.null(args$max_height)) {
  max_height <- args$max_height
  write(paste('Points above configured plot ceiling:', sum(df$raw_lp > max_height)), stderr())
} else {
  max_height <- max(df$lp, -log10(args$significance), -log10(args$suggestive))
}

# see http://github.com/FINNGEN/saige-pipelines/blob/master/scripts/qqplot.R
transform_lp <- function(x) {
  selected <- x > args$loglog_pval
  x[selected] <- args$loglog_pval * log10(x[selected]) / log10(args$loglog_pval)
  x
}
max_loglog <- transform_lp(max_height)
df$lp <- transform_lp(df$lp)
significance_y <- transform_lp(-log10(args$significance))
suggestive_y <- transform_lp(-log10(args$suggestive))
write(paste('Reference lines (raw P / plotted Y):', args$significance, significance_y, ';', args$suggestive, suggestive_y), stderr())
tick_pos <- round(seq(1, max_loglog, length.out = 10))
tick_lab <- sapply(tick_pos, function(x) { round(ifelse(x < args$loglog_pval, x, args$loglog_pval^(x/args$loglog_pval))) })

if (length(unique(df$chrom)) == 1) {
  p <- ggplot(df, aes(x = pos/1e6, y = lp)) +
    geom_hline(yintercept = significance_y, color = 'red', lty = 'longdash',linewidth = 0.3) +
    geom_hline(yintercept = suggestive_y, color = 'gray30', lty = 'longdash',linewidth = 0.4) +
    geom_point(size = 1/8) +
    scale_x_continuous(paste('Chromosome', unique(df$chrom), '(Mbp position)'), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', breaks = tick_pos, labels = tick_lab, expand = c(.01,.01)) +
    theme_bw(base_size = args$fontsize)
  
} else {
  p <- ggplot(df, aes(x = chrompos/1e6, y = lp, color = as.numeric(chrom)%%2 == 1)) +
    geom_hline(yintercept = significance_y, color = 'red', lty = 'longdash',linewidth = 0.3) +
    geom_hline(yintercept = suggestive_y, color = 'gray30', lty = 'longdash',linewidth = 0.4) +
    geom_point(size = 1/8) +
    scale_x_continuous(NULL, breaks = (cumsum(args$spacing * 1e6 + chrlen) - chrlen/2 - args$spacing * 1e6) / 1e6, labels = names(chrlen), guide = guide_axis(n.dodge = args$axis_label_rows), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', breaks = tick_pos, labels = tick_lab, expand = c(.01,.01)) +
    scale_color_manual(guide = 'none', values = c('FALSE' = 'dodgerblue', 'TRUE' = 'gray')) +
    theme_bw(base_size = args$fontsize)
}

caption <- args$caption_template
caption_values <- list(plotted_records = nrow(df), input_records = args$input_records, minimum_lp = args$min_lp, minimum_af = args$min_af, loglog_threshold = args$loglog_pval, empty_chromosomes = paste(setdiff(names(chrlen), unique(as.character(df$chrom))), collapse = ', '))
for (key in names(caption_values)) caption <- gsub(paste0('{', key, '}'), as.character(caption_values[[key]]), caption, fixed = TRUE)
write(paste('Plot caption:', caption), stderr())
p <- p + labs(caption = caption) +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(colour="black"),
    axis.line.y = element_line(colour="black"),
    plot.caption = element_text(size = args$caption_fontsize, hjust = 0)
  )

# this errors out when df is empty
if (args$csq) {
  df$coding <- FALSE
  for (annotation in strsplit(args$coding_terms, ',', fixed = TRUE)[[1]]) {
    df$coding <- df$coding | grepl(annotation, df$consequence, fixed = TRUE)
  }
  p <- p + geom_point(data = df[df$coding, ], size = 1/3, color = 'red')
}

if (!is.null(args$cytoband) && length(unique(df$chrom)) == 1) {
  cyto_height <- (max_loglog - args$min_lp) / args$cyto_ratio
  p <- p  +
    geom_rect(data = df_cyto[df_cyto$chrom == unique(df$chrom) & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = args$min_lp - cyto_height, ymax = args$min_lp, color = 'black', linewidth = 1/4, show.legend = FALSE) +
    scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
  if (exists('df_cen')) {
    p <- p + geom_polygon(data = df_cen[df_cen$chrom == unique(df$chrom),], aes(x = x/1e6, y = args$min_lp + cyto_height * y, shape = NULL, group = name), color = 'black', fill = 'red', linewidth = 1/8)
  }
} else {
  cyto_height <- 0
}

xlim <- FALSE
if ( !is.null(args$region) ) {
  if (grepl(':', args$region) && grepl('-', args$region)) {
    left <- as.numeric(gsub('^.*:', '', gsub('-.*$', '', args$region)))
    right <- as.numeric(gsub('^.*-', '', args$region))
    xlim <- TRUE
  }
}

if (!is.null(args$max_height) && xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6, ylim = c(transform_lp(args$min_lp) - cyto_height, max_loglog))
} else if (xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6)
} else if (!is.null(args$max_height)) {
  p <- p + coord_cartesian(ylim = c(transform_lp(args$min_lp) - cyto_height, max_loglog))
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = args$dpi)
}
print(p)
invisible(dev.off())
# Preserve the exact plotted values and untransformed LP for audit.
fwrite(df, args$plot_data, sep = '\t', na = 'NA')
