#!/usr/bin/env Rscript

# P-site calibration and periodicity, the Gate 3 check, computed with riboWaltz — the published tool
# the PRD names. riboWaltz infers a coherent offset for every footprint length from the pile-up of
# read ends at annotated start codons, so nothing is hardcoded. Its input is a transcriptome-
# coordinate alignment, which places the P-site in transcript space and never has to reason about
# exon junctions.
#
# One sample per invocation. A transcriptome alignment carries every isoform a footprint hits — the
# secondary records outnumber the primary ones roughly thirty to one — so the whole library will not
# fit in memory at once. Each run keeps the primary alignment of each read, calibrates, writes its
# tables and exits, which returns the memory before the next sample starts. `--combine` merges the
# per-sample tables afterwards.

suppressPackageStartupMessages({
  library(riboWaltz)
  library(Rsamtools)
  library(data.table)
  library(optparse)
})

options <- parse_args(OptionParser(option_list = list(
  make_option("--bam", help = "One *.toTranscriptome.out.bam footprint alignment"),
  make_option("--sample", help = "Name to record this library under"),
  make_option("--gtf", help = "The GENCODE annotation the alignments were made against"),
  make_option("--outdir", default = "results/reads/qc/psite", help = "[default %default]"),
  make_option("--lengths", default = "26:34", help = "Footprint length range [default %default]"),
  make_option("--combine", action = "store_true", default = FALSE,
              help = "Merge the per-sample tables already in --outdir"),
  make_option("--evidence", default = "docs/figures",
              help = "Where the tracked Gate 3 evidence lands [default %default]")
)))

dir.create(options$outdir, recursive = TRUE, showWarnings = FALSE)
part <- function(kind, sample) file.path(options$outdir, sprintf("%s.%s.tsv", sample, kind))

if (options$combine) {
  gather <- function(kind) {
    files <- list.files(options$outdir, pattern = sprintf("\\.%s\\.tsv$", kind), full.names = TRUE)
    stopifnot(length(files) > 0)
    rbindlist(lapply(files, fread))
  }
  offsets <- gather("offsets")
  frames <- gather("frame")
  fwrite(offsets, file.path(options$outdir, "psite_offsets.tsv"), sep = "\t")
  fwrite(frames, file.path(options$outdir, "frame_by_length.tsv"), sep = "\t")
  fwrite(gather("region"), file.path(options$outdir, "region_distribution.tsv"), sep = "\t")
  fwrite(gather("shift"), file.path(options$outdir, "offset_shift_check.tsv"), sep = "\t")
  fwrite(gather("metaprofile"), file.path(options$outdir, "metaprofile.tsv"), sep = "\t")

  frame_cds <- frames[, .(n = sum(n)), by = .(sample, frame)]
  frame_cds[, fraction := n / sum(n), by = sample]
  fwrite(frame_cds, file.path(options$outdir, "frame_cds.tsv"), sep = "\t")

  # The offset is reported from both ends. The familiar ~12 nt figure is the distance from the 5'
  # end; the 3' figure is what riboWaltz optimises and the two are not interchangeable.
  dominant <- offsets[, .SD[which.max(total_percentage)], by = sample]
  summary <- merge(
    dominant[, .(sample, dominant_length = length,
                 offset_from_5 = corrected_offset_from_5,
                 offset_from_3 = corrected_offset_from_3)],
    frame_cds[frame == 0, .(sample, cds_frame0 = round(fraction, 3))],
    by = "sample"
  )
  fwrite(summary, file.path(options$outdir, "psite_summary.tsv"), sep = "\t")

  # Gate 3 is a judgement a person makes, so the evidence behind it is tracked rather than left in
  # the regenerable results tree: the table that was read, and the periodicity that was looked at.
  dir.create(options$evidence, recursive = TRUE, showWarnings = FALSE)
  fwrite(summary, file.path(options$evidence, "gate3_psite_summary.tsv"), sep = "\t")

  profile <- gather("metaprofile")
  profile[, region := factor(region,
                             levels = c("Distance from start (nt)", "Distance from stop (nt)"),
                             labels = c("start codon", "stop codon"))]
  figure <- ggplot2::ggplot(profile[!is.na(region)],
                            ggplot2::aes(x = distance, y = scaled_count)) +
    ggplot2::geom_line(linewidth = 0.3, colour = "#c0392b") +
    ggplot2::facet_grid(sample ~ region, scales = "free", switch = "y") +
    ggplot2::labs(x = "distance from codon (nt)", y = "P-site frequency") +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1))
  ggplot2::ggsave(file.path(options$evidence, "gate3_periodicity.png"),
                  figure, width = 9, height = 10, dpi = 150)

  print(summary)
  quit(save = "no")
}

stopifnot(!is.null(options$bam), !is.null(options$sample), !is.null(options$gtf))
length_range <- eval(parse(text = options$lengths))

# Building the annotation walks the whole GTF, so it is built once and reused by later samples.
cache <- file.path(options$outdir, "annotation.rds")
annotation <- if (file.exists(cache)) {
  readRDS(cache)
} else {
  built <- create_annotation(gtfpath = options$gtf)
  saveRDS(built, cache)
  built
}

# One alignment per read, so a footprint is counted once rather than once per isoform it happens to
# fit. STAR's primary placement is a pragmatic assignment for calibration, not a claim about which
# isoform the footprint truly came from; offsets and frame are read relative to whichever coding
# sequence it was assigned, and isoform identity is not inferred here.
staged <- file.path(tempdir(), "psite_bam")
dir.create(staged, showWarnings = FALSE)
primary <- file.path(staged, paste0(options$sample, ".bam"))
filterBam(
  options$bam, primary, index = character(0), indexDestination = FALSE,
  param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE))
)

reads <- bamtolist(bamfolder = staged, annotation = annotation)
reads <- length_filter(reads, length_filter_mode = "custom", length_range = length_range)
offsets <- psite(reads, flanking = 6, extremity = "auto")
fwrite(offsets, part("offsets", options$sample), sep = "\t")

reads <- psite_info(reads, offsets)
sample_reads <- reads[[options$sample]]

# Three-frame distribution over the coding sequence, per length. Reported in full — the fractions
# speak for themselves, without a pass threshold imposed here.
frame_len <- sample_reads[psite_region == "cds",
                          .(sample = options$sample, n = .N),
                          by = .(length, frame = psite_from_start %% 3)]
fwrite(frame_len, part("frame", options$sample), sep = "\t")

# Where the P-sites fall: 5' UTR, CDS, 3' UTR. A footprint library sits overwhelmingly in the CDS.
region <- sample_reads[, .(sample = options$sample, n = .N), by = .(region = psite_region)]
region[, fraction := n / sum(n)]
fwrite(region, part("region", options$sample), sep = "\t")

# An internal consistency check, not held-out validation: it re-scores the same reads the offset was
# inferred from. A resolved offset should still beat its immediate neighbours on coding-frame
# enrichment; a flat result would mean the offset is not determined by the data.
dominant <- offsets[which.max(total_percentage)]
in_cds <- sample_reads[length == dominant$length & psite_region == "cds"]
shift <- rbindlist(lapply(-2:2, function(k) {
  data.table(sample = options$sample, length = dominant$length, shift = k,
             frame0 = mean(((in_cds$psite_from_start + k) %% 3) == 0))
}))
fwrite(shift, part("shift", options$sample), sep = "\t")

# Start and stop metaprofiles. Three-nucleotide periodicity has to be visible, not only implied by a
# frame fraction, so the profile is kept as counts and drawn.
meta <- metaprofile_psite(reads, annotation, sample = setNames(list(options$sample), options$sample),
                          utr5l = 25, cdsl = 50, utr3l = 25)
fwrite(meta$count_dt, part("metaprofile", options$sample), sep = "\t")

# The returned list carries both a table named plot_dt and the drawing itself, so the plot is picked
# by class rather than by name.
drawn <- Filter(function(x) inherits(x, "ggplot"), meta)
ggplot2::ggsave(file.path(options$outdir, paste0(options$sample, ".metaprofile.png")),
                drawn[[1]], width = 10, height = 4, dpi = 120)

cat(sprintf("%s: %d nt dominant, offset %g from 5' / %g from 3'\n", options$sample,
            dominant$length, dominant$corrected_offset_from_5, dominant$corrected_offset_from_3))
