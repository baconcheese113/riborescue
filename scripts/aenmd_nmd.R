#!/usr/bin/env Rscript

# Run aenmd (Klonowski et al., Bioinformatics 2023) over a table of nonsense SNVs and write its
# per-transcript NMD-escape rules to a TSV. aenmd annotates every Ensembl v105 transcript a variant
# overlaps; the Python side keeps only the MANE Select transcript, so the model tier reads NMD on the
# same molecule as the rule tier. Output is one row per variant x transcript, the rules aenmd fired.

suppressPackageStartupMessages({
  library(optparse)
  library(aenmd)
  library(GenomicRanges)
  library(S4Vectors)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--variants", type = "character", help = "TSV with chrom,pos,ref,alt columns"),
  make_option("--out", type = "character", help = "output TSV of per-transcript aenmd rules")
)))

vars <- read.delim(opt$variants, stringsAsFactors = FALSE, colClasses = "character")
vars <- vars[nchar(vars$ref) == 1L & nchar(vars$alt) == 1L, ] # aenmd's NMD rules score SNVs here
vars$pos <- as.integer(vars$pos)
message(sprintf("read %d SNVs", nrow(vars)))

# A minimal, coordinate-sorted VCF in NCBI (no-'chr') naming, which BSgenome.Hsapiens.NCBI.GRCh38 —
# aenmd's genome — expects. aenmd's own parser then reads it into the GRanges the pipeline wants.
ord <- order(factor(vars$chrom, levels = c(1:22, "X", "Y", "MT")), vars$pos)
vars <- vars[ord, ]
vcf_path <- tempfile(fileext = ".vcf")
con <- file(vcf_path, "w")
writeLines(c(
  "##fileformat=VCFv4.2",
  paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", sep = "\t")
), con)
writeLines(
  paste(vars$chrom, vars$pos, ".", vars$ref, vars$alt, ".", ".", ".", sep = "\t"), con
)
close(con)

vcf <- aenmd:::parse_vcf_VariantAnnotation(vcf_path)
vcf_rng <- vcf$vcf_rng
message(sprintf("parsed %d variants", length(vcf_rng)))

vcf_rng_fil <- process_variants(vcf_rng)
message(sprintf("%d variants pass aenmd filtering (coding, non-splice, stop-generating)", length(vcf_rng_fil)))

ann <- annotate_nmd(vcf_rng_fil, rettype = "gr")
res <- ann$res_aenmd
message(sprintf("%d variant x transcript annotations", nrow(res)))

out <- data.frame(
  key = ann$key,
  transcript = res$transcript,
  is_ptc = res$is_ptc,
  is_last = res$is_last,
  is_penultimate = res$is_penultimate,
  is_css_proximal = res$is_cssProximal,
  is_single = res$is_single,
  is_407plus = res$is_407plus,
  stringsAsFactors = FALSE
)
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("wrote %s", opt$out))
