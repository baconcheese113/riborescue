#!/usr/bin/env Rscript

# Per-drug replicate reliability, the ceiling every held-out score is read against.
#
# A model at r-squared 0.80 against a ceiling of 0.89 is a different claim from the same number
# against a ceiling of 1.0, and the ceiling is not one number: the paper's 0.89 is a pan-drug figure,
# while its own text puts CC90009's inter-replicate correlation far below that. Reported per drug,
# from the authors' own two replicates of every variant.
#
# Two things decide what the number is, and getting either wrong inverts the conclusion.
#
# **The correlation is taken on `fitness_rep`.** `RT` and `RT_binomial` are merged values, identical
# across a variant's two replicate rows, so correlating either returns exactly 1 and reports perfect
# reproducibility. The assertion below fires if a future export makes the response per-replicate.
#
# **The ceiling is the reliability of the mean of two replicates, not of one.** The modelled response
# is the merged value — `RT_binomial` is an exact linear function of the two replicates' mean, which
# this script checks rather than assumes — and a mean of two measurements is more reliable than
# either. The single-measure correlation is corrected by Spearman-Brown, `2r / (1 + r)`. Without that
# correction three of the six drugs come out *above* their own ceiling, which is not a finding about
# those drugs but an error about which quantity was being predicted.

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

options <- parse_args(OptionParser(option_list = list(
  make_option("--samples", default = "data/toledano/treated_samples.rds",
              help = "The authors' per-replicate measurements [default %default]"),
  make_option("--out", default = "tests/fixtures/oracle/reliability.tsv",
              help = "Where to write the per-drug ceiling [default %default]")
)))

samples <- as.data.table(readRDS(options$samples))
stopifnot(all(c("treatment", "replicate", "identifier", "fitness_rep") %in% names(samples)))
stopifnot(!anyDuplicated(samples[, .(treatment, replicate, identifier)]))

wide <- dcast(samples, treatment + identifier ~ replicate,
              value.var = c("fitness_rep", "RT_binomial"))
setnames(wide, c("fitness_rep_1", "fitness_rep_2"), c("rep1", "rep2"))
paired <- wide[!is.na(rep1) & !is.na(rep2) & !is.na(RT_binomial_1)]

# The response is the merged value, and the merged value is the replicates' mean on another scale.
stopifnot(isTRUE(all.equal(paired$RT_binomial_1, paired$RT_binomial_2)))
stopifnot(abs(cor((paired$rep1 + paired$rep2) / 2, paired$RT_binomial_1)) > 1 - 1e-6)

reliability <- paired[, {
  single <- cor(rep1, rep2)
  .(variants = .N, r = single, single_replicate = single^2,
    ceiling = 2 * single / (1 + single))
}, by = treatment][order(treatment)]

fwrite(reliability, options$out, sep = "\t")
print(reliability)
