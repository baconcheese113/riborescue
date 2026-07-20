#!/usr/bin/env Rscript

# The reproduction oracle: the authors' own drug-specific readthrough model, run against their
# published data, exported as golden fixtures. Everything downstream is tested for parity against
# these artifacts rather than against the summary statistics printed in the paper — identical
# predictions on identical folds is the only evidence that counts.
#
# The model, the eligibility filter, the seed and the partition scheme are transcribed from
# Fig4_extdataFig4.Rmd of lehner-lab/Stop_codon_readthrough and must not be "improved" here.

suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(jsonlite)
  library(optparse)
})

FORMULA <- RT_binomial ~ 0 + stop_type + down_123nt + up_123nt + stop_type:down_123nt
SEED <- 721
ROUNDS <- 10
TRAIN_FRACTION <- 0.9

# The six drugs the paper models: those triggering >1% readthrough for >3% of PTCs.
DRUGS <- c("CC90009", "Clitocine", "DAP", "G418", "SJ6986", "SRI")

options <- parse_args(OptionParser(option_list = list(
  make_option("--input", default = file.path(
    Sys.getenv("RIBORESCUE_DATA", "data"), "toledano", "treated_samples.rds"
  ), help = "The published readthrough measurements [default %default]"),
  make_option("--out", default = file.path("tests", "fixtures", "oracle"),
    help = "Where the golden fixtures are written [default %default]")
)))
input <- options$input
outdir <- options$out

if (!file.exists(input)) {
  stop(sprintf("%s is absent; fetch it with `riborescue fetch toledano_treated_samples`", input))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
treated_samples <- readRDS(input)

# The authors' eligibility filter: one replicate, adequately covered, non-viral, measured.
eligible_for <- function(drug) {
  treated_samples[
    treatment == drug & replicate == 2 & reads_allbins > 15 & viral == "no" & !is.na(RT)
  ]
}

metrics <- data.table()
for (drug in DRUGS) {
  eligible <- eligible_for(drug)
  set.seed(SEED)
  partitions <- createDataPartition(
    eligible$RT_binomial, p = TRAIN_FRACTION, list = FALSE, times = ROUNDS
  )

  features <- eligible[, .(
    row = .I, stop_type, up_123nt, down_123nt, RT_binomial
  )]
  fwrite(features, file.path(outdir, sprintf("features_%s.tsv.gz", drug)), sep = "\t")

  folds <- rbindlist(lapply(seq_len(ncol(partitions)), function(round) {
    data.table(round = round, row = partitions[, round])
  }))
  fwrite(folds, file.path(outdir, sprintf("folds_%s.tsv.gz", drug)), sep = "\t")

  coefficients <- data.table()
  predictions <- data.table()
  r2 <- numeric()
  for (round in seq_len(ncol(partitions))) {
    train <- eligible[partitions[, round], ]
    test <- eligible[-partitions[, round], ]
    model <- glm(FORMULA, train, family = "binomial")
    predicted <- predict(model, test, type = "response")

    coefficients <- rbind(coefficients, data.table(
      round = round, term = names(coef(model)), estimate = unname(coef(model))
    ))
    predictions <- rbind(predictions, data.table(
      round = round,
      row = setdiff(seq_len(nrow(eligible)), partitions[, round]),
      predicted = unname(predicted),
      observed = test$RT_binomial
    ))
    r2 <- c(r2, cor(predicted, test$RT_binomial)^2)
  }
  fwrite(coefficients, file.path(outdir, sprintf("coefficients_%s.tsv.gz", drug)), sep = "\t")
  fwrite(predictions, file.path(outdir, sprintf("predictions_%s.tsv.gz", drug)), sep = "\t")
  metrics <- rbind(metrics, data.table(drug = drug, round = seq_along(r2), r2 = r2))
  cat(sprintf("%-10s n=%d  mean r2=%.4f\n", drug, nrow(eligible), mean(r2)))
}

fwrite(metrics, file.path(outdir, "metrics.tsv"), sep = "\t")

provenance <- list(
  source = "lehner-lab/Stop_codon_readthrough, Fig4_extdataFig4.Rmd",
  input = basename(input),
  input_md5 = as.character(tools::md5sum(input)),
  formula = paste(deparse(FORMULA), collapse = " "),
  family = "binomial",
  seed = SEED,
  rounds = ROUNDS,
  train_fraction = TRAIN_FRACTION,
  drugs = DRUGS,
  r_version = R.version.string,
  caret = as.character(packageVersion("caret")),
  data_table = as.character(packageVersion("data.table"))
)
write_json(provenance, file.path(outdir, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
