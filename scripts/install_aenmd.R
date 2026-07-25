#!/usr/bin/env Rscript

# Install aenmd (Klonowski et al., Bioinformatics 2023; kostkalab/aenmd, BSD-2-Clause) and its
# Ensembl-v105 GRCh38 transcript-model data package from pinned commits. The heavy Bioconductor
# dependencies come from bioconda and are captured in pixi.lock; this script adds only the two
# GitHub-only packages and never upgrades a bioconda-provided dependency (upgrade = "never").

aenmd_sha <- "7aea66992504ed28d2d6d4001e4be5bec9183aef"
data_sha <- "2d5db4a4d532740c14215a3eb41a1e75f573007e"

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_pinned <- function(spec, ref) {
  remotes::install_github(spec, ref = ref, upgrade = "never", dependencies = NA, quiet = FALSE)
}

if (!requireNamespace("aenmd.data.ensdb.v105", quietly = TRUE)) {
  install_pinned("kostkalab/aenmd_data/aenmd.data.ensdb.v105", data_sha)
}
if (!requireNamespace("aenmd", quietly = TRUE)) {
  install_pinned("kostkalab/aenmd", aenmd_sha)
}

suppressPackageStartupMessages({
  library(aenmd)
  library(aenmd.data.ensdb.v105)
})
cat(
  "aenmd", as.character(packageVersion("aenmd")),
  "+ data", as.character(packageVersion("aenmd.data.ensdb.v105")), "ready\n"
)
