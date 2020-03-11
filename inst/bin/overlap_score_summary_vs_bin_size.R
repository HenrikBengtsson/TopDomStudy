#!/usr/bin/env Rscript
##-----------------------------------------------------------------------------
## Usage:
##
##  Rscript overlap_score_summary_vs_bin_size.R
##
## Parameters:
##  * chromosomes
##  * bin_sizes
##  * rhos
##  * nsamples
##  * window_size
##  * weights
##  * domain_length
##
## Input:
##  * human,HAP1,unique,chr*.rds files in R package 'TopDomStudy', i.e.
##    in folder system.file("compiledData", package="TopDomStudy")
##
## Output:
##  * figures/
##  * overlapScoreData/
##  * overlapScoreSummary/
##-----------------------------------------------------------------------------
library(TopDomStudy)
library(future) ## sources ./.future.R, if it exists

## Allow for 3-GiB objects to be exported during parallelizing
options(future.globals.maxSize = 3 * 1024^3)

bin_sizes <- c(5, 6, 8, 10, 12, 15, 20, 40, 60, 80, 100) * 1e3
bin_sizes <- c(         10,         20, 40, 60, 80, 100) * 1e3
bin_sizes <- c(   6, 8, 10, 12, 15, 20, 40, 60, 80, 100) * 1e3
bin_sizes <- c(   6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 80, 100) * 1e3
## FIXME: chromosome = "22", rho = 0.01, bin_size = 10000, nsamples = 1L
## gives an error
rhos <- c(0.01, 0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50)[-1]
rhos <- c(      0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50) ## Skip rho=0.02 due to chr=12 memory constraints
rhos <- c(      0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50) ## Skip rho=0.02 due to chr=12 memory constraints
reference_rhos <- c("50%", "same")[1]

domain_length <- NULL
#domain_length <- c(300e3, 1000e3)

done <- overlap_score_summary_vs_bin_size(
  dataset        = "human,HAP1",
  chromosomes    = rev(c("12", "16", "22")),
  bin_sizes      = rev(bin_sizes),
  rhos           = rhos,
  reference_rhos = reference_rhos,
  window_size    = 5L,
  weights        = c("by_length", "uniform")[1],
  domain_length  = domain_length,
  nsamples       = 50L,
  verbose        = TRUE
)
print(done)
