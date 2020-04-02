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
cmdArg <- R.utils::cmdArg

## Allow for 3-GiB objects to be exported during parallelizing
options(future.globals.maxSize = 3 * 1024^3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parse command-line options
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chromosomes <- cmdArg(chromosomes = c("12", "16", "22"))

bin_sizes <- cmdArg(bin_sizes = c(6, 8, 10, 12, 15, 20,
                                  30, 40, 50, 60, 80, 100) * 1e3)

## FIXME: chromosome="22", rho=0.01, bin_size=10000, nsamples=1L gives an error
## Skip rho=0.02 due to chr=12 memory constraints
rhos <- cmdArg(rhos = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.08,
                        0.10, 0.12, 0.14, 0.16, 0.18, 0.20,
			0.25, 0.30, 0.40, 0.50))

choices <- c("50%", "same")
reference_rhos <- match.arg(cmdArg(reference_rhos = choices[1]), choices)

window_size <- cmdArg(window_size = 5L)

choices <- c("by_length", "uniform")
weights <- match.arg(cmdArg(weights = choices[1]), choices)

#domain_length <- "ref_len_iqr"  ## WARNING: Requires that domain_length = NULL has been run before
#domain_length <- c(300e3, 1000e3)
domain_length <- cmdArg(domain_length = NULL)

nsamples <- cmdArg(nsamples = 100L)

verbose <- cmdArg(verbose = TRUE)

choices <- c("bin_size", "fraction")
vs <- match.arg(cmdArg(vs = choices[1]), choices)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
FUN <- switch(vs,
  bin_size = overlap_score_summary_vs_bin_size,
  fraction = overlap_score_summary_vs_fraction
)

done <- FUN(
  dataset        = "human,HAP1",
  chromosomes    = chromosomes,
  bin_sizes      = bin_sizes,
  rhos           = rhos,
  reference_rhos = reference_rhos,
  window_size    = window_size,
  weights        = weights,
  domain_length  = domain_length,
  nsamples       = nsamples,
  verbose        = verbose
)
print(done)
