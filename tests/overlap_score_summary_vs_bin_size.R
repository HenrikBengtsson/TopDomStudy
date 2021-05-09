source("incl/start.R")

## NOTE:
## These tests consumes a large amount of RAM (> 20 GiB)
## /HB 2021-05-09

library(future)
plan(list(
  chr_bin_rho = multisession, ## overlap_score_summary_grid()
  mono_chr    = sequential,   ## +-- overlap_scores_partitions()
  samples     = sequential    ##     - " -
))

for (reference_rhos in c("50%", "same")) {

  for (weights in c("uniform", "by_length")) {
    pathnames <- overlap_score_summary_vs_bin_size(
      dataset        = dataset,
      chromosomes    = chromosomes,
      bin_sizes      = bin_sizes,
      rhos           = rhos,
      reference_rhos = reference_rhos,
      window_size    = 5L,
      weights        = weights,
      domain_length  = c(500e3, 1000e3),
      nsamples       = nsamples,
      verbose        = TRUE
    )
    print(pathnames)
  }

}

source("incl/end.R")
