pathnames <- overlap_score_summary_grid(
  dataset       = "human,HAP1",
  chromosomes   = "22",
  bin_sizes     = 100*1e3,
  rhos          = c(0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 10L,
  verbose       = TRUE
)
print(pathnames)

gg <- gg_overlap_score_summary_vs_fraction(
  dataset       = "human,HAP1",
  chromosome    = "22",
  bin_size      = 100*1e3,
  rhos          = c(0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 10L,
  verbose       = TRUE
)
print(gg)
