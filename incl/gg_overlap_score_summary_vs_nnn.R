## These calls with produce folders in the working directory

gg1 <- gg_overlap_score_summary_vs_fraction(
  dataset       = "human,HAP1",
  chromosome    = "22",
  bin_size      = 100*1e3,
  rhos          = c(0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 4L,
  line_col      = "red",
  verbose       = TRUE
)
print(gg1)

gg2 <- gg_overlap_score_summary_vs_bin_size(
  dataset       = "human,HAP1",
  chromosome    = "22",
  bin_sizes     = c(20, 40, 60, 80, 100)*1e3,
  rho           = 0.50,
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 4L,
  line_col      = "blue",
  verbose       = TRUE
)
print(gg2)
