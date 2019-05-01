library(TopDomStudy)
library(future)

plan(list(
  chr_bin  = sequential,
  rho      = sequential,
  monochr  = sequential,
  samples  = multiprocess
))

for (weights in c("uniform", "by_length")) {
  done <- overlap_score_summary_vs_fraction(
    dataset       = "human,HAP1",
    chromosomes   = "22",
    bin_sizes     = 100e3,
    rhos          = c(0.05, 0.20, 0.50),
    window_size   = 5L,
    weights       = weights,
    domain_length = c(500e3, 1000e3),
    nsamples      = 10L,
    verbose       = TRUE
  )
  print(done)
}
