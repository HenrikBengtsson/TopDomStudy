library(TopDomStudy)
library(future)

plan(list(
  chr_bin  = sequential,
  rho      = sequential,
  monochr  = sequential,
  samples  = multiprocess
))

for (weights in c("uniform", "by_length")) {
  pathnames <- overlap_score_summary_vs_fraction(
    dataset       = "human,HAP1",
    chromosomes   = "22",
    bin_sizes     = c(50e3, 100e3, 200e3),
    rhos          = 0.25,
    window_size   = 5L,
    weights       = weights,
    domain_length = c(500e3, 1000e3),
    nsamples      = 5L,
    verbose       = TRUE
  )
  print(pathnames)
}
