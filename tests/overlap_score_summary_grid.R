library(TopDomStudy)
library(future)

plan(list(
  chr_bin_rho = multiprocess, ## overlap_score_summary_grid()
  mono_chr    = sequential,   ## +-- overlap_scores_partitions()
  samples     = sequential    ##     - " -
))

chromosomes <- getOption("TopDomStudy.tests.chromosomes", "22")
bin_sizes <- getOption("TopDomStudy.tests.bin_sizes", c(50e3, 100e3, 200e3))
rhos <- getOption("TopDomStudy.tests.rhos", c(0.10, 0.25, 0.50))
nsamples <- getOption("TopDomStudy.tests.nsamples", 5L)

for (weights in c("uniform", "by_length")) {
  pathnames <- overlap_score_summary_grid(
    dataset       = "human,HAP1",
    chromosomes   = chromosomes,
    bin_sizes     = bin_sizes,
    rhos          = rev(rhos)[1],
    window_size   = 5L,
    weights       = weights,
    domain_length = c(500e3, 1000e3),
    nsamples      = nsamples,
    verbose       = TRUE
  )
  print(pathnames)
}
