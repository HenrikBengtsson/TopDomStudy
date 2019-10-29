%\VignetteIndexEntry{TopDomStudy: Calculating TopDom Overlap Scores and their Summaries}
%\VignetteAuthor{Henrik Bengtsson}
%\VignetteEngine{TopDomStudy::self}
%\VignetteEncoding{UTF-8}
%\VignetteKeyword{R}
%\VignetteKeyword{package}
%\VignetteKeyword{vignette}


In this document we show how to precalculate and summarize TopDom overlap across combination of fraction and bin-size parameter settings.


## Precalculate and Summarize TopDom Overlap Scores

To precalculate and summarize TopDom overlap scores across the full 10-by-10 grid of 10 different sample fractions and 10 different bin sizes for chromosomes 12, 16, and 22, run:

```r
library(TopDomStudy)
## Allow for 2-GiB objects to be exported to parallel workers
options(future.globals.maxSize = 2*1024^3)

done <- overlap_score_summary_grid(
  dataset       = "human,HAP1",
  chromosomes   = c("12", "16", "22"),
  bin_sizes     = c(6, 8, 10, 12, 15, 20, 40, 60, 80, 100) * 1e3,
  rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 50L,
  verbose       = TRUE
)
print(done)
, , 0.02

   6000                                                                                                                                      
# 12 "overlapScoreSummary/human,HAP1,chr=12,cells_by_half,avg_score,bin_size=6000,fraction=0.020,window_size=5,weights=uniform,nsamples=50.rds"
# 16 "overlapScoreSummary/human,HAP1,chr=16,cells_by_half,avg_score,bin_size=6000,fraction=0.020,window_size=5,weights=uniform,nsamples=50.rds"
# 22 "overlapScoreSummary/human,HAP1,chr=22,cells_by_half,avg_score,bin_size=6000,fraction=0.020,window_size=5,weights=uniform,nsamples=50.rds"
# [ ... another 3*98 entries ... ]
# 12 "overlapScoreSummary/human,HAP1,chr=12,cells_by_half,avg_score,bin_size=100000,fraction=0.500,window_size=5,weights=uniform,nsamples=50.rds"
# 16 "overlapScoreSummary/human,HAP1,chr=16,cells_by_half,avg_score,bin_size=100000,fraction=0.500,window_size=5,weights=uniform,nsamples=50.rds"
# 22 "overlapScoreSummary/human,HAP1,chr=22,cells_by_half,avg_score,bin_size=100000,fraction=0.500,window_size=5,weights=uniform,nsamples=50.rds"
```

The overlap scores and their summaries are written to file to folders `./overlapScoreData/` and `./overlapScoreSummaries/`, respectively.  Because of this file cache, already processed grid points will be skipped if rerun (in the same R session or in a future R session).

The above takes approximately 24-30 hours and up to 20 GiB of RAM to complete when running in parallel using two cores (with `future::plan("multicore", workers=2)`) on a Lenovo Thinkpad Carbon X1 (gen 6) with 16 GiB of RAM.  The majority of the processing time, and memory, is consumed on Chr 12 at the higher resolutions (bin sizes <= 10 kb).  Running sequentially (default) will lower the memory requirements to approximately 10 GiB of RAM.


In addition to the above computationally expensive 10-by-10 grid, we will also need a a few additional calcuations;

```r
library(TopDomStudy)
done <- unlist(lapply(c(5, 50, 100), FUN = function(nsamples) {
  overlap_score_summary_grid(
    dataset       = "human,HAP1",
    chromosomes   = "22",
    bin_sizes     = 100 * 1e3,
    rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
    window_size   = 5L,
    weights       = "uniform",
    domain_length = NULL,
    nsamples      = nsamples,
    verbose       = TRUE
  )
}))
print(done)
```

```r
library(TopDomStudy)
done <- unlist(lapply(c("uniform", "by_length"), FUN = function(weights) {
  overlap_score_summary_grid(
    dataset       = "human,HAP1",
    chromosomes   = "22",
    bin_sizes     = 100 * 1e3,
    rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
    window_size   = 5L,
    weights       = weights,
    domain_length = NULL,
    nsamples      = 50L,
    verbose       = TRUE
  )
}))
print(done)
```

```r
library(TopDomStudy)
done <- unlist(lapply(c(5, 10, 20, 40), FUN = function(window_size) {
  overlap_score_summary_grid(
    dataset       = "human,HAP1",
    chromosomes   = "22",
    bin_sizes     = 100 * 1e3,
    rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
    window_size   = window_size,
    weights       = "uniform",
    domain_length = NULL,
    nsamples      = 50L,
    verbose       = TRUE
  )
}))
print(done)
```



## Producing figures

To produce figures, as PNG images in `./figures/`, presenting the relationship between overlap-score summaries as a function of sample fraction used, run:

```r
library(TopDomStudy)

done <- overlap_score_summary_vs_fraction(
  dataset       = "human,HAP1",
  chromosomes   = c("12", "16", "22"),
  bin_sizes     = c(6, 8, 10, 12, 15, 20, 40, 60, 80, 100) * 1e3,
  rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 50L,
  fig_path      = "figures/",
  verbose       = TRUE
)
print(done)
```

To produce figures, as PNG images in `./figures/`, presenting the relationship between overlap-score summaries as a function of sample bin size used, run:

```r
library(TopDomStudy)

done <- overlap_score_summary_vs_bin_size(
  dataset       = "human,HAP1",
  chromosomes   = c("12", "16", "22"),
  bin_sizes     = c(6, 8, 10, 12, 15, 20, 40, 60, 80, 100) * 1e3,
  rhos          = c(0.02, 0.04, 0.05, 0.06, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50),
  window_size   = 5L,
  weights       = "uniform",
  domain_length = NULL,
  nsamples      = 50L,
  fig_path      = "figures/",
  verbose       = TRUE
)
print(done)
```

The produces figures use the default gray [ggplot2] theme (`theme_set(theme_gray())`).  To use produce more publication friendly figures, call `ggplot2::theme_set(cowplot::theme_cowplot())` prior to the above.


## References

1. Ramani, V., Deng, X., Qiu, R., Gunderson, K. L., Steemers, F. J., Disteche, C. M., … Shendure, J. (2017). Massively multiplex single-cell Hi-C. Nature methods, 14(3), 263–266. doi:10.1038/nmeth.4155, [PMC5330809](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330809/)


[ggplot2]: https://cran.r-project.org/package=ggplot2
