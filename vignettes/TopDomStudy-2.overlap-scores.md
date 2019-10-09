%\VignetteIndexEntry{TopDomStudy: TopDom Overlap Score Summaries}
%\VignetteAuthor{Henrik Bengtsson}
%\VignetteEngine{TopDomStudy::self}
%\VignetteEncoding{UTF-8}
%\VignetteKeyword{R}
%\VignetteKeyword{package}
%\VignetteKeyword{vignette}


In this document we show how to precalculate TopDom overlap-score summaries for a combination of fractions and bin sizes.


## Overlap Scores versus Fraction

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
  verbose       = TRUE
)
print(done)
```

_Comment_: The above takes approximately NNN hours to complete.



## Overlap Scores versus Bin Size

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
  verbose       = TRUE
)
print(done)
```

_Comment_: The above takes approximately NNN hours to complete.


## References

1. Ramani, V., Deng, X., Qiu, R., Gunderson, K. L., Steemers, F. J., Disteche, C. M., … Shendure, J. (2017). Massively multiplex single-cell Hi-C. Nature methods, 14(3), 263–266. doi:10.1038/nmeth.4155, [PMC5330809](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5330809/)
