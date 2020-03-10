#' Coerces `hic_bin()` Data into list of TopDomData Objects
#' 
#' @param data A named list.
#'
#' @return A named list of length `length(data)`.
#' 
#' @export
as_TopDomData <- function(data) {
  stopifnot(is.list(data))
  chrs <- names(data)
  stopifnot(!is.null(chrs))

  res <- vector("list", length = length(data))
  
  for (kk in seq_along(data)) {
    chr <- chrs[kk]
    counts <- data[[kk]]
    bins <- attr(counts, "bins", exact = TRUE)
    stopifnot(!is.null(bins))
    attr(counts, "bins") <- NULL

    ## Calculate bin midpoints
    n <- length(bins$chr_a)
    loci <- data.frame(
      chr        = rep(chr, times = n - 1L),
      from.coord = (bins$chr_a[-1L] + bins$chr_a[-n]) / 2,
      to.coord   = (bins$chr_b[-1L] + bins$chr_b[-n]) / 2,
      stringsAsFactors = FALSE
    )

    stopifnot(nrow(loci) == nrow(counts), nrow(counts) == ncol(counts))
    res_kk <- structure(list(bins = loci, counts = counts), class = "TopDomData")
    res[[kk]] <- res_kk
  }
  
  names(res) <- names(data)
  
  res
}
