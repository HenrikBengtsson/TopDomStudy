#' @importFrom utils file_test
extract_domain_length_limits <- function(dataset, chromosome = ".*", bin_size = ".*", nsamples, weights = c("by_length", "uniform"), verbose = FALSE) {
  weights <- match.arg(weights)
  
  path <- "overlapScoreSummary"
  stopifnot(file_test("-d", path))
  pattern <- sprintf("%s,chr=%s,cells_by_half,avg_score-vs-fraction,bin_size=%s,test=[.0-9]+,nsamples=%d,weights=%s.rds", dataset, chromosome, bin_size, nsamples, weights)
##  if (verbose) mprint(pattern)
  pathnames <- dir(path, pattern = pattern, full.names = TRUE)
##  if (verbose) mprint(pathnames)
  stopifnot(length(pathnames) > 0L)

  data <- lapply(pathnames, FUN = function(pathname) {
    summary <- read_rds(pathname)
    ref_len_qs <- vapply(grep("^ref_len_q", names(summary), value = TRUE), FUN = function(name) mean(summary[[name]], na.rm = TRUE), FUN.VALUE = NA_real_)
    summary_mu <- cbind(summary[1L, c("chromosome", "bin_size", "fraction")], nsamples = nrow(summary), as.list(ref_len_qs))
    summary_mu
  })
  data <- do.call(rbind, args = data)
  data <- data[order(data$bin_size, data$fraction), ]
  data <- by(data, INDICES = list(data$chromosome, data$bin_size), FUN = identity)
  
  limits <- lapply(data, FUN = function(data) {
    is_num <- vapply(data, FUN = is.numeric, FUN.VALUE = FALSE)
    cbind(data[1L, !is_num, drop = FALSE], as.list(colMeans(data[, is_num, drop = FALSE], na.rm = TRUE)))
  })
  limits <- do.call(rbind, args = limits)
  rownames(limits) <- NULL
  limits <- limits[order(limits$chromosome, limits$bin_size), ]
  limits
}
