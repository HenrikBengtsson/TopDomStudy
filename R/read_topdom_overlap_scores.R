#' Read TopDom Regions from TopDomStudy Files
#'
#' @param pathname (character string) A \file{*.rds} file.
#'
#' @param format (character string) The returned format.
#'
#' @return A data frame with TopDom domains with one overlap score and
#' one length per domain.
#'
#' @importFrom utils file_test
#' @importFrom tibble as_tibble
#' @export
read_topdom_overlap_scores <- function(pathname, format = c("tibble", "data.frame")) {
  format <- match.arg(format)
  stopifnot(file_test("-f", pathname))

  fraction <- gsub(".*,test=([0-9.]+),.*", "\\1", basename(pathname))
  stopifnot(nzchar(fraction))
  fraction <- as.numeric(fraction)
  stopifnot(is.numeric(fraction), length(fraction) == 1L, is.finite(fraction),
            fraction > 0, fraction <= 1/2)
  
  seed <- gsub(".*,seed=([a-z0-9]+).*", "\\1", basename(pathname))
  stopifnot(nzchar(seed))
  seed <- eval(parse(text = sprintf("0x%s", seed)))
  stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed))
  
  data <- read_rds(pathname)
  
  config <- attributes(data)[c("chromosome", "bin_size", "min_cell_size", "window_size", "partition_by", "seed")]
  config$seed <- seed
  config$fraction <- fraction
  config <- config[c("chromosome", "bin_size", "fraction", "min_cell_size", "window_size", "partition_by", "seed")]

  ## TODO: Drop support for '^fraction=..." /HB 2020-03-10
  if (any(grepl("^fraction=0.5$", names(data)))) {
    .Deprecated(msg = "The 'fraction=0.5' name is deprecated; please use 'test=0.5' instead")
  }
  test_partition <- grep("^(fraction|test)=0.5$", names(data))
  stopifnot(length(test_partition) == 1L)
  topdom <- data[[test_partition]]

  if (inherits(topdom, "try-error")) {
    td <- data.frame(best_score = double(0L), best_length = integer(0L))
    config <- as.data.frame(config)
    config <- config[integer(0L), ]
  } else {
    td <- topdom[[1]][,1:2]
  }
  data <- cbind(chr = config[[1L]], td, config[-1L])
  stopifnot(is.data.frame(data), nrow(data) == nrow(td))
  
  if (format == "tibble") {
    data <- as_tibble(data)
  }
  
  data
}
