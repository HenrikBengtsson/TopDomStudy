#' Read TopDom Data from Overlap Score Data
#'
#' @param pathname (character string) A \file{*,topdom.rds} file.
#'
#' @param format (character string) The returned format.
#'
#' @return A data frame with one TopDom domain per row.
#'
#' @importFrom utils file_test
#' @importFrom tibble as_tibble
#' @export
read_topdom_domains <- function(pathname, format = c("tibble", "data.frame")) {
  format <- match.arg(format)
  stopifnot(file_test("-f", pathname))

  seed <- gsub(".*,seed=([a-z0-9]+),.*", "\\1", pathname)
  stopifnot(nzchar(seed))
  seed <- eval(parse(text = sprintf("0x%s", seed)))
  stopifnot(is.numeric(seed), length(seed) == 1L, is.finite(seed))
  
  data <- readRDS(pathname)

  topdom <- data[["fraction=0.5"]]
  td <- topdom$domain
  config <- attributes(data)[c("bin_size", "min_cell_size", "window_size", "partition_by", "seed")]
  config$seed <- seed
  data <- cbind(td, config)
  stopifnot(is.data.frame(data), nrow(data) == nrow(td))
  
  if (format == "tibble") {
    data <- as_tibble(data)
  }
  
  data
}
