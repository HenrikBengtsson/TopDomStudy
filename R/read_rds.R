#' Robustly Reads an RDS File
#'
#' @param pathname RDS file to read.
#'
#' @return The \R object read.
#'
#' @details
#' Uses [base::readRDS] internally but gives a more informative error message
#' on failure.
#'
#' @importFrom utils file_test
#' @export
#' @keywords internal
read_rds <- function(pathname) {
  if (!file_test("-f", pathname)) {
    stop(sprintf("No such file: %s", sQuote(pathname)))
  }
  tryCatch({
    readRDS(pathname)
  }, error = function(ex) {
    msg <- conditionMessage(ex)
    msg <- sprintf("readRDS() failed to read file %s (%.0f bytes). The reason was: %s",
                   sQuote(pathname), file.size(pathname), msg)
    ex$message <- msg
    stop(ex)
  })
}
