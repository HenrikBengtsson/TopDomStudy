#' @export
read_rds <- function(pathname) {
  tryCatch({
    readRDS(pathname)
  }, error = function(ex) {
    msg <- conditionMessage(ex)
    msg <- sprintf("readRDS() failed to read file %s (%.0f bytes). The reason was: %s",
                   sQuote(pathname), file.size(pathname), msg)
    ex$message <- msg
    signalCondition(ex)
  })
}
