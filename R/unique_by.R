#' @param x A [base:data.frame]
#'
#' @param by The column, as an integer or by name, to be used
#' to split up `x` in chunks.
#' If `NA` (default), then the first factor, string, or
#' integer column will be used.  If neither type exists, then
#' the first column will be used.
#'
#' @param \ldots Additional arguments passed to [base::unique]
#' per chunk.
#'
#' @return A [base:data.frame] with duplicated rows drop.
#'
#' @details
#' `unique_by(x)` is typically more memory efficient than
#' `unique(x)`.
#'
#' @export
unique_by <- function(x, by = NA_integer_, ...) {
  stopifnot(is.data.frame(x))
  if (nrow(x) <= 1L) return(x)
  if (ncol(x) == 0L) return(x)
  rownames_org <- rownames(x)
  stopifnot(length(by) == 1L)
  if (is.na(by)) {
    by <- which(vapply(x, FUN = is.factor, FUN.VALUE = FALSE))[1]
    if (is.na(by)) by <- which(vapply(x, FUN = is.character, FUN.VALUE = FALSE))[1]
    if (is.na(by)) by <- which(vapply(x, FUN = is.integer, FUN.VALUE = FALSE))[1]
    if (is.na(by)) by <- 1L
  } else if (is.character(by)) {
    stopifnot(by %in% colnames(x))
  } else {
    stopifnot(1L <= by, by <= ncol(x))
  }
  
  by_values <- unique(x[[by]])
  if (length(by_values) == nrow(x)) return(x)
  
  ux <- vector("list", length = length(by_values))
  for (ii in seq_along(by_values)) {
    rows <- which(x[[by]] == by_values[ii])
    x_ii <- x[rows, , drop = FALSE]
    x <- x[-rows, , drop = FALSE]
    ux[[ii]] <- unique(x_ii, ...)
    x_ii <- NULL
  }
  x <- x_by <- rows <- NULL
  x <- do.call(rbind, args = ux)
  ux <- NULL
  rows <- match(rownames_org, table = rownames(x))
  rows <- rows[!is.na(rows)]
  x <- x[rows, , drop = FALSE]
  x
}
