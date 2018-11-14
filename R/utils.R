#' Assigns a value to a variable, unless it already exists
#'
#' @param x (RHS) The variables to be assigned.
#'
#' @param value (LHS) The expression whose value to assign to the variables.
#'
#' @return (invisibly) the value of the variable.
#'
#' @export
`%<-?%` <- function(x, value) {
  target <- substitute(x)
  expr <- substitute(value)

  target_name <- as.character(target)
  envir <- parent.frame(1L)
  if (exists(target_name, envir = envir, inherits = FALSE)) {
    value <- get(target_name, envir = envir, inherits = FALSE)
    return(invisible(value))
  }
  value <- eval(expr, envir = envir)
  assign(target_name, value = value, envir = envir, inherits = FALSE,
         immediate = TRUE)
  invisible(value)
}

## See https://github.com/HenrikBengtsson/r-ideas/issues/74


## Faster than stopifnot()
#' @export
stop_if_not <- function(...) {
  res <- list(...)
  for (ii in seq_along(res)) {
    res_ii <- .subset2(res, ii)
    if (length(res_ii) != 1L || is.na(res_ii) || !res_ii) {
        mc <- match.call()
        call <- deparse(mc[[ii + 1]], width.cutoff = 60L)
        if (length(call) > 1L) call <- paste(call[1L], "....")
        stop(sprintf("%s is not TRUE", sQuote(call)),
             call. = FALSE, domain = NA)
    }
  }
  
  NULL
}


## https://github.com/HenrikBengtsson/r-ideas/issues/75
#' @export
Try <- function(fcn) function(...) try(fcn(...), silent = TRUE)


## https://github.com/eddelbuettel/digest/issues/84
#' @importFrom digest digest
#' @export
crc32 <- local({
  digest <- digest::digest
  
  function(x) {
    s <- digest(x, algo = "crc32")
    npad <- 8L - nchar(s)
    if (npad > 0L) s <- paste(c(rep("0", npad), s), collapse = "")
    s
  }
})


	