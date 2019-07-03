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


#' Asserts the Truth of R Expressions
#'
#' @param \dots Zero or more \R expressions to be asserted to be TRUE.
#'
#' @return Nothing.
#'
#' @details
#' A bare bone, faster version of [base::stopifnot].
#'
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
  invisible()
}


#' Tweaks a Function to be Evaluted using try()
#'
#' @param fcn A function.
#'
#' @return A function that calls `try(fcn(...), silent = TRUE)`.
#'
#' @details
#' The function  [base::stopifnot].
#'
#' @export
## https://github.com/HenrikBengtsson/r-ideas/issues/75
Try <- function(fcn) function(...) try(fcn(...), silent = TRUE)


#' CRC32 Checksum of an Object
#'
#' @param x An \R object.
#'
#' @return An eight-character string.
#'
#' @export
#' @importFrom digest digest
#' @export
crc32 <- function(x) digest(x, algo = "crc32")


mprintf <- function(..., appendLF = TRUE) message(sprintf(...), appendLF = appendLF)

#' @importFrom utils capture.output
mprint <- function(..., appendLF = TRUE) {
  message(paste(capture.output(print(...)), collapse = "\n"), appendLF = appendLF)
}

#' @importFrom utils capture.output str
mstr <- function(..., appendLF = TRUE) {
  message(paste(capture.output(str(...)), collapse = "\n"), appendLF = appendLF)
}


## WORKAROUND: This will create a dummy progressor() until
## the 'progressr' package is publicly available / installed.
import_progressor <- function() {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor
  } else {
    progressor <- function(...) {
      function(...) invisible()
    }
  }
}
