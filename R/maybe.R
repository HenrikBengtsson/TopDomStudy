#' An R Expression that Maybe Evaluated Later
#'
#' @param expr,substitute An \R expression (substituted by default).
#'
#' @param envir The environment where the expression should be evaluated
#' and from where global variables should be identified and frozen.
#'
#' @return A lazy [future::Future].
#'
#' @importFrom future sequential
#' @export
maybe <- function(expr, substitute = TRUE, envir = parent.frame) {
  if (substitute) expr <- substitute(expr)
  f <- sequential(expr, substitute = FALSE, envir = envir, lazy = TRUE)
  class(f) <- c("MaybeFuture", class(f))
}
