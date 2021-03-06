% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maybe.R
\name{maybe}
\alias{maybe}
\title{An R Expression that Maybe Evaluated Later}
\usage{
maybe(expr, substitute = TRUE, envir = parent.frame())
}
\arguments{
\item{expr, substitute}{An \R expression (substituted by default).}

\item{envir}{The environment where the expression should be evaluated
and from where global variables should be identified and frozen.}
}
\value{
A lazy \link[future:Future-class]{future::Future}.
}
\description{
An R Expression that Maybe Evaluated Later
}
\keyword{internal}
