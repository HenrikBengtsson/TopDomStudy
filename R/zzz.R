.onLoad <- function(libname, pkgname) {
  ## AD HOC: Only register vignette engine startup::selfonly during R CMD build
  if (nzchar(Sys.getenv("R_CMD"))) register_vignette_engine(pkgname)
}


