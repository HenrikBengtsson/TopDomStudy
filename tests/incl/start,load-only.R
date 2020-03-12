## Record original state
ovars <- ls()
oenvs <- oenvs0 <- Sys.getenv()
oopts0 <- options()

## Default options for tests
oopts <- options()

## Default test set
dataset <- "human,HAP1"
chromosomes <- getOption("TopDomStudy.tests.chromosomes", "22")
bin_sizes <- getOption("TopDomStudy.tests.bin_sizes", c(50e3, 100e3, 200e3))
rhos <- getOption("TopDomStudy.tests.rhos", c(0.10, 0.20, 0.30, 0.40, 0.50))
nsamples <- getOption("TopDomStudy.tests.nsamples", 5L)

