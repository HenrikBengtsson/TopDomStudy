library("TopDomStudy") ## https://github.com/HenrikBengtsson/TopDomStudy
library("progressr")   ## https://github.com/HenrikBengtsson/progressr
library("progress")
library("R.utils")

options(width = 140)
progressr::handlers("progress")

## NOTES: GSM2438426_ML4,human takes ~14 GiB RAM and ~25 mins to finish

## ---------------------------------------------------------------------
## Setup data
## ---------------------------------------------------------------------
known_samples <- c(
  "GSM2254215_ML1",
  "GSM2254216_ML2",
  "GSM2254217_ML3",  
  "GSM2254218_PL1",
  "GSM2254219_PL2",
  "GSM2438426_ML4"  ## Non-supported file format
)
known_organisms <- c("human", "mouse")

samples <- cmdArg(samples = c("GSM2254215_ML1", "GSM2254216_ML2", "GSM2254218_PL1", "GSM2254219_PL2"))
organisms <- cmdArg(organisms = "human")
stopifnot(all(samples %in% known_samples), all(organisms %in% known_organisms))


## ---------------------------------------------------------------------
## Process
## ---------------------------------------------------------------------
res <- TopDomStudy::compile_by_organism(
  samples = samples,
  organisms = organisms,
  path = file.path("hicData", "GSE84920"),
  path_dest = "compiledData"
)
print(res)

