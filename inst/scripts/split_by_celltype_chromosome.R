options(width = 140)
progressr::handlers("progress")

## NOTES: 6 human + 2 mouse cell types takes ~40 minutes to complete
known_celltypes <- list(
  human = c("HAP1", "HeLa", "GM12878", "K562", "Asynchronous", "Nocadazole"),
  mouse = c("MEF", "Patski")
)
known_chromosomes <- 1:23

celltypes <- list(
  human = c("HAP1")
)
chromosomes <- c(12, 16, 22)

progressr::with_progress({
  res <- TopDomStudy::split_by_celltype_chromosome(celltypes = celltypes, chromosomes = chromosomes, path = "compiledData")
})
print(res)
