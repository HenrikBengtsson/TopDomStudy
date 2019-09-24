options(width = 140)
progressr::handlers("progress")

## NOTES: 6 human + 2 mouse cell types takes ~40 minutes to complete
celltypes <- list(
  human = c("HAP1", "HeLa", "GM12878", "K562", "Asynchronous", "Nocadazole"),
  mouse = c("MEF", "Patski")
)

progressr::with_progress({
  res <- TopDomStudy::split_by_celltype(celltypes = celltypes, path = "compiledData")
})
print(res)
