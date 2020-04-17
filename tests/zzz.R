## Post-testing cleanup
dirs <- c("figures", "overlapScoreData", "overlapScoreSummary", "topdomData")
for (dir in dirs) {
  if (!utils::file_test("-d", dir)) next
  unlink(dir, recursive = TRUE)
}

