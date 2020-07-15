library(ggplot2)

## HiC read-pair data
dataset <- "human,HAP1"
chromosome <- 22
pathname <- system.file("compiledData", sprintf("%s,unique,chr=%s.rds", dataset, chromosome), package="TopDomStudy", mustWork=TRUE)
reads <- read_rds(pathname)

## Run TopDom on test and reference partitions and calcluate overlap scores
pathname <- overlap_scores_partitions(reads, dataset=dataset, partition_by="cells_by_half", min_cell_size=2L, nsamples=1L, seed=TRUE, mainseed=0xBEEF, reference_rho=1/2, rho=0.50, bin_size=10000, chrs="22")[[1]]

## TopDom overlap scores and lengths
tdos <- read_rds(pathname)$test[[1]]
tdos <- na.omit(tdos)

gg <- ggplot(tdos) + aes(x=best_length, y=best_score)
gg <- gg + labs(x="TAD length (kbp)", y="TAD overlap score")
gg <- gg + scale_y_continuous(limits=c(0,1), labels=scales::percent_format())
gg <- gg + scale_x_continuous(limits=c(0,NA), labels=scales::unit_format(unit="", scale=1e-3))
gg <- gg + geom_point() + geom_smooth(method="loess", formula=y~x)
print(gg)

## ggsave("overlap_scores_vs_length.pdf", gg, width=8, height=0.618*8)
