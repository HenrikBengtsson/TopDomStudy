# Notes

The files within here, were created as follows.

Extract HAP1 data for chromosomes 12, 16, and 22:

```sh
dir.create("compiledData")
data <- readRDS("compiledData/human,HAP1,unique.rds")
for (chr in c(12, 16, 22)) {
  data_chr <- subset(data, chr_a == chr & chr_b == chr_a)
  saveRDS(data_chr, file.path("compiledData", sprintf("human,HAP1,unique,chr=%d.rds", chr)))
}
```

```sh
mkdir -p inst/hicData
src=hicData/GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt.gz
dest=inst/${src/.txt.gz/,chr1_vs_chr1,lines=1-10000}.txt
zcat ../$src | head -20000 | grep -P '\t+chr1\t.*\t+chr1\t' | head -10000 > $dest
```
