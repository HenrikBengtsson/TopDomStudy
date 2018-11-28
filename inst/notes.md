The files within here, was created using:

```sh
data <- readRDS("compiledData/human,HAP1,unique.rds")
data22 <- subset(data, chr_a == 22, chr_b = chr_a)
saveRDS(data22, "human,HAP1,unique,chr=22.rds")
```

```sh
mkdir -p inst/hicData
src=hicData/GSE35156_GSM862720_J1_mESC_HindIII_ori_HiC.nodup.hic.summary.txt.gz
dest=inst/${src/.txt.gz/,chr1_vs_chr1,lines=1-10000}.txt
zcat ../$src | head -20000 | grep -P '\t+chr1\t.*\t+chr1\t' | head -10000 > $dest
```
