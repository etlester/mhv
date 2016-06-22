


rm(list=ls())

library(ggplot2)
library(readr)

sample_vector <- 1:8




mhv <- read_tsv("~/Documents/Hesselberth_Lab/mhv/data/28S_rRNA/mhv.sample1.mm28S.bed.gz", col_names = c('chrom','start','end', "ndupreads", "nreads"))

total_reads <- sum(mhv$nreads)
normalized_reads <- mhv$nreads/total_reads
ext_mhv <- cbind(mhv, normalized_reads)

colnames(ext_mhv) <- c('chrom','start','end', "ndupreads", "nreads","normreads")

ggplot(ext_mhv, aes(x = start, y = normreads)) + geom_line()