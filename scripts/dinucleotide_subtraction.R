rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

#this script will plot bar charts of dinucleotide cleavge subtraction pairs

sample_type = "mhv" #sample types: mhv, mm18S, mm28S

load(paste0('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/',
              sample_type,'/dinuc/', sample_type,'.dinuc_df.RData'))

#trim to only variables you want
df_small <- df %>%
  select(samplen, dinuc, freq_adj_norm_treads)
names(df_small) <- c('sample', 'dinuc', 'reads')


#generate 8 subtraction pairs: s1s3,s2s4,s5s7,s6s8,s8s6,s7s5,s4s2,s3s1

#make s1s3, s2s4,s3s5,s4s6,s5s7,s6s8
i <- 1
while(i <=6){
  samplea <- i
  sampleb <- i+2
  a <- filter_(df_small, "sample" == "samplea")
  b <- filter_(df_small, "sample" == "sampleb")
  n <- paste0("s",samplea,'s',sampleb)
  d <- left_join(a,b, by = 'dinuc')
  e <- mutate_(d, diff = "reads.x" - "reads.y")
  assign(n,e)
  i <- i+1
}   



#make s8s6, s7s5, s6s4, s5s3, s4s2,s3s1
i <- 8
while(i >= 3){
  samplea <- i
  sampleb <- i-2
  a <- filter(df_small, sample == samplea)
  b <- filter(df_small, sample == sampleb)
  n <- paste0("s",samplea,'s',sampleb)
  d <- inner_join(a,b, by = 'dinuc')
  d <- d %>% mutate(diff = freq_adj_norm_treads.x - freq_adj_norm_treads.y)
  assign(n,d)
  i <- i-1
}  

#select only the relevant subtraction pairs and put them into a df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)





