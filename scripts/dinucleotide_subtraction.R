rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

#this script will make a histogram of the most common dinucleotide cleavage sites
#for a sample type

sample_type = "mhv" #sample types: mhv, mm18S, mm28S