
rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

#alter this line to choose what sample type to plot
sample_type = "mhv" #sample types: mhv, mm18S, mm28S
lab_cutoff = .01

##y value cutoffs and nudges for labeling plots
if(sample_type == "mhv"){
  lab_cutoff <- .01
  xnudge <- -500
  ynudge <- .005
 } else if(sample_type == "mm18S"){
    lab_cutoff <- .0075
    xnudge <- 0
    ynudge <- .01
 } else {
  lab_cutoff <- .01
  xnudge <- 0
  ynudge <- 0
  }


#while loop to load all + strand bedgraph samples
setwd("~/Documents/Hesselberth_Lab/mhv/data")
nsamples = 8
i <- 1
while(i <= nsamples){
#create the object name
  o <-  paste0(sample_type, i)
#read the file and name the columns  
  f <- read_tsv(paste0(sample_type, "/mhv.sample", i,".", sample_type, ".+.bg"), 
                  col_names = c('condition','start','end', "nreads", "norm", "samplen"))

#change chrom name
  f$condition <- o
  f$samplen <- i
#normalize
  tot_reads <- sum(f$nreads)
  f$norm <- f$nreads/tot_reads
#save tot_reads to mhvUMI_reads*
  t <- paste0(sample_type, "UMI_reads",i)
  assign(t,tot_reads)
#sample num
  
#assign the file to the object name
  assign(o,f) 
 i <- i+ 1
}


#####need to figure out what Daphne is normalizing to
#she seems to be hard coding numbers from somewhere

#put all mhv files into a single data frame
if(sample_type == "mhv"){
  tdf <- rbind(mhv1, mhv2, mhv3, mhv4, mhv5, mhv6, mhv7, mhv8)
}
if(sample_type == "mm18S"){
  tdf <- rbind(mm18S1, mm18S2, mm18S3, mm18S4, mm18S5, mm18S6, mm18S7, mm18S8)
}
if(sample_type == "mm28S"){
  tdf <- rbind(mm28S1, mm28S2, mm28S3, mm28S4, mm28S5, mm28S6, mm28S7, mm28S8)
}





ggplot(tdf, aes(x = start, y = norm)) + geom_line() + 
  facet_wrap(~condition, ncol = 2) +
  labs(x = "Nucleotide", y = paste ('% cDNA reads aligned to', sample_type),
       title = paste("Cleavage sites in", sample_type, "rRNA")) +
      geom_text(data = subset(tdf, norm > lab_cutoff), 
                aes(start, norm, label = start, alpha = .7),
                nudge_x = xnudge, nudge_y = ynudge, angle = 45)



##plot hpi variable in two different color lines
##add tag for 9 hpi and 12hpi for color coding

tdf <- tbl_df(tdf)

tdf9 <- filter(tdf, samplen == 1 | samplen == 3 | samplen == 5 | samplen == 7 ) %>%
  mutate(hpi = 'hpi 9')
tdf12 <- filter(tdf, samplen == 2 | samplen == 4 | samplen == 6 | samplen == 8 ) %>%
  mutate(hpi = 'hpi 12')

tdf_hpi <- rbind(tdf9, tdf12)

##group sample types for plotting 9 ad 12 hours on top of one another
tdfa <- filter(tdf_hpi, samplen == 1 | samplen == 2) %>%
  mutate(sample_group = 'WT MHV | WT RNaseL')
tdfb <- filter(tdf_hpi, samplen == 3 | samplen == 4) %>%
  mutate(sample_group = 'WT MHV | RNaseL -/-')
tdfc <- filter(tdf_hpi, samplen == 5 | samplen == 6) %>%
  mutate(sample_group = 'NS2 Mut | WT RNase L')
tdfd <- filter(tdf_hpi, samplen == 7 | samplen == 8) %>%
  mutate(sample_group = 'NS2 Mut | RNaseL -/-')

tdf_hpi <- rbind(tdfa,tdfb,tdfc,tdfd)


ggplot(tdf_hpi, aes(x = start, y = norm, col = hpi)) + 
  geom_bar(position="dodge", stat = 'identity') + 
  facet_wrap(~sample_group, ncol = 1) +
  labs(x = 'nucleotide', y = paste('%cDNA reads aligned to', sample_type), 
       title = paste('Cleavage sites in', sample_type,"RNA" )) +
        geom_text(data = subset(tdf_hpi, norm > lab_cutoff), 
                  aes(start, norm, label = start, alpha = .7),
                  nudge_x = xnudge, nudge_y = ynudge, angle = 45)



