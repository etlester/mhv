

rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)





#while loop to load all + strand bedgraph samples
setwd("~/Documents/Hesselberth_Lab/mhv/data")
sample_type1 <- "mhv" #sample types: mhv, mm18S, mm28S
sample_type2 <- "mm18S"
sample_type3 <- 'mm28S'


treads_df <- c()

#load mhv
sample_type <- sample_type1
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
  #put total reads in vector
  sample_group <- sample_type
  treads_entry <- c(sample_group, i, tot_reads)
  treads_df <- rbind(treads_df, treads_entry)
  
  
  
  #assign the file to the object name
  assign(o,f) 
  i <- i+ 1
}

#load mm18S
sample_type <- sample_type2
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
  #put total reads in vector
  sample_group <- sample_type
  treads_entry <- c(sample_group, i, tot_reads)
  treads_df <- rbind(treads_df, treads_entry)
  
  #assign the file to the object name
  assign(o,f) 
  i <- i+ 1
}


#load mm28S
sample_type <- sample_type3
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
  #put total reads in vector
  sample_group <- sample_type
  treads_entry <- c(sample_group, i, tot_reads)
  treads_df <- rbind(treads_df, treads_entry)
  
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

Daphnes_mhv_data <- c(19416, 8470, 14220, 16826, 7188, 7418, 7051, 8075)
Daphnes_mm18S_data <- c(47629, 20971,31270,35293,22685,23222,13481,17885)
Daphnes_mm28S_data <- c(224478,90231,132454,152496,110639,116313,73160,84090)



Daphnes_data <- c(Daphnes_mhv_data,Daphnes_mm18S_data,Daphnes_mm28S_data)
treads_df <- cbind(treads_df, Daphnes_data)

colnames(treads_df) <- c('sample_group','sample', 'total_reads', 'daphnes_reads')

treads_df <- tbl_df(treads_df)
treads_df <- mutate(treads_df, total_reads = as.integer(total_reads),
                        daphnes_reads = as.integer(daphnes_reads))

#want to plot the difference between my data and daphne's in terms of aligned reads

treads_df <- mutate(treads_df, diff = total_reads - daphnes_reads)



#plot of total number of my aligned reads per sample type
ggplot(treads_df, aes(x = sample_group, y = total_reads, fill = sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = "Sample Group", y = paste ('Number of aligned reads'),
       title = "Aligned reads per sample type \n (my data)")
       

#plot of total number of Daphne's aligned reads per sample type
ggplot(treads_df, aes(x = sample_group, y = daphnes_reads, fill = sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = "Sample Group", y = paste ('Number of aligned reads'),
       title = "Aligned reads per sample type \n (Daphne's data)")

##normalize difference data by average number of reads per sample group
mhvnorm <- filter(treads_df, sample_group == 'mhv') %>%
  mutate(mean_num_reads = mean(total_reads + daphnes_reads)/2) %>%
  mutate(normalized_diff = diff/mean_num_reads)
mm18Snorm <- filter(treads_df, sample_group == 'mm18S') %>%
  mutate(mean_num_reads = mean(total_reads + daphnes_reads)/2) %>%
  mutate(normalized_diff = diff/mean_num_reads)
mm28Snorm  <- filter(treads_df, sample_group == 'mm28S') %>%
  mutate(mean_num_reads = mean(total_reads + daphnes_reads)/2) %>%
  mutate(normalized_diff = diff/mean_num_reads)

#combine back into a big dataset
treads_df_norm <- rbind(mhvnorm, mm18Snorm, mm28Snorm)

      
##plot of difference between my reads and daphne's reads
ggplot(treads_df_norm, aes(x = sample_group, y = normalized_diff, fill = sample)) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(x = "Sample Group", y = paste ('Number of aligned reads / \n mean(# reads per sample group)'),
       title = "Difference (my data - Daphne's) 
       \n in aligned reads per sample type")



