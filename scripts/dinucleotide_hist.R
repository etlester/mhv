rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

#this script will make a histogram of the most common dinucleotide cleavage sites
#for a sample type

sample_type = "mhv" #sample types: mhv, mm18S, mm28S

#while loop to load all + strand dinucleotide samples
setwd("~/Documents/Hesselberth_Lab/mhv/data")
nsamples = 8
tdf <- list()
i <- 1
while(i <= nsamples){
  #create the object name
  o <-  paste0(sample_type, i)
  #read the file and name the columns  
  f <- read_tsv(paste0(sample_type, "/dinuc/dinuc.norm.mhv.sample", i,".", sample_type, ".+.txt"), 
                col_names = c('condition','start','end', "nreads", "norm", "dinuc"))
  
  #change chrom name
  f$condition <- o
  f$samplen <- i
  #convert to tbl_df
  f < tbl_df(f)
  
  #add sample number as column
  f <- mutate(f, samplen = i)
  
  
#add NS2 condition column (1 = WT, 0 = Mutant)
  if(f$samplen == 1 | f$samplen == 2 |f$samplen == 3 |f$samplen == 4){
    f <- mutate(f, ns2 = 1)
  } else{
    f <- mutate(f, ns2 = 0)
  }
  
#add RNaseL condition column (1 = WT, 0 = -/-)
  if(f$samplen == 1 | f$samplen == 2 |f$samplen == 5 |f$samplen == 6){
    f <- mutate(f, rnasel = 1)
  } else{
    f <- mutate(f, rnasel = 0)
  }
  
#add HPI as column
  if(f$samplen %% 2 == 0 ){
    f<- mutate(f, hpi = 12)
  } else {
    f <- mutate(f, hpi = 9)
  }
  
  #assign the file to the object name
  assign(o,f) 
  
  
  #add to big data frame
  tdf <- bind_rows(tdf, f)
  
  i <- i+ 1
}
#reorder columns
tdf <- tdf[c("condition", "samplen", "start",'end', 'nreads','norm','dinuc','ns2','rnasel','hpi')]
#make sure that tdf is a tbl_df 
tdf <- tbl_df(tdf)

####the above block will generate and organize all the data needed for analysis###


##load dinucleotide frequencies from the FASTA reference files

#normalize by dinucleotide frequency in the original FASTA
if(sample_type == 'mhv'){
  dinuc_freq <- read_tsv('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/mhv/mhv_dinuc_freq.txt')
} else if(sample_type == 'mm18S'){
  dinuc_freq <- read_tsv('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/mm18S/mm18S_dinuc_freq.txt')
} else {
  dinuc_freq <- read_tsv('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/mm28S/mm28S_dinuc_freq.txt')
}
dinuc_freq <- select(dinuc_freq, nuc, count, freq)
names(dinuc_freq) <- c('dinuc','count', 'freq')

#create sorted list of dinucleotides
dinuc <- tdf %>% 
  #select(samplen, dinuc, nreads) %>%
  group_by(samplen, dinuc) %>%
  summarise(total_reads = sum(nreads)) %>%
  arrange(samplen, desc(total_reads))

#left join dinucleotide table with frequency table
dinuc_big <- left_join(dinuc, dinuc_freq, by = 'dinuc')
#add column of tota reads for each dinuc * the frequency of each dinuc 
#from the reference file. this column is call freq_norm_reads
dinuc_big <- mutate(dinuc_big, freq_norm_reads = total_reads * freq)

#sum the number of freq_norm_reads and divide each dinucleotide by that to get
#the normalized value, which will be plotted
final_dinuc <- dinuc_big %>%
  group_by(samplen) %>%
  mutate(norm_treads = freq_norm_reads/sum(freq_norm_reads))

#plot using geom_bar()
ggplot(final_dinuc, aes(x= dinuc, y = norm_treads)) + 
geom_bar(stat = "identity") +
facet_wrap(~samplen, ncol = 2) +
  labs(title = paste(sample_type, 'Dinucleotide Cleavage'),
       x = 'Dinucleotide',
       y = 'Cleavage frequency \n normalized to reference dinucleotide frequency')










