


rm(list=ls())

library(ggplot2)
library(readr)

#alter this to choose alignment type
#sample types: mhv, mm18S, mm28S
alignment_type = "mm28S" 




#while loop to load all + strand bedgraph samples
setwd("~/Documents/Hesselberth_Lab/mhv/data")
nsamples = 8
i <- 1
while(i <= nsamples){
  #create the object name
  condition <- paste0(alignment_type, i)
  o <-  paste0("sample", i)
  #read the file and name the columns  
  f <- read_tsv(paste0(alignment_type, "/mhv.sample", i,".", alignment_type, ".+.bg"), 
                col_names = c('condition','start','end', "nreads", "norm", "samplen"))
  
  #change condition name
  f$condition <- condition
  f$samplen <- i
  #normalize
  tot_reads <- sum(f$nreads)
  f$norm <- f$nreads/tot_reads
  #save tot_reads to mhvUMI_reads*
  t <- paste0(alignment_type, "UMI_reads",i)
  assign(t,tot_reads)
  #sample num
  
  #assign the file to the object name
  assign(o,f) 
  i <- i+ 1
}

#make subtraction plots
#options
# 1. bin data into bins of 10
# 2. inner_join data such that only rows in both sets are included
# 3. full_join data such that all values, all rows are included





#put all mhv files into a single data frame
tdf <- rbind(sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8)



tdf_small <- select(tdf, start, norm, samplen)

#make s1s3, s2s4,s3s5,s4s6,s5s7,s6s8
i <- 1
while(i <=6){
samplea <- i
sampleb <- i+2
a <- filter(tdf_small, samplen == samplea)
b <- filter(tdf_small, samplen == sampleb)
n <- paste0("s",samplea,'s',sampleb)
df <- inner_join(a,b, by = 'start')
df <- mutate(df, diff = norm.x - norm.y, subtraction = n)
assign(n,df)
i <- i+1
}   


#make s8s6, s7s5, s6s4, s5s3, s4s2,s3s1
i <- 8
while(i >= 3){
  samplea <- i
  sampleb <- i-2
  a <- filter(tdf_small, samplen == samplea)
  b <- filter(tdf_small, samplen == sampleb)
  n <- paste0("s",samplea,'s',sampleb)
  df <- inner_join(a,b, by = 'start')
  df <- mutate(df, diff = norm.x - norm.y, subtraction = n)
  assign(n,df)
  i <- i-1
}  



#select only the relevant subtraction pairs and put them into a df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)

#plot the subtraction pairs
ggplot(sub_df, aes(x = start, y = diff)) + geom_line() + 
  facet_wrap(~subtraction, ncol = 2) +
  labs(x = "Nucleotide", y = paste ('difference in % cDNA reads \n aligned to', alignment_type),
     title = paste("Subtractive analysis of", alignment_type))

#this is eventually how I want data plotted 
#data frame needs hpi column and set condition types
#ggplot(tdf, aes(x = start, y = norm, fill = hpi)) + geom_bar(position="dodge") +facet_wrap(~condition, ncol = 1)



