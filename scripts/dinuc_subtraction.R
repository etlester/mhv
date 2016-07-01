rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

sample_type = "mm28S" #sample types: mhv, mm18S, mm28S
freq_adj <- 1 #1 if you want freq adj, 0 if you don't

load(paste0('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/',
            sample_type,'/dinuc/', sample_type,'.dinuc_df.RData'))


#trim to only variables you want
if(freq_adj == 1){
df_small <- df %>%
  select(samplen, dinuc, freq_adj_norm_treads) #for selecting freq_adj_norm_treads
} else if(freq_adj == 0){
  df_small <- df %>%
  select(samplen, dinuc, norm_treads) #for selecting non freq_adj_norm_treads
}
names(df_small) <- c('samplen', 'dinuc', 'reads')


#make s1s3, s2s4,s3s5,s4s6,s5s7,s6s8
i <- 1
while(i <=6){
s1 <- i
s2 <- i +2
##pull out first subtraction pair
a <- filter(df_small, samplen == s1) 
#order dinucleotides in alphabetical order
  ordera_vector <- order(a$dinuc)
   ordered_a<- a[ordera_vector,] 
#pull out second subtraction pair
b <- filter(df_small, samplen == s2)
  orderb_vector <- order(b$dinuc)
  ordered_b<- b[orderb_vector,]
n <- paste0("s",s1,'s',s2)
diff <- ordered_a$reads - ordered_b$reads
d <- data.frame(ordered_a,ordered_b,diff,n)
names(d) <- c('sample1', 'dinuc1','reads1','sample2','dinuc2','reads2','diff','sub_pair')
assign(n,d)
i <- i +1
}

#make s8s6, s7s5, s6s4, s5s3, s4s2,s3s1
i <- 8
while(i >= 3){
  s1 <- i 
  s2 <- i-2
  ##pull out first subtraction pair
  a <- filter(df_small, samplen == s1) 
  #order dinucleotides in alphabetical order
  ordera_vector <- order(a$dinuc)
  ordered_a<- a[ordera_vector,] 
  #pull out second subtraction pair
  b <- filter(df_small, samplen == s2)
  orderb_vector <- order(b$dinuc)
  ordered_b<- b[orderb_vector,]
  n <- paste0("s",s1,'s',s2)
  diff <- ordered_a$reads - ordered_b$reads
  d <- data.frame(ordered_a,ordered_b,diff,n)
  names(d) <- c('sample1', 'dinuc1','reads1','sample2','dinuc2','reads2','diff','sub_pair')
  assign(n,d)
  i <- i -1
}

#select only the relevant subtraction pairs and put them into a df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)

#plot frequency adjusted data using geom_bar()
sub_df %>%
  #filter(samplen == 1) %>%
  ggplot(aes(x= dinuc1, y = diff)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~sub_pair, ncol = 2) +
  labs(title = paste(sample_type, 'Dinucleotide Subtraction \n dinucleotide frequancy adjusted'),
       x = 'Dinucleotide Pair',
       y = 'Difference in % cDNA Reads \n Normalized to Reference Dinucleotide Frequency') +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, vjust=0)) 

ggsave(paste0(
  '/Users/evanlester/Documents/Hesselberth_Lab/mhv/figures/dinuc_cleavage_subtraction/freq_adj/',
  sample_type, '.dinuc_sub_freq_adj.pdf'))



##############do the same thing with non-freq adjusted data

df_small <- df %>%
  select(samplen, dinuc, norm_treads) #for selecting non freq_adj_norm_treads

names(df_small) <- c('samplen', 'dinuc', 'reads')


#make s1s3, s2s4,s3s5,s4s6,s5s7,s6s8
i <- 1
while(i <=6){
  s1 <- i
  s2 <- i +2
  ##pull out first subtraction pair
  a <- filter(df_small, samplen == s1) 
  #order dinucleotides in alphabetical order
  ordera_vector <- order(a$dinuc)
  ordered_a<- a[ordera_vector,] 
  #pull out second subtraction pair
  b <- filter(df_small, samplen == s2)
  orderb_vector <- order(b$dinuc)
  ordered_b<- b[orderb_vector,]
  n <- paste0("s",s1,'s',s2)
  diff <- ordered_a$reads - ordered_b$reads
  d <- data.frame(ordered_a,ordered_b,diff,n)
  names(d) <- c('sample1', 'dinuc1','reads1','sample2','dinuc2','reads2','diff','sub_pair')
  assign(n,d)
  i <- i +1
}

#make s8s6, s7s5, s6s4, s5s3, s4s2,s3s1
i <- 8
while(i >= 3){
  s1 <- i 
  s2 <- i-2
  ##pull out first subtraction pair
  a <- filter(df_small, samplen == s1) 
  #order dinucleotides in alphabetical order
  ordera_vector <- order(a$dinuc)
  ordered_a<- a[ordera_vector,] 
  #pull out second subtraction pair
  b <- filter(df_small, samplen == s2)
  orderb_vector <- order(b$dinuc)
  ordered_b<- b[orderb_vector,]
  n <- paste0("s",s1,'s',s2)
  diff <- ordered_a$reads - ordered_b$reads
  d <- data.frame(ordered_a,ordered_b,diff,n)
  names(d) <- c('sample1', 'dinuc1','reads1','sample2','dinuc2','reads2','diff','sub_pair')
  assign(n,d)
  i <- i -1
}

#select only the relevant subtraction pairs and put them into a df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)

#plot frequency adjusted data using geom_bar()
sub_df %>%
  #filter(samplen == 1) %>%
  ggplot(aes(x= dinuc1, y = diff)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~sub_pair, ncol = 2) +
  labs(title = paste(sample_type, 'Dinucleotide Subtraction'),
       x = 'Dinucleotide Pair',
       y = 'Difference in % cDNA Reads') +
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=90, vjust=0)) 




#save non-freq adj plot
  ggsave(paste0(
    '/Users/evanlester/Documents/Hesselberth_Lab/mhv/figures/dinuc_cleavage_subtraction/non_freq_adj/',
    sample_type, '.dinuc_sub.pdf'))





