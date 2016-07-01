rm(list=ls())

library(ggplot2)
library(readr)
library(dplyr)

sample_type = "mm28S" #sample types: mhv, mm18S, mm28S
keep_negatives <- 0 #1 = keep negative values 
                    #0 = throw out negative values and set them to

load(paste0('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/',
            sample_type,'/dinuc/', sample_type,'.dinuc_df.RData'))


#trim to only the variables you want
df_small <- select(df, samplen, dinuc, freq_adj_norm_treads, norm_treads)
names(df_small) <- c('samplen', 'dinuc', 'freq_adj_reads', 'reads')


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
n1 <- paste0("s",s1,'s',s2)
d1 <- data.frame(ordered_a,ordered_b,n1)
names(d1) <- c('sample1', 'dinuc1','freq_adj_reads1','reads1', 'sample2','dinuc2','freq_adj_reads2','reads2','sub_pair')
assign(n1,d1)
#do the same thing such that you put b on the left and a on the right
n2 <- paste0("s",s2,'s',s1)
d2 <- data.frame(ordered_b,ordered_a,n2)
names(d2) <- c('sample1', 'dinuc1','freq_adj_reads1','reads1', 'sample2','dinuc2','freq_adj_reads2','reads2','sub_pair')
assign(n2,d2)
i <- i +1
}

#select only the relevant subtraction pairs and put them into a df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)
sub_df_adj <- mutate(sub_df, diff_adj = freq_adj_reads1 - freq_adj_reads2)
sub_df <- mutate(sub_df, diff = reads1 - reads2)

#throw out negative values to compare to Daphne's plots
if(keep_negatives == 0){
  for(i in 1:nrow(sub_df)){
    if(sub_df[i,'diff'] < 0){
      sub_df[i,'diff'] <- 0
      }
    }
  for(i in 1:nrow(sub_df_adj)){
    if(sub_df_adj[i,'diff_adj'] < 0){
      sub_df_adj[i,'diff_adj'] <- 0
    }
  }
  }


#plot frequency adjusted data using geom_bar()
sub_df_adj %>%
  ggplot(aes(x= dinuc1, y = diff_adj)) + 
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





