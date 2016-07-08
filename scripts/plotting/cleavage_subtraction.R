


rm(list=ls())

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

#alter this to choose alignment type
#sample types: mhv, mm18S, mm28S
alignment_type = "mhv" 
join_type <- "left" #join type can either be "left" or "inner"

if(alignment_type == 'mhv'){
  label_cutoff <- .003
  } else if(alignment_type == 'mm18S'){
    label_cutoff <- .003
  } else if(alignment_type == 'mm28S'){
    label_cutoff <- .003
  }



#collapsed section loads all + strand bedgraph samples and dinucleotides into an object
#named tdf
setwd("~/Documents/Hesselberth_Lab/mhv/data")
nsamples = 8
tdf <- list()
i <- 1
while(i <= nsamples){
  #create the object name
  o <-  paste0(alignment_type, i)
  #read the file and name the columns  
  f <- read_tsv(paste0(alignment_type, "/dinuc/dinuc.norm.mhv.sample", i,".", alignment_type, ".+.txt"), 
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
#select columns you care about
tdf_small <- select(tdf, samplen, end, norm, dinuc)


#make s1s3, s2s4,s3s5,s4s6,s5s7,s6s8
i <- 1
while(i <=6){
  s1 <- i
  s2 <- i +2
  ##pull out first subtraction pair
  a <- filter(tdf_small, samplen == s1) 
  #order dinucleotides in alphabetical order
  ordera_vector <- order(a$end)
  ordered_a<- a[ordera_vector,] 
  #pull out second subtraction pair
  b <- filter(tdf_small, samplen == s2)
  orderb_vector <- order(b$end)
  ordered_b<- b[orderb_vector,]
  n1 <- paste0("s",s1,'s',s2)
  if(join_type == "left") {
  d1 <- left_join(ordered_a, ordered_b, by = "end")
  }else if(join_type == "inner"){
  d1 <- inner_join(ordered_a, ordered_b, by = "end")
}
  d1 <- data.frame(d1,n1)
  names(d1) <- c('sample1', 'end','norm_reads1','dinuc1','sample2','norm_reads2','dinuc2','sub_pair')
  assign(n1,d1)
  #do the same thing such that you put b on the left and a on the right
  n2 <- paste0("s",s2,'s',s1)
  if(join_type == "left") {
    d2 <- left_join(ordered_b, ordered_a, by = "end")
  }else if(join_type == "inner"){
    d2 <- inner_join(ordered_b, ordered_a, by = "end")
  }
  d2 <- data.frame(d2,n2)
  names(d2) <- c('sample1', 'end','norm_reads1','dinuc1','sample2','norm_reads2','dinuc2','sub_pair')
  assign(n2,d2)
  i <- i +1
}


#select only the relevant subtraction pairs and put them into an object called sub_df
sub_df <- rbind(s1s3, s2s4, s3s1, s4s2, s5s7, s6s8, s7s5, s8s6)

#remove NAs from norm_reads2
sub_df[is.na(sub_df)] <- 0

#calculate subtraction differences
sub_df <- sub_df %>%
  mutate(diff = norm_reads1 - norm_reads2)

#set diff negative values to 0
index <- sub_df$diff < 0
sub_df$diff[index] <- 0



#look up table for labeling graphs
lut <- c("s1s3" = "s1s3 (WT MHV|WT RNaseL|9hpi - WT MHV|RNaseL -/-|9hpi)\n RNaseL activity inhibited by NS2",
         "s2s4" = "s2s4 (WT MHV|WT RNaseL|12hpi - WT MHV|RNaseL -/-|12hpi)\n RNaseL activity inhibited by NS2",
         "s3s1" = "s3s1 (WT MHV|RNaseL -/-|9hpi - WT MHV|WT RNaseL|9hpi)\n not RNaseL activity with NS2 present",
         "s4s2" = "s4s2 (WT MHV|RNaseL -/-|12hpi - WT MHV|WT RNaseL|12hpi)\n not RNaseL activity with NS2 present",
         "s5s7" = "s5s7 (NS2 Mut|WT RNaseL|9hpi - NS2 Mut|RNaseL -/-|9hpi)\n RNaseL activity",
         "s6s8" = "s6s8 (NS2 Mut|WT RNaseL|9hpi - NS2 Mut|RNaseL -/-|12hpi)\n RNaseL activity",
         "s7s5" = "s7s5 (NS2 Mut|RNaseL -/-|9hpi - NS2 Mut|WT RNaseL|9hpi)\n not RNaseL activity",
         "s8s6" = "s8s6 (NS2 Mut|RNaseL -/-|12hpi - NS2 Mut|WT RNaseL|9hpi)\n not RNaseL activity")

sub_df$sub_pair <- lut[sub_df$sub_pair]


xnudge <- 0
ynudge <- .0025
t <- theme(text = element_text(size=10),
           axis.text.x = element_text(angle=90, vjust=0),
           strip.text.x = element_text(size = 9))


#plot subtraction pairs
sub_df %>%
  ggplot(aes(x= end, y = diff)) + 
  geom_line(alpha = .5) +
  facet_wrap(~sub_pair, ncol = 2, scales = "free") +
  labs(title = paste(alignment_type, 'Subtractive Analysis',join_type,'join'),
       x = 'Nuceotide',
       y = 'Difference in % cDNA Reads') +
  t +
  geom_text_repel(data = subset(sub_df, diff > label_cutoff), 
            aes(end, diff, label = paste0(dinuc1, end)), size = 3)


#code fragment that can put in superscripts
##(bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')'))
#link: http://stackoverflow.com/questions/27445452/superscript-and-subscript-axis-labels-in-ggplot2


#save subtraction pair plot
ggsave(paste0(
'/Users/evanlester/Documents/Hesselberth_Lab/mhv/figures/cleavage_subtraction/', join_type,'_join/',
  alignment_type,'_',join_type, '_sub.pdf'), width = 10.5, height = 8, units = "in")


###save sub_df and tdf as a references for the Rmarkdown document
save(sub_df,file = paste0('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/',
                      alignment_type,'/', alignment_type,'.cleavage_subtraction.RData'))

lut <- c('1' = "WT MHV | WT RNaseL | 9hpi",
         '2' = "WT MHV | WT RNaseL | 12hpi",
         '3' = "WT MHV | RNaseL -/- | 9hpi",
         '4' = "WT MHV | RNaseL -/- | 12hpi",
         '5' = "NS2 Mut | WT RNaseL | 9hpi",
         '6' = "NS2 Mut | WT RNaseL | 12hpi",
         '7' = "NS2 Mut | RNaseL -/- | 9hpi",
         '8' = "NS2 Mut | RNaseL -/- | 12hpi")

tdf$samplen <- lut[tdf$samplen]

save(tdf,file = paste0('/Users/evanlester/Documents/Hesselberth_Lab/mhv/data/',
                      alignment_type,'/', alignment_type,'.cleavage_location.RData'))






