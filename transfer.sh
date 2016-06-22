#script for downloading data from tesla to my computer

#im changing things to see it git notices
#changing some more things

#copy mhv data
scp \
lestere@amc-tesla:~/mhv/results/20160620/1.mhv/mhv.sample*.mhv.bed.gz \
~/Documents/Hesselberth_Lab/mhv/data/mhv

#copy 28S data
scp \
lestere@amc-tesla:~/mhv/results/20160620/2.28S_rRNA/mhv.sample*.mm28S.bed.gz \
 ~/Documents/Hesselberth_Lab/mhv/data/28S_rRNA

#copy 18S data
scp \
lestere@amc-tesla:~/mhv/results/20160620/3.18S_rRNA/mhv.sample*.mm18S.bed.gz \
 ~/Documents/Hesselberth_Lab/mhv/data/18S_rRNA




