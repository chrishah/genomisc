#!/bin/bash

file=$1
cutoff=$2
id=$3
num=$4

prefix=$(echo $file | sed 's/.txt$//g')

if [ $# -ne 4 ]; then
	echo -e "\nusage: ./plot_pope.sh <PATH:infile.txt> <cutoff:float> <string:ID> <INT:number_of_replicates> \n"
	exit
fi

echo $file
echo $cutoff
echo $id
echo $prefix
echo $num

echo "#R script
setwd('$(pwd)')
data <- read.delim(file = '$file', header = T, sep = '\\\t')
rel_rank_cutoff <- $cutoff
top <- subset(data, avg_rank_rel >= rel_rank_cutoff, select = c(chrom,bp,SNPID,std_rank,avg_rank_rel))

svg(filename = '$prefix-top-$cutoff.svg')
plot(data\$avg_rank_rel, data\$std_rank/max(data\$std_rank), pch=19, cex=0.2, ylab = 'relative rank standard deviation', xlab = 'average relative rank', main = 'Transformed average Bayesfactor ranks ($id) across $num runs', xlim=c(0,1), ylim=c(0,1))
points(top\$avg_rank_rel, top\$std_rank/max(data\$std_rank), col = 'red', pch=19, cex=0.5)

abline(v=rel_rank_cutoff, lty = 3, col='red')
dev.off()

write.table(top, file='$prefix-top-$cutoff.tsv', sep = '\\\t', col.names = T, quote = F, row.names = F)" > plot.R

Rscript plot.R
mv plot.R $prefix.R
