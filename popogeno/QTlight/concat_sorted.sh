#!/bin/bash

for rep in {1..10}
do
	echo -e "\nProcessing replicate $rep:\n"
	out=$(($rep+10))
	echo -e "*.xtx files"
	for file in $(ls -1 replicate_$rep/Diplo-1M-replicate-$rep-* | grep "xtx$"); do cat "$file" | perl -ne 'chomp; @c=split("/"); @a=split("\t", $c[-1]); @b=split("-", $a[0]); $b[1] =~ s/.txt//; $temp=sprintf("%s-%07d.txt",$b[0],$b[1]); $a[0]=$temp; $out=join("\t", @a); print "$out\n"'; done | sort -n > rep_$out.xtx

	echo -e "*.bf files"
	for file in $(ls -1 replicate_$rep/Diplo-1M-replicate-$rep-* | grep "bf$"); do cat "$file" | perl -ne 'chomp; @c=split("/"); @a=split("\t", $c[-1]); @b=split("-", $a[0]); $b[1] =~ s/.txt//; $temp=sprintf("%s-%07d.txt",$b[0],$b[1]); $a[0]=$temp; $out=join("\t", @a); print "$out\n"'; done | sort -n > rep_$out.bf

	count=$(cat rep_$out.bf | wc -l)
	echo -e "counted $count SNPs\n"
done
