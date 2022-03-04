#! /bin/bash

# Author: Ajeet Singh
# Email: singh.ajeet@nipgr.ac.in

# This script add the gene_name or transcript_name in GTF file, if missing, by using respective ids

gtf=$1
#get the entries with transcript_id
grep -v "\bgene\b\|^#" $gtf|while read i
do
	if [[ ! "$i" =~ "transcript_name" ]] && [[ "$i" =~ "gene_name" ]];then
		transcript_id=$(echo "$i" | grep -o "transcript_id.*" |cut -f2 -d " ")
		echo "$i"|sed "s/transcript_id/transcript_id $transcript_id transcript_name/g"

	elif [[ ! "$i" =~ "gene_name" ]] && [[ "$i" =~ "transcript_name" ]];then
		gene_id=$(echo "$i" | grep -o "gene_id.*" |cut -f2 -d " ")
		echo "$i"|sed "s/gene_id/gene_id $gene_id gene_name/g"

        elif [[ ! "$i" =~ "gene_name" ]] && [[ ! "$i" =~ "transcript_name" ]];then
                gene_id=$(echo "$i" | grep -o "gene_id.*" |cut -f2 -d " ")
                transcript_id=$(echo "$i" | grep -o "transcript_id.*" |cut -f2 -d " ")
                echo "$i"|sed "s/gene_id/gene_id $gene_id gene_name/g;s/transcript_id/transcript_id $transcript_id transcript_name/g"
	else
		echo "$i"
	fi
done | awk -F\\t '!seen[$0]++' OFS='\t'
