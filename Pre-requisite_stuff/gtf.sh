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
done | awk -F\\t '!seen[$0]++' OFS='\t' > genome_transcript.gtf

# This part is for gtf exclusively for FusionBloom

#$2=Give cdna file as input
if [[ ! -z "$2" ]];then
	grep -o "transcript_id.*" genome_transcript.gtf|cut -f2 -d " "|sort -u|sed 's/"\|";//g'|sort >tid_gtf
	grep ">" $2|cut -f1 -d " "|sed 's/>//g' |sort >tid_cdna
	comm -12 tid_cdna tid_gtf > comm_tid
	grep -wf comm_tid genome_transcript.gtf > comm.gtf
	sort -t $'\t' -k1,1V -k4,4n -k5,5n comm.gtf > fusion_bloom.gtf
	bgzip fusion_bloom.gtf
	tabix fusion_bloom.gtf.gz
	rm tid_gtf tid_cdna comm_tid comm.gtf
fi
	

