#! /bin/bash

# Author: Ajeet Singh
# Email: singh.ajeet@nipgr.ac.in

# This script is designed for:
#  1. prepare a BED files for all brekpoints in FuMa overlap results
#  2. get the fasta file from prepared BED file

#REF: https://doi.org/10.1093/nar/gkz1223 
#To identify motifs present in the 5′ and 3′ parental genes, we used the
#Gapped Local Alignment of Motifs (GLAM2) tool from the MEME SUITE with default parameters
#to find enriched motifs in the 200 bp upstream and downstream sequences from the 
#breakpoint position. Further, we used the Tomtom tool from MEME SUITE with default 
#parameters on the list of identified motifs and scanned a database of RNA binding protein motifs.

SRR=$1
fuma=$SRR/FuMA_overlap_no_if.txt
squid=$SRR/squid_out_sv.txt
mapsplice=$SRR/MapSplice_out/fusions_candidates.txt
starfusion=$SRR/STARFusion-out/star-fusion.fusion_predictions.abridged.coding_effect.tsv

#read file from second line(no header);separate multiple fusion ids into next line with intial tab; replace vacant field with previous value in a column, file here separted with tab

#Squid MapSplice STAR-Fusion Trinity Fusion-Bloom

#######start from 4th col; Squid##########
trm4=$(head -1 ${fuma}|cut -f4)
if [ "$trm4" == "Squid" ];then
	sed -n "2,$"p ${fuma} |\
	awk -F\\t '$4!=""{print $1"_"$2"\t"$4}'|\
	sed 's/,/\n\t/g'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(".*.=chr","",$2);gsub("-chr","\n\t",$2)}1' OFS='\t'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(":","\t",$2)}1' OFS='\t' \
	> four.1xxeeu
	
	#check if file is not empty and assign strand for squid results
	ch=$(find $pwd -empty -name "four.1xxeeu"|wc -l)
	if [ "$ch" -eq 0 ];then
		cat four.1xxeeu|while read i
		do
			printf "$i\t";
			gan=$(echo "$i"|cut -f2);
			gab=$(echo "$i"|cut -f3);
			grep "$gan" ${squid} |grep "$gab"|wc -l;
		done >four.a.2xxeeu
		awk -F\\t '{if($NF==0)gsub($NF,"+",$NF);else{gsub($NF,"-",$NF)}}1' OFS='\t' four.a.2xxeeu > four.b.2xxeeu
	fi
	rm four.a.2xxeeu four.1xxeeu
fi
#######5th col; MapSplice###############
trm5=$(head -1 ${fuma}|cut -f5)
if [ "$trm5" == "MapSplice" ];then
	sed -n "2,$"p ${fuma} |\
	awk -F\\t '$5!=""{print $1"_"$2"\t"$5}'|\
	sed 's/,/\n\t/g'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(".*.=chr","",$2);gsub("-chr","\n\t",$2)}1' OFS='\t'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(":","\t",$2)}1' OFS='\t' \
	> five.1xxeeu
	
	#check if file is not empty and assign strand for MapSplice results
	ch=$(find $pwd -empty -name "five.1xxeeu"|wc -l)
	if [ "$ch" -eq 0 ];then
		cat five.1xxeeu|while read i
		do
			gan=$(echo "$i"|cut -f2);
			gab=$(echo "$i"|cut -f3);
			grep -w "$gan" ${mapsplice} |awk -F\\t -v a="$gab" -v z="$i" '{split($6,b,""); for(i=1; i<=NF; ++i){if($i ~ a"$")if(i == 2 )print z,b[1]; else {print z,b[2]}} }' OFS='\t';
		done > five.b.2xxeeu
	fi
	rm five.1xxeeu
fi
			
#####6th col; STAR-Fusion###############
trm6=$(head -1 ${fuma}|cut -f6)
if [ "$trm6" == "STAR-Fusion" ];then
	sed -n "2,$"p ${fuma} |\
	awk -F\\t '$6!=""{print $1"_"$2"\t"$6}'|\
	sed 's/,/\n\t/g'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(".*.=chr","",$2);gsub("-chr","\n\t",$2)}1' OFS='\t'|\
	awk '/^\t/{$0=(x)substr($0,"\\t"length(x)+1)}{x=$1}1'|\
	awk -F\\t '{gsub(":","\t",$2)}1' OFS='\t' \
	> six.1xxeeu
	
	#check if file is not empty and assign strand for STAR-Fusion results
	ch=$(find $pwd -empty -name "six.1xxeeu"|wc -l)
	if [ "$ch" -eq 0 ];then
		cat six.1xxeeu|while read i
		do
			printf "$i\t";
			gan=$(echo "$i"|cut -f2);
			gab=$(echo "$i"|cut -f3);
			grep -w "$gan" ${starfusion}|grep -wo "$gab:."|sort -u|sed "s/$gab://g" ;
		done > six.b.2xxeeu
	fi
	rm six.1xxeeu
fi

# cat all and sort unique and create up down 200bp BED file for 5' and 3'
#downstream bases 200, from 5'end for fusion junction
cat *b.2xxeeu|sort -u|awk -F $'\t' '{if($4=="-"){print $2,$3-1,$3+200,$1,1,$4;} else if($4=="+")print $2,$3-200,$3,$1,1,$4}' OFS='\t' > 5prime.2xxeeu
#upstream 200 bases from 3' end for fusion junction
cat *b.2xxeeu|sort -u|awk -F $'\t' '{if($4=="-"){print $2,$3-200,$3,$1,1,$4;} else if($4=="+")print $2,$3-1,$3+200,$1,1,$4}' OFS='\t' > 3prime.2xxeeu
#replace start or end if start from -ve to 0
awk -F\\t '{if($2~/^-/)gsub($2,"0",$2);else if($3~/^-/)gsub($3,"0",$3)}1' OFS='\t' 5prime.2xxeeu > 3xxeeu5.bed
awk -F\\t '{if($2~/^-/)gsub($2,"0",$2);else if($3~/^-/)gsub($3,"0",$3)}1' OFS='\t' 3prime.2xxeeu > 3xxeeu3.bed

#cat 5prime.2xxeeu 3prime.2xxeeu |awk -F\\t '{if($2~/^-/)gsub($2,"0",$2);else if($3~/^-/)gsub($3,"0",$3)}1' OFS='\t' > 3xxeeu.bed

rm *b.2xxeeu 5prime.2xxeeu 3prime.2xxeeu

#get fasta file
#bedtools getfasta -fi ~/Documents/chickpea/genome_files/GCF_000331145.1_ASM33114v1_genomic.fna -bed 3xxeeu.bed -fo $1/$1\_up_down_200bp_5_3_prime.fa -name -s
bedtools getfasta -fi ~/Documents/chickpea/genome_files/GCF_000331145.1_ASM33114v1_genomic.fna -bed 3xxeeu5.bed -fo $1/$1\_up_200bp_5_prime.fa -name -s
bedtools getfasta -fi ~/Documents/chickpea/genome_files/GCF_000331145.1_ASM33114v1_genomic.fna -bed 3xxeeu3.bed -fo $1/$1\_down_200bp_3_prime.fa -name -s

rm 3xxeeu5.bed 3xxeeu3.bed

#/home/nipgr/software/meme-5.3.3/bin/glam2 n -o zch all_up_down_200bp_5_3_prime.fa
#/home/nipgr/software/meme-5.3.3/bin/tomtom glam2.meme glam2.meme
#tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 glam2.meme db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme db/EUKARYOTE/jolma2013.meme
