#! /bin/bash

Genome="PATH TO respective_genome_file ############"

pdir="/mnt/nipgr_user/motif_databases/CISBP-RNA/plant"
plant_motif_dbs="${pdir}/Oryza_sativa.dna_encoded.meme ${pdir}/Cucumis_sativus.dna_encoded.meme ${pdir}/Glycine_max.dna_encoded.meme ${pdir}/Prunus_persica.dna_encoded.meme ${pdir}/Oryza_indica.dna_encoded.meme ${pdir}/Populus_trichocarpa.dna_encoded.meme ${pdir}/Manihot_esculenta.dna_encoded.meme ${pdir}/Zea_mays.dna_encoded.meme ${pdir}/Ricinus_communis.dna_encoded.meme ${pdir}/Arabidopsis_lyrata.dna_encoded.meme ${pdir}/Medicago_truncatula.dna_encoded.meme ${pdir}/Sorghum_bicolor.dna_encoded.meme ${pdir}/Arabidopsis_thaliana.dna_encoded.meme ${pdir}/Malus_x_domestica.dna_encoded.meme ${pdir}/Lotus_japonicus.dna_encoded.meme ${pdir}/Vitis_vinifera.dna_encoded.meme ${pdir}/Physcomitrella_patens.dna_encoded.meme"

glam2="/home/nipgr/software/meme-5.3.3/bin/glam2"
tomtom="/home/nipgr/software/meme-5.3.3/bin/tomtom"

resline=''

FusionName () {
        echo "$resline"|cut -f2,3|awk '{print $1"_"$2}'
}

line_to_coordinate_read () {
        echo "$resline"|cut -f1,2,3,4 --complement\
                | sed 's/^\s\{1,100\}//g'|sed 's/,/\t/g'\
                | sed 's/chr\|UID//g'| sed 's/\t/\n/g' \
                | sed '/^$/d' | sed 's/.*=//g' | sed 's/:\|-/\t/g'
}

CordToBED () {
        echo "$1" | awk -F\\t '{print $2,$3-100,$3+100,$1"_1","1"}' OFS='\t'
        echo "$1" | awk -F\\t '{print $4,$5-100,$5+100,$1"_2","1"}' OFS='\t'
}

TypeSeparate () {
        FuNa=$(FusionName "$resline")
        line_to_coordinate_read "$resline" | sed "s/^/$FuNa\t/g" \
                |sort -u | while read corrd
do
        CordToBED "$corrd"
done
}

cat $1|grep -v "Left-genes"| while read resline; do TypeSeparate "$resline"; done > uuq.bed

grep "_1" uuq.bed | sed '/\t-[0-9]/d' > gene1.bed
grep "_2" uuq.bed | sed '/\t-[0-9]/d' > gene2.bed

for bed in `ls gene1.bed gene2.bed`
do
	bname=$(basename $bed .bed)
	bedtools getfasta -fi ${Genome} -bed ${bname}.bed -name -fo ${bname}.fa
	${glam2} -2 n -o ${bname}_glam2_output ${bname}.fa
	${tomtom} -norc -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10 ${bname}_glam2_output/glam2.meme ${plant_motif_dbs}
done
rm uuq.bed
