#! /bin/bash

resline=''
genome=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

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
        echo "$1" | awk -F\\t '{print $2,$3-200,$3,$1,"1"}' OFS='\t' > cord.bed
        echo "$1" | awk -F\\t '{print $4,$5,$5+200,$1,"1"}' OFS='\t' >> cord.bed
	bedtools getfasta -fi $genome -bed cord.bed -tab | cut -f2 |tr '\n' '\t' \
		| sed 's/\t$/\n/g'|sed 's/\t//g'
}

TypeSeparate () {
	FuNa=$(FusionName "$resline")
	line_to_coordinate_read "$resline" | sed "s/^/$FuNa\t/g" \
		|sort -u | while read corrd
do
	echo "$corrd"
	fna=$(echo "$corrd"|cut -f1)
	chr1=$(echo "$corrd"|cut -f2)
	chr2=$(echo "$corrd"|cut -f4)
	pos1=$(echo "$corrd"|cut -f3)
	pos2=$(echo "$corrd"|cut -f5)
	if [[ "$chr1" -ne "$chr2" ]];then
		CordToBED "$corrd" | sed "1 i >$fna"  >> interType
	elif [[ "$chr1" -eq "$chr2" ]];then
		gap=$(echo "$pos2-$pos1"|bc)
		if [[ "$gap" -ge 200000 ]];then
			CordToBED "$corrd" | sed "1 i >$fna" >> distalType
		elif [[ "$gap" -lt 200000 ]];then
			CordToBED "$corrd" | sed "1 i >$fna" >> proximalType
		fi
	fi
done

}


cat $1| while read resline; do TypeSeparate "$resline"; done


