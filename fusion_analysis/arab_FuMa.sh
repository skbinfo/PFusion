#! /bin/bash
#singh.ajeet@nipgr.ac.in
#contributor malini.nk@nipgr.ac.in
GENE_ANNOTATION="/home/nipgr/Documents/arabidopsis/fusion_genomefiles/fuma_gtf.bed"
#fuma-gencode-gtf-to-bed /home/nipgr/Documents/chickpea/genome_files/GCF_000331145.1_ASM33114v1_genomic.gtf > $GENE_ANNOTATION

ERIC_IN=$1/fusion_eric_arabidopsis.results.filtered.tsv
#ERIC_IN=$1/fusion_eric_arabidopsis.results.total.txt
SQUID_IN=$1/squid_out_sv.txt
SQUID_OUT=$1/squid_fuma_uuqZ1Q.bedpe
#MAPSPLICE_IN=$1/fusions_candidates.txt
#MAPSPLICE_OUT=$1/fusions_candidates_uuqZ1Q.txt
STAR_IN=$1/star-fusion.fusion_predictions.tsv
TRINITY_IN=$1/TrinityFusion-UC.fusion_predictions.tsv
#FUSION_BLOOM_IN=$1/sv.bedpe
#FUSION_BLOOM_OUT=$1/sv_uuqZ1Q.bedpe
FUMA_IN=""
LABEL=""

if [ -f "$ERIC_IN" ]; then
FUMA_IN+="EricScript:ericscript:$ERIC_IN"
LABEL+="EricScript:TAIR10"
else
echo "ERIC output does not exist for $1"
fi

if [ -f "$SQUID_IN" ]; then
grep -vi "#" $SQUID_IN|awk -F\\t -v 'col=1' '{gsub(".","UID"col,$7);col+=1}1' OFS='\t' >$SQUID_OUT
FUMA_IN+=" Squid:chimerascan:$SQUID_OUT"
LABEL+=" Squid:TAIR10"
else
echo "SQUID output does not exist for $1"
fi

#if [ -f "$MAPSPLICE_IN" ]; then
#awk -F "\t" '{split($1,a,"~");split($6,b,"");OFS="\t";$1="gene_1\tgene_2\t"a[1]"\t"$2"\t"b[1]"\t"a[2]"\t"$3"\t"b[2] FS $1;}1' $MAPSPLICE_IN|cut --complement -f9,10,11,14 > $MAPSPLICE_OUT
#cut -f4 --complement $MAPSPLICE_IN >$MAPSPLICE_OUT
#FUMA_IN+=" MapSplice:tophat-fusion_post_result:$MAPSPLICE_OUT"
#FUMA_IN+=" MapSplice:ericscript:$MAPSPLICE_OUT"
#LABEL+=" MapSplice:TAIR10"
#else
#echo "MapSplice output does not exist for $1"
#fi

if [ -f "$STAR_IN" ];then
FUMA_IN+=" STAR-Fusion:star-fusion_final:$STAR_IN"
LABEL+=" STAR-Fusion:TAIR10"
else
echo "STAR-Fusion output does not exist for $1"
fi

if [ -f "$TRINITY_IN" ];then
FUMA_IN+=" Trinity:trinity-gmap:$TRINITY_IN"
LABEL+=" Trinity:TAIR10"
else
echo "Trinity-fusion output does not exist for $1"
fi

#if [ -f "$FUSION_BLOOM_IN" ]; then
#grep -vi "#" $FUSION_BLOOM_IN >$FUSION_BLOOM_OUT
#FUMA_IN+=" Fusion-Bloom:chimerascan:$FUSION_BLOOM_OUT"
#LABEL+=" Fusion-Bloom:TAIR10"
#else
#echo "FusionBloom output does not exist for $1"
#fi

echo "${FUMA_IN}"
echo ""
echo "${LABEL}"
echo ""
#exact gene match
fuma \
    -a TAIR10:${GENE_ANNOTATION} \
    -m egm \
    --strand-specific-matching \
    --no-acceptor-donor-order-specific-matching \
    -s ${FUMA_IN} \
    -l ${LABEL} \
    -o $1/FuMA_egm.txt \
    --verbose

#overlap
fuma \
    -a TAIR10:${GENE_ANNOTATION} \
    -m overlap \
    --strand-specific-matching \
    --no-acceptor-donor-order-specific-matching \
    -s ${FUMA_IN} \
    -l ${LABEL} \
    -o $1/FuMA_overlap.txt \
    --verbose
find $1/ -iname "*uuqZ1Q*" -exec rm -f "{}" \;
