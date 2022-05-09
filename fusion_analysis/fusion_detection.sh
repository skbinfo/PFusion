#! /bin/bash

# Author Ajeet Singh
# Email: singh.ajeet@nipgr.ac.in
# Contributor: Malini Nemalikanti, Rohan Bhardwaj

# This script detects the fusion from RNA-Seq samples 
# Usage: fusion_analysis.sh SAMPLE_ID[SRR_NUM]


srr_id="$1"
group="$2"

echo $'\n'"RUN STARTED FOR $srr_id"$'\n'
#Tools
CTAT_BUILD="perl /mnt/storage/sklab202/software/STAR-Fusion-v1.10.0/ctat-genome-lib-builder/prep_genome_lib.pl"
StarFusion="/mnt/storage/sklab202/software/STAR-Fusion-v1.10.0/STAR-Fusion"
Star="/mnt/storage/sklab202/software/STAR-2.7.9a/bin/Linux_x86_64/STAR"
SQUID="/mnt/storage/sklab202/software/squid-v1.5_linux_x86_64/squid"
Mapsplice="/mnt/storage/sklab202/software/MapSplice-v2.2.1/mapsplice.py"
FusionBloom="/mnt/storage/sklab202/software/pavfinder_env/bin/fusion-bloom"
TrinityFusion="/mnt/storage/sklab202/software/TrinityFusion-v0.3.5/TrinityFusion"
fastp="/mnt/storage/sklab202/software/fastp"
fasterq_dump="/mnt/storage/sklab202/software/sratoolkit.2.11.0-centos_linux64/bin/fasterq-dump"
ericscript="perl /mnt/storage/sklab202/software/EricScript-Plants/EricScript-Plants/ericscript.pl"

#dir
CTAT_BUILD_DIR="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/CTAT_Library/ctat_genome_lib_build_dir"
star_index_dir="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/CTAT_Library/ctat_genome_lib_build_dir/ref_genome.fa.star.idx"
mapsplice_chromosome_files="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/mapsplice_fasta_files"
mapsplice_index_files="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/index_files_mapsplice"

#files
genome="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/genome.fa" # created as symlink from /mnt/ISOLON/CHIP_Seq_Rice
gtf="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/genome.gtf" # created as symlink from /mnt/ISOLON/CHIP_Seq_Rice
dfam="/mnt/ISOLON/Fusion_Oryza/domain_DB/Dfam_curatedonly.hmm" #  https://www.dfam.org/releases/Dfam_3.5/families/Dfam_curatedonly.hmm.gz
pfam="/mnt/ISOLON/Fusion_Oryza/domain_DB/Pfam-A.hmm" #ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
fusionbloom_profile="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_genome_files/fusion-bloom_profile"
run_dir="/mnt/ISOLON/Fusion_Oryza/oryza_${group}_fusion_result/$srr_id"

#Check for results

if [ -d "$run_dir/STARFusion-out" ] && [ -f "$run_dir/squid_out_sv.txt" ] &&  [ -d "$run_dir/mapsplice_out" ] && [ -d "$run_dir/FusionBloom_out" ] && [ -d "$run_dir/TrinityFusion_out" ] && [ -d "$run_dir/ericscript_out" ]; then
	echo $'\n'"Results already present for $srr_id"$'\n'
	exit 1
fi

#fastp_files

if [ -f "$run_dir/reads/${srr_id}_1_trim.fastq" ] && [ -f "$run_dir/reads/${srr_id}_2_trim.fastq" ];then
	echo $'\n'"fastp files exist for $srr_id"$'\n'

else
	echo $'\n'"###### RUNNING fasterq-dump ######"$'\n'
	$fasterq_dump -S /mnt/ISOLON/Fusion_Oryza/oryza_${group}_RNA_Seq_data/$srr_id/${srr_id}.sra -O $run_dir/reads
                if [ "$?" -ne "0" ]; then
                        echo $'\n'"fasterq-dump failed for $srr_id at line $LINENO"$'\n'
			echo $'\n'"Exit"$'\n'
                        exit 1
                fi

	echo $'\n'"###### RUNNING fastp ######"$'\n' 
	$fastp --detect_adapter_for_pe -i $run_dir/reads/${srr_id}_1.fastq -I $run_dir/reads/${srr_id}_2.fastq -o $run_dir/reads/${srr_id}_1_trim.fastq -O $run_dir/reads/${srr_id}_2_trim.fastq -h $run_dir/reads/${srr_id}.html -j $run_dir/reads/${srr_id}.json -w 16
                if [ "$?" -ne "0" ]; then
                        echo $'\n'"fastp failed for $srr_id at line $LINENO"$'\n'
                        echo $'\n'"Exit"$'\n'
                        exit 1
                fi

	rm -f $run_dir/reads/${srr_id}_1.fastq $run_dir/reads/${srr_id}_2.fastq

	left_fq_filename=$run_dir/reads/${srr_id}_1_trim.fastq
	right_fq_filename=$run_dir/reads/${srr_id}_2_trim.fastq
fi

#CTAT_Library

if [ -d "${CTAT_BUILD_DIR}" ];then
        echo $'\n'"CTAT Genome exists at ${CTAT_BUILD_DIR}"$'\n'
        if [ -d "${star_index_dir}" ];then
                echo $'\n'"STAR index exists at ${star_index_dir}"$'\n'
        else
                echo $'\n'"Please generate STAR index for genome"$'\n'
        fi
else
        ${CTAT_BUILD} --genome_fa ${genome} --gtf ${gtf} --dfam_db ${dfam} --output_dir ${CTAT_BUILD_DIR} --pfam_db ${pfam} --CPU 110 --gmap_build
fi

#Tools

echo $'\n'"RUNNING FUSION_DETECTION_TOOLSET for $1"$'\n'

echo $'\n'"###### RUNNING STAR-FUSION ######"$'\n'
if [ -d "$run_dir/STARFusion-out" ]; then
	echo $'\n'"STARFusion-out exists for $srr_id"$'\n'
else
	$StarFusion --left_fq ${left_fq_filename} --right_fq ${right_fq_filename} --genome_lib_dir ${CTAT_BUILD_DIR} --CPU 100 --output_dir $run_dir/STARFusion-out --FusionInspector validate --examine_coding_effect
		if [ "$?" -ne "0" ]; then
			 echo $'\n'"STARFusion run failed for $srr_id at line $LINENO"$'\n'
			rm -rf $run_dir/STARFusion-out
			echo $'\n'"Exit"$'\n'
  			exit 1
		fi
	echo $'\n'"Star-Fusion run completed"$'\n'
#--denovo_reconstruct
# If FusionInspector 'inspect' mode is invoked, then only the fusion-evidence reads are de novo assembled. If FusionInspector 'validate' mode is selected, then all reads aligned to the fusion gene contigs are assembled.
fi

echo $'\n' "###### RUNNING SQUID ######"$'\n'
if [ -f "$run_dir/squid_out_sv.txt" ]; then
	echo $'\n'"squid_out exists for $srr_id"$'\n'
else
	$Star --runThreadN 100 --genomeDir ${star_index_dir} --readFilesIn ${left_fq_filename} ${right_fq_filename} --outFileNamePrefix $run_dir/star_run. --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --chimOutType SeparateSAMold
	echo $'\n'RUNNING SAMTOOLS$'\n'
	samtools view -Shb --threads 100 $run_dir/star_run.Chimeric.out.sam -o $run_dir/star_run.Chimeric.out.bam
	echo $'\n'RUNNING SQUID$'\n'
	$SQUID -b $run_dir/star_run.Aligned.sortedByCoord.out.bam -c $run_dir/star_run.Chimeric.out.bam -o $run_dir/squid_out
	rm -f $run_dir/star_run*
		if [ "$?" -ne "0" ]; then
         		echo $'\n'"Squid run failed for $srr_id at line $LINENO"$'\n'
        		rm -f $run_dir/squid_out_sv.txt
        		exit 1
        	fi
        echo $'\n'"Squid run completed"$'\n'

fi

echo $'\n'"###### RUNNING MAPSPLICE2 ######"$'\n'
if [ -d "$run_dir/mapsplice_out" ]; then
	echo $'\n'"mapsplice_out exists for $srr_id"$'\n'
else
	python2 $Mapsplice -1 $left_fq_filename -2 $right_fq_filename  --fusion -c $mapsplice_chromosome_files -x $mapsplice_index_files/${group}_index -p 2 -o $run_dir/mapsplice_out 2>$run_dir/mapsplice_log.txt
	rm -f $run_dir/mapsplice_out/alignments.sam
		if [ "$?" -ne "0" ]; then
                        echo $'\n'"Mapsplice run failed for $srr_id at line $LINENO"$'\n'
                        rm -rf $run_dir/mapsplice_out
                        exit 1
                fi

	echo $'\n'"Mapsplice2 run completed for $srr_id"$'\n'
fi

echo $'\n'"###### RUNNING FUSION-BLOOM ######"$'\n'
if [ -d "$run_dir/FusionBloom_out" ]; then
	echo $'\n'"FusionBloom_out exists for $srr_id"$'\n'
else 
	echo $'\n'"Calculating read length"$'\n'
	read1_length=$(sed -n '2~4'p $left_fq_filename|head -100|awk '{print length($1)}'|sort -n|tail -n 1)
	read2_length=$(sed -n '2~4'p $right_fq_filename|head -100|awk '{print length($1)}'|sort -n|tail -n 1)
	read_length=$(printf "$read1_length\n$read2_length\n"|sort -u|sort -n|tail -n 1)
	echo $'\n'"Read length is ${read_length}"$'\n'
	source $fusionbloom_profile
$FusionBloom profile=$fusionbloom_profile left=$left_fq_filename right=$right_fq_filename readlen=$read_length outdir=$run_dir/FusionBloom_out name=at_fusion
		if [ "$?" -ne "0" ]; then
                        echo $'\n'"Fusion Bloom run failed for $srr_id at line $LINENO"$'\n'
                        rm -rf $run_dir/FusionBloom_out
                        exit 1
                fi

	echo $'\n'"FusionBloom run completed"$'\n'
fi

echo $'\n'"###### RUNNING TRINITY-FUSION ######"$'\n'
if [ -d "$run_dir/TrinityFusion_out" ]; then
	echo $'\n'"TrinityFusion_out exists for $srr_id"$'\n'
else 
	#for arabidopsis: export CTAT_GENOME_LIB=$CTAT_GENOME_LIB:/mnt/storage/sklab202/genome_files/CTAT_Library/ctat_genome_lib_build_dir/
	$TrinityFusion --left_fq ${left_fq_filename} --right_fq ${right_fq_filename} --chimeric_junctions $run_dir/STARFusion-out/Chimeric.out.junction --aligned_bam $run_dir/STARFusion-out/Aligned.out.bam --genome_lib_dir ${CTAT_BUILD_DIR} --output_dir $run_dir/TrinityFusion_out --CPU 55
		if [ "$?" -ne "0" ]; then
                        echo $'\n'"Trinity-Fusion run failed for $srr_id at line $LINENO"$'\n'
                        rm -rf $run_dir/TrinityFusion_out
                        exit 1
                fi
	

	echo $'\n'"Trinity-Fusion run completed"$'\n'
fi

if [ -d "$run_dir/ericscript_out" ]; then
        echo $'\n'"ericscript_out exists for $srr_id"$'\n'
else
	if [ "$group" == "indica" ];then
		#perl /mnt/storage/sklab202/software/EricScript-Plants/EricScript-Plants/ericscript.pl --downdb --refid oryza_indica
		${ericscript} -p 100 --refid oryza_indica -name ${srr_id} -o $run_dir/ericscript_out ${left_fq_filename} ${right_fq_filename} > $run_dir/ericscript_log.txt
	elif [ "$group" == "japonica" ];then
		#perl /mnt/storage/sklab202/software/EricScript-Plants/EricScript-Plants/ericscript.pl --downdb --refid oryza_sativa
		${ericscript} -p 100 --refid oryza_sativa -name ${srr_id} -o $run_dir/ericscript_out ${left_fq_filename} ${right_fq_filename} > $run_dir/ericscript_log.txt
	fi
	rm -rf $run_dir/ericscript_out/aln
                if [ "$?" -ne "0" ]; then
                        echo $'\n'"Ericscript-Plants run failed for $srr_id at line $LINENO"$'\n'
                        rm -rf $run_dir/ericscript_out
                        exit 1
                fi

        echo $'\n'"Ericscript-Plants run completed"$'\n'
fi
rm -f $left_fq_filename $right_fq_filename
