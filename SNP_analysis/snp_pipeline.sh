#RNA-Seq pipeline script

sgi=/home/nipgr/Documents/chickpea/genome_files/STAR_FUSION_CHICKPEA_LIB/ref_genome.fa.star.idx
fastp=/home/nipgr/software/fastp
gatk=/home/nipgr/software/gatk-4.2.0.0/gatk
bcftools=/home/nipgr/software/bcftools/bin/bcftools
genome=/home/nipgr/Documents/chickpea/genome_files/GCF_000331145.1_ASM33114v1_genomic.fna

#Aschochyta
#Control
for file in $(<Aschochyta/Aschochyta_Infected_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Aschochyta/control/ "{}" \;
$fastp --detect_adapter_for_pe -i Aschochyta/control/${file}_1.fastq -I Aschochyta/control/${file}_2.fastq -o Aschochyta/control/${file}_1_trim.fastq -O Aschochyta/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Aschochyta/control/${file}_1_trim.fastq Aschochyta/control/${file}_2_trim.fastq --outFileNamePrefix Aschochyta/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1     

rm Aschochyta/control/${file}_1.fastq Aschochyta/control/${file}_2.fastq Aschochyta/control/${file}_1_trim.fastq Aschochyta/control/${file}_2_trim.fastq
rm Aschochyta/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Aschochyta/control/${file}Aligned.sortedByCoord.out.bam -O Aschochyta/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Aschochyta -SM control

$gatk MarkDuplicates -I Aschochyta/control/${file}Aligned.sortedByCoord.RG.out.bam -O Aschochyta/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Aschochyta/control/${file}_output_metrics.txt

done

#Treated
for file in $(<Aschochyta/Aschochyta_Infected.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Aschochyta/infected/ "{}" \;
$fastp --detect_adapter_for_pe -i Aschochyta/infected/${file}_1.fastq -I Aschochyta/infected/${file}_2.fastq -o Aschochyta/infected/${file}_1_trim.fastq -O Aschochyta/infected/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Aschochyta/infected/${file}_1_trim.fastq Aschochyta/infected/${file}_2_trim.fastq --outFileNamePrefix Aschochyta/infected/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Aschochyta/infected/${file}_1.fastq Aschochyta/infected/${file}_2.fastq Aschochyta/infected/${file}_1_trim.fastq Aschochyta/infected/${file}_2_trim.fastq
rm Aschochyta/infected/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Aschochyta/infected/${file}Aligned.sortedByCoord.out.bam -O Aschochyta/infected/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Aschochyta -SM infected

$gatk MarkDuplicates -I Aschochyta/infected/${file}Aligned.sortedByCoord.RG.out.bam -O Aschochyta/infected/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Aschochyta/infected/${file}_output_metrics.txt

done

#CO2_stress
#Control
for file in $(<CO2_stress/CO2_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O CO2_stress/control/ "{}" \;
$fastp --detect_adapter_for_pe -i CO2_stress/control/${file}_1.fastq -I CO2_stress/control/${file}_2.fastq -o CO2_stress/control/${file}_1_trim.fastq -O CO2_stress/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn CO2_stress/control/${file}_1_trim.fastq CO2_stress/control/${file}_2_trim.fastq --outFileNamePrefix CO2_stress/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm CO2_stress/control/${file}_1.fastq CO2_stress/control/${file}_2.fastq CO2_stress/control/${file}_1_trim.fastq CO2_stress/control/${file}_2_trim.fastq
rm CO2_stress/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I CO2_stress/control/${file}Aligned.sortedByCoord.out.bam -O CO2_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB CO2_stress -SM control

$gatk MarkDuplicates -I CO2_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -O CO2_stress/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M CO2_stress/control/${file}_output_metrics.txt

done

#stress
for file in $(<CO2_stress/CO2_stress.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O CO2_stress/stress/ "{}" \;
$fastp --detect_adapter_for_pe -i CO2_stress/stress/${file}_1.fastq -I CO2_stress/stress/${file}_2.fastq -o CO2_stress/stress/${file}_1_trim.fastq -O CO2_stress/stress/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn CO2_stress/stress/${file}_1_trim.fastq CO2_stress/stress/${file}_2_trim.fastq --outFileNamePrefix CO2_stress/stress/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm CO2_stress/stress/${file}_1.fastq CO2_stress/stress/${file}_2.fastq CO2_stress/stress/${file}_1_trim.fastq CO2_stress/stress/${file}_2_trim.fastq
rm CO2_stress/stress/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I CO2_stress/stress/${file}Aligned.sortedByCoord.out.bam -O CO2_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB CO2_stress -SM stress

$gatk MarkDuplicates -I CO2_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -O CO2_stress/stress/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M CO2_stress/stress/${file}_output_metrics.txt

done

#Drought_stress
#Control
for file in $(<Drought_stress/Drought_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Drought_stress/control/ "{}" \;
$fastp --detect_adapter_for_pe -i Drought_stress/control/${file}_1.fastq -I Drought_stress/control/${file}_2.fastq -o Drought_stress/control/${file}_1_trim.fastq -O Drought_stress/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Drought_stress/control/${file}_1_trim.fastq Drought_stress/control/${file}_2_trim.fastq --outFileNamePrefix Drought_stress/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Drought_stress/control/${file}_1.fastq Drought_stress/control/${file}_2.fastq Drought_stress/control/${file}_1_trim.fastq Drought_stress/control/${file}_2_trim.fastq
rm Drought_stress/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Drought_stress/control/${file}Aligned.sortedByCoord.out.bam -O Drought_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Drought_stress -SM control

$gatk MarkDuplicates -I Drought_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -O Drought_stress/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Drought_stress/control/${file}_output_metrics.txt

done

#Stress
for file in $(<Drought_stress/Drought_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Drought_stress/stress/ "{}" \;
$fastp --detect_adapter_for_pe -i Drought_stress/stress/${file}_1.fastq -I Drought_stress/stress/${file}_2.fastq -o Drought_stress/stress/${file}_1_trim.fastq -O Drought_stress/stress/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Drought_stress/stress/${file}_1_trim.fastq Drought_stress/stress/${file}_2_trim.fastq --outFileNamePrefix Drought_stress/stress/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Drought_stress/stress/${file}_1.fastq Drought_stress/stress/${file}_2.fastq Drought_stress/stress/${file}_1_trim.fastq Drought_stress/stress/${file}_2_trim.fastq
rm Drought_stress/stress/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Drought_stress/stress/${file}Aligned.sortedByCoord.out.bam -O Drought_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Drought_stress -SM stress

$gatk MarkDuplicates -I Drought_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -O Drought_stress/stress/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Drought_stress/stress/${file}_output_metrics.txt

done

#Foc1
#Control
for file in $(<Foc1/Foc1_Infected_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Foc1/control/ "{}" \;
$fastp --detect_adapter_for_pe -i Foc1/control/${file}_1.fastq -I Foc1/control/${file}_2.fastq -o Foc1/control/${file}_1_trim.fastq -O Foc1/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Foc1/control/${file}_1_trim.fastq Foc1/control/${file}_2_trim.fastq --outFileNamePrefix Foc1/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Foc1/control/${file}_1.fastq Foc1/control/${file}_2.fastq Foc1/control/${file}_1_trim.fastq Foc1/control/${file}_2_trim.fastq
rm Foc1/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Foc1/control/${file}Aligned.sortedByCoord.out.bam -O Foc1/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Foc1 -SM control

$gatk MarkDuplicates -I Foc1/control/${file}Aligned.sortedByCoord.RG.out.bam -O Foc1/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Foc1/control/${file}_output_metrics.txt

done

#Infected
for file in $(<Foc1/Foc1_Infected.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Foc1/infected/ "{}" \;
$fastp --detect_adapter_for_pe -i Foc1/infected/${file}_1.fastq -I Foc1/infected/${file}_2.fastq -o Foc1/infected/${file}_1_trim.fastq -O Foc1/infected/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Foc1/infected/${file}_1_trim.fastq Foc1/infected/${file}_2_trim.fastq --outFileNamePrefix Foc1/infected/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Foc1/infected/${file}_1.fastq Foc1/infected/${file}_2.fastq Foc1/infected/${file}_1_trim.fastq Foc1/infected/${file}_2_trim.fastq
rm Foc1/infected/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Foc1/infected/${file}Aligned.sortedByCoord.out.bam -O Foc1/infected/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Foc1 -SM infected

$gatk MarkDuplicates -I Foc1/infected/${file}Aligned.sortedByCoord.RG.out.bam -O Foc1/infected/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Foc1/infected/${file}_output_metrics.txt

done

#Fusarium
#Control
for file in $(<Fusarium/Fusarium_Infected_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Fusarium/control/ "{}" \;
$fastp --detect_adapter_for_pe -i Fusarium/control/${file}_1.fastq -I Fusarium/control/${file}_2.fastq -o Fusarium/control/${file}_1_trim.fastq -O Fusarium/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Fusarium/control/${file}_1_trim.fastq Fusarium/control/${file}_2_trim.fastq --outFileNamePrefix Fusarium/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Fusarium/control/${file}_1.fastq Fusarium/control/${file}_2.fastq Fusarium/control/${file}_1_trim.fastq Fusarium/control/${file}_2_trim.fastq
rm Fusarium/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Fusarium/control/${file}Aligned.sortedByCoord.out.bam -O Fusarium/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Fusarium -SM control

$gatk MarkDuplicates -I Fusarium/control/${file}Aligned.sortedByCoord.RG.out.bam -O Fusarium/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Fusarium/control/${file}_output_metrics.txt

done

#Infected
for file in $(<Fusarium/Fusarium_Infected.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O Fusarium/infected/ "{}" \;
$fastp --detect_adapter_for_pe -i Fusarium/infected/${file}_1.fastq -I Fusarium/infected/${file}_2.fastq -o Fusarium/infected/${file}_1_trim.fastq -O Fusarium/infected/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn Fusarium/infected/${file}_1_trim.fastq Fusarium/infected/${file}_2_trim.fastq --outFileNamePrefix Fusarium/infected/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm Fusarium/infected/${file}_1.fastq Fusarium/infected/${file}_2.fastq Fusarium/infected/${file}_1_trim.fastq Fusarium/infected/${file}_2_trim.fastq
rm Fusarium/infected/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I Fusarium/infected/${file}Aligned.sortedByCoord.out.bam -O Fusarium/infected/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB Fusarium -SM infected

$gatk MarkDuplicates -I Fusarium/infected/${file}Aligned.sortedByCoord.RG.out.bam -O Fusarium/infected/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M Fusarium/infected/${file}_output_metrics.txt

done

#H.armigera
#Control
for file in $(<H.armigera/H.armigera_infected_and_wounded_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O H.armigera/control/ "{}" \;
$fastp --detect_adapter_for_pe -i H.armigera/control/${file}_1.fastq -I H.armigera/control/${file}_2.fastq -o H.armigera/control/${file}_1_trim.fastq -O H.armigera/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn H.armigera/control/${file}_1_trim.fastq H.armigera/control/${file}_2_trim.fastq --outFileNamePrefix H.armigera/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm H.armigera/control/${file}_1.fastq H.armigera/control/${file}_2.fastq H.armigera/control/${file}_1_trim.fastq H.armigera/control/${file}_2_trim.fastq
rm H.armigera/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I H.armigera/control/${file}Aligned.sortedByCoord.out.bam -O H.armigera/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB H.armigera -SM control

$gatk MarkDuplicates -I H.armigera/control/${file}Aligned.sortedByCoord.RG.out.bam -O H.armigera/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M H.armigera/control/${file}_output_metrics.txt

done

#Infected
for file in $(<H.armigera/H.armigera_infected_and_wounded.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O H.armigera/infected/ "{}" \;
$fastp --detect_adapter_for_pe -i H.armigera/infected/${file}_1.fastq -I H.armigera/infected/${file}_2.fastq -o H.armigera/infected/${file}_1_trim.fastq -O H.armigera/infected/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn H.armigera/infected/${file}_1_trim.fastq H.armigera/infected/${file}_2_trim.fastq --outFileNamePrefix H.armigera/infected/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm H.armigera/infected/${file}_1.fastq H.armigera/infected/${file}_2.fastq H.armigera/infected/${file}_1_trim.fastq H.armigera/infected/${file}_2_trim.fastq
rm H.armigera/infected/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I H.armigera/infected/${file}Aligned.sortedByCoord.out.bam -O H.armigera/infected/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB H.armigera -SM infected

$gatk MarkDuplicates -I H.armigera/infected/${file}Aligned.sortedByCoord.RG.out.bam -O H.armigera/infected/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M H.armigera/infected/${file}_output_metrics.txt

done

#NaCl_stress
#Control
for file in $(<NaCl_stress/NaCl_Control.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O NaCl_stress/control/ "{}" \;
$fastp --detect_adapter_for_pe -i NaCl_stress/control/${file}_1.fastq -I NaCl_stress/control/${file}_2.fastq -o NaCl_stress/control/${file}_1_trim.fastq -O NaCl_stress/control/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn NaCl_stress/control/${file}_1_trim.fastq NaCl_stress/control/${file}_2_trim.fastq --outFileNamePrefix NaCl_stress/control/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm NaCl_stress/control/${file}_1.fastq NaCl_stress/control/${file}_2.fastq NaCl_stress/control/${file}_1_trim.fastq NaCl_stress/control/${file}_2_trim.fastq
rm NaCl_stress/control/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I NaCl_stress/control/${file}Aligned.sortedByCoord.out.bam -O NaCl_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB NaCl_stress -SM control

$gatk MarkDuplicates -I NaCl_stress/control/${file}Aligned.sortedByCoord.RG.out.bam -O NaCl_stress/control/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M NaCl_stress/control/${file}_output_metrics.txt

done

#Stress
for file in $(<NaCl_stress/NaCl_80mM.txt)
do
find /home/nipgr/Documents/chickpea/ -name "$file\.sra" -exec fasterq-dump -3 -p -O NaCl_stress/stress/ "{}" \;
$fastp --detect_adapter_for_pe -i NaCl_stress/stress/${file}_1.fastq -I NaCl_stress/stress/${file}_2.fastq -o NaCl_stress/stress/${file}_1_trim.fastq -O NaCl_stress/stress/${file}_2_trim.fastq -w 16

STAR --runThreadN 40 --genomeDir ${sgi} --readFilesIn NaCl_stress/stress/${file}_1_trim.fastq NaCl_stress/stress/${file}_2_trim.fastq --outFileNamePrefix NaCl_stress/stress/${file} --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 1

rm NaCl_stress/stress/${file}_1.fastq NaCl_stress/stress/${file}_2.fastq NaCl_stress/stress/${file}_1_trim.fastq NaCl_stress/stress/${file}_2_trim.fastq
rm NaCl_stress/stress/*.[ot][ua][tb]

$gatk AddOrReplaceReadGroups -I NaCl_stress/stress/${file}Aligned.sortedByCoord.out.bam -O NaCl_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -ID ${file} -PL ILLUMINA -PU unit1 -LB NaCl_stress -SM stress

$gatk MarkDuplicates -I NaCl_stress/stress/${file}Aligned.sortedByCoord.RG.out.bam -O NaCl_stress/stress/${file}Aligned.sortedByCoord.RG.dedupped.out.bam --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M NaCl_stress/stress/${file}_output_metrics.txt

done

################### MERGE BAM FILES ################### 
#Aschochyta/Aschochyta_Control
samtools merge Aschochyta_Control_merged.bam `ls Aschochyta/control/*RG.dedu*.bam`
#Aschochyta/Aschochyta_Infected
samtools merge Aschochyta_Infected_merged.bam `ls Aschochyta/infected/*RG.dedu*.bam`
#CO2_stress/CO2_Control
samtools merge CO2_Control_merged.bam `ls CO2_stress/control/*RG.dedu*.bam`
#CO2_stress/CO2_stress
samtools merge CO2_stress_merged.bam `ls CO2_stress/stress/*RG.dedu*.bam`
#Drought_stress/Drought_Control
samtools merge Drought_Control_merged.bam `ls Drought_stress/control/*RG.dedu*.bam`
#Drought_stress/Drought_Stress
samtools merge Drought_Stress_merged.bam `ls Drought_stress/stress/*RG.dedu*.bam`
#Foc1/Foc1_Infected_Control
samtools merge Foc1_Control_merged.bam `ls Foc1/control/*RG.dedu*.bam`
#Foc1/Foc1_Infected
samtools merge Foc1_Infected_merged.bam `ls Foc1/infected/*RG.dedu*.bam`
#Fusarium/Fusarium_Infected_Control
samtools merge Fusarium_Control_merged.bam `ls Fusarium/control/*RG.dedu*.bam`
#Fusarium/Fusarium_Infected
samtools merge Fusarium_Infected_merged.bam `ls Fusarium/infected/*RG.dedu*.bam`
#H.armigera/H.armigera_infected_and_wounded_Control
samtools merge H.armigera_control_merged.bam `ls H.armigera/control/*RG.dedu*.bam`
#H.armigera/H.armigera_infected_and_wounded
samtools merge H.armigera_infected_merged.bam `ls H.armigera/infected/*RG.dedu*.bam`
#NaCl_stress/NaCl_80mM
samtools merge NaCl_80mM_merged.bam `ls NaCl_stress/control/*RG.dedu*.bam`
#NaCl_stress/NaCl_Control
samtools merge NaCl_control_merged.bam `ls NaCl_stress/stress/*RG.dedu*.bam`

################### call variant ###################
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} Aschochyta_Control_merged.bam Aschochyta_Infected_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o Aschochyta_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} CO2_Control_merged.bam CO2_stress_merged.bam | ${bcftools} call -mv -Ob |${bcftools} view -i '%QUAL>=20' -Ov -o CO2_stress_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} Drought_Control_merged.bam Drought_Stress_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o Drought_stress_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} Foc1_Control_merged.bam Foc1_Infected_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o Foc1_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} Fusarium_Control_merged.bam Fusarium_Infected_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o Fusarium_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} H.armigera_control_merged.bam H.armigera_infected_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o H.armigera_calls.vcf
${bcftools} mpileup -E -Oua  -d10000000 -f ${genome} NaCl_control_merged.bam NaCl_80mM_merged.bam | ${bcftools} call -mv -Ob | ${bcftools} view -i '%QUAL>=20' -Ov -o NaCl_calls.vcf
