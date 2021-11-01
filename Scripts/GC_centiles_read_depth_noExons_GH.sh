#!/bin/bash -l

# Usage: Bash GC_centiles_read_depth_noExons.sh
# ==========================================================================================================
# A pipeline to determine the read coverage for GC centiles
# ==========================================================================================================
# Jesper Boman                      28 aug 2020
# ==========================================================================================================
# Extra info:
# See gBGC_main_pipeline.txt for information on how dependent files were created.
# 
# Script dependencies:
# annotateGaps.py - Made by Valentina Peona
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Updated 26/8 - 2020 for noExons filtering instead of noGenes and changed so bases in the reference genome with no mapping reads count as well

population="swe_sin"
skipn=$(wc -l ${population}_noExons_GC_scafWind | cut -f1 -d ' ' | awk '{print $1%100}')
SNPs_per_C=$(wc -l ${population}_noExons_GC_scafWind | cut -f1 -d ' ' | awk '{print int($1/100)}')

#Create GC centile SNP lists
awk -v sk="$skipn" 'NR>sk{print $0}' ${population}_noExons_GC_scafWind | awk -v SNPs="$SNPs_per_C" -v pop="$population" 'NR==1{i=1}{print >> pop"_noExons_GC"i; if(NR%SNPs==0){i=i+1}}'



#In main directory of project
module load bioinfo-tools biopython/1.73-py3
python3 annotateGaps.py ancestral_state_inference/leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked.fa >> leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_masked_ranges.bed


module load bioinfo-tools BEDTools/2.27.1
bedtools complement -i leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_masked_ranges.bed \
	-g <(awk '{print $1 "\t" $2}' ancestral_state_inference/leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked.fa.fai)\
	 > leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_nonMasked_ranges.bed




bamdir= #Add directory to BAM files here 
gc1kb= #Add GC1 kb files here, can be created using bedtools nuc on a bed file containing windows e.g. 1 kb

mkdir ${population}_genome_cov_files
cd ${population}_genome_cov_files

while IFS= read -r ind
do
        bedtools genomecov -ibam ${bamdir}/${ind}.recal.bam -d |awk '{if($3 != 0)print $1 "\t" $2-1 "\t" $2 "\t" $3}' > ${ind}_genome_cov.bed &

done < "../${population}_inds"
wait

cd ..


for ((c=1; c<=100; c++))
do

if [ $c -eq 1 ]
then
        headGC=0
        tailGC=$(tail -n1 ${population}_noExons_GC${c} | cut -f10)
elif [ $c -eq 100 ]
then
        headGC=$(head -n1 ${population}_noExons_GC${c} | cut -f10)
        tailGC=1
else
        headGC=$(head -n1 ${population}_noExons_GC${c} | cut -f10)
        tailGC=$(tail -n1 ${population}_noExons_GC${c} | cut -f10)
fi





awk -v headGC="$headGC" -v tailGC="$tailGC" 'FNR==NR{if($1 !~ /mtANA/  && $10 >= headGC && $10 <= tailGC){a[$11]}} FNR!=NR{if ($4 in a){print $1 "\t" $2 "\t" $3}}'  ${population}_noExons_GC_scafWind $gc1kb  > 1kb_wind_temp.bed

bedtools intersect -a 1kb_wind_temp.bed -b ../leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_nonMasked_ranges.bed > ranges_temp.bed
	
range_len=$(awk '{len+=($3-$2)}END{print len}' ranges_temp.bed)

while IFS= read -r ind
do
	bedtools intersect -a ${population}_genome_cov_files/${ind}_genome_cov.bed -b ranges_temp.bed -wa |awk -v ind="$ind" -v c="$c" -v len="$range_len"  '{sum+=$4} END{print ind "\t" "C_" c "\t" sum "\t" NR "\t" sum/len}' >> ${population}_centile_cov_noExons.results &
sleep 1m

done < "${population}_inds"
wait


done