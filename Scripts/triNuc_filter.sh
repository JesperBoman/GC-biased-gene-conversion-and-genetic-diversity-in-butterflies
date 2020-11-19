#!/bin/bash -l
#SBATCH -J trinuc_filter
#SBATCH -o trinuc_filter.out
#SBATCH -e trinuc_filter.err
#SBATCH --mail-user #Add email
#SBATCH --mail-type=ALL
#SBATCH -t 00-05:00:00
#SBATCH -A snic2020-5-20
#SBATCH -p core
#SBATCH -n 10

GENOME="leptidea_ancestral_JuveSpecial_unMasked.fa"

ml bioinfo-tools samtools/1.9

#This code parallellizes the process of finding ancestral CpG-prone sites

nrows=$(wc -l ancestral_state_list_outTough_JuveSpecial_noExons | cut -f1 -d ' ')

for ((i=0; i<=9; i++))
do
lower=$((i* nrows/ 10))
upper=$(((i+1) * nrows /10))
awk -v genome=$GENOME -f triNuc_filter.awk <(awk -v lower=$lower -v upper=$upper '{if (NR > lower && NR <= upper)print $0}' ancestral_state_list_outTough_JuveSpecial_noExons ) > temp${i} &
done
wait

for ((i=0; i<=9; i++))
do
awk '{print $0}' temp${i} >> ancestral_state_list_outTough_noExons_CpG_annot
rm temp${i}
done
