#!/bin/bash -l
#SBATCH -J ancestral_genome
#SBATCH -o ancestral_genome.out
#SBATCH -e ancestral_genome.err
#SBATCH --mail-user #add email
#SBATCH --mail-type=ALL
#SBATCH -t 00-02:00:00
#SBATCH -A snic2020-5-20
#SBATCH -p core 
#SBATCH -n 2

#This script creates an ancestral genome sequence (.fasta) based on a list of inferred ancestral alleles 

GENOME="../N.Backstrom_leptidea.scf.1.4.fasta"

ml bioinfo-tools BEDTools/2.27.1

SNP_list="ancestral_state_list_outTough_ALT_sorted_JuveSpecial_biAllelic"
grep -w --no-group-separator "A" $SNP_list > snps_A

n=1
awk '{print $1 "\t" $2-1 "\t" $2}' snps_A > snpsA.bed


bedtools maskfasta -fi $GENOME \
 -bed snpsA.bed \
 -fo leptidea_ancestral${n}.fa \
 -mc A

rm snps_A
rm snpsA.bed

o=$((n+1))
grep -w --no-group-separator "T" $SNP_list > snps_T
awk '{print $1 "\t" $2-1 "\t" $2}' snps_T > snpsT.bed

bedtools maskfasta -fi leptidea_ancestral${n}.fa \
 -bed snpsT.bed \
 -fo leptidea_ancestral${o}.fa \
 -mc T

rm leptidea_ancestral${n}.fa
rm snps_T
rm snpsT.bed

n=$o
o=$((n+1))
grep -w --no-group-separator "C" $SNP_list > snps_C
awk '{print $1 "\t" $2-1 "\t" $2}' snps_C > snpsC.bed

bedtools maskfasta -fi leptidea_ancestral${n}.fa \
 -bed snpsC.bed \
 -fo leptidea_ancestral${o}.fa \
 -mc C

rm leptidea_ancestral${n}.fa
rm snps_C
rm snpsC.bed

n=$o
grep -w --no-group-separator "G" $SNP_list > snps_G
awk '{print $1 "\t" $2-1 "\t" $2}' snps_G > snpsG.bed

bedtools maskfasta -fi leptidea_ancestral${n}.fa \
 -bed snpsG.bed \
 -fo leptidea_ancestral_JuveSpecial_unMasked.fa \
 -mc G

rm leptidea_ancestral${n}.fa
rm snps_G
rm snpsG.bed
