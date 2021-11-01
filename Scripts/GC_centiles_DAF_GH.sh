#!/bin/bash -l

# Usage: Bash GC_centiles_DAF_GH.sh
# ==========================================================================================================
# A pipeline to create derived allele frequency spectra for a range (100 windows) of GC values based on local GC content
# The pipeline also synthesises other important information for the GC centiles such as GC content, nucleotide diversity (pi), 
# CDS density, the number of unmasked bases (L) etc.
# ==========================================================================================================
# Jesper Boman                      27 aug 2020
# ==========================================================================================================
# Extra info:
# See gBGC_main_pipeline.txt for information on how dependent files were created.
# 
# Script dependencies:
# DAF_gap_fill.awk
# transpose.awk 
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #


population="swe_sin"
skipn=$(wc -l ${population}_noExons_GC_scafWind | cut -f1 -d ' ' | awk '{print $1%100}')
SNPs_per_C=$(wc -l ${population}_noExons_GC_scafWind | cut -f1 -d ' ' | awk '{print int($1/100)}')

awk -v sk="$skipn" 'NR>sk{print $0}' ${population}_noExons_GC_scafWind | awk -v SNPs="$SNPs_per_C" -v pop="$population" 'NR==1{i=1}{print >> pop"_noExons_GC"i; if(NR%SNPs==0){i=i+1}}'

awk 'FNR==NR{a[$11]; next} !($1 in a)' ${population}_noExons_GC_scafWind ${population}_mm1_1kb_pi_GC_noN_noExons > variant_windows_with_no_pol_sites 
awk 'FNR==NR{a[$1]; next} !($1 in a)' ${population}_mm1_1kb_pi_GC_noN_noExons ../leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_1kb_windows_GC1kb_callable_bp > invariant_windows


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




#Different scritps for pi, L and GC
awk 'FNR==NR{a[$11]; next} $1 in a' ${population}_noExons_GC${c} ${population}_mm1_1kb_pi_GC_noN_noExons | awk '{sum_of_sites_pi+=$2; L+=$4; n+=$6; gc+=$7} END{print n "\t" sum_of_sites_pi "\t"  sum_of_sites_pi/L "\t" L "\t" gc}' > pi_n_L_pol_anchor

pi_n_L_pol_anchor=$(cat pi_n_L_pol_anchor)

awk -v headGC="$headGC" -v tailGC="$tailGC" 'NR>1{if($1 !~ /mtANA/ && $5 >= headGC && $5 <= tailGC){pi_sum+=$2; L_sum+=$4; n+=$6; gc+=$7}} END{print n "\t" pi_sum "\t" pi_sum/L_sum "\t" L_sum "\t" gc}' variant_windows_with_no_pol_sites > pi_n_L_unanchor

pi_n_L_unanchor=$(cat pi_n_L_unanchor)
 

awk -v headGC="$headGC" -v tailGC="$tailGC" 'NR>1{if($1 !~ /mtANA/  && $3 >= headGC && $3 <= tailGC){L_sum+=$2; gc+=$4}} END{print L_sum "\t" gc}' invariant_windows > invar_length_temp

invar_L=$(cat invar_length_temp)


#CDS
awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}; FNR==NR{a[$1, ceil($2/1000)]; next}FNR!=NR{if($3=="CDS"){b[$1, ceil($4/1000)]=($5-$4)}} END{for (wind in a){sum+=b[wind]; c++};print sum "\t" c}' ${population}_noExons_GC${c} ../leptidea_sinapis_rc1.gff > CDS_density_pol_temp

cds_pol=$(cat CDS_density_pol_temp)


awk -v headGC="$headGC" -v tailGC="$tailGC" 'NR>1{if($1 !~ /mtANA/ && $3 >= headGC && $3 <= tailGC){split($1, a, "-"); print a[1] "\t" a[2]}}' ../leptidea_ancestral_JuveSpecial_Exons_repeats_HardMasked_1kb_windows_GC1kb_callable_bp > all_windows_temp 
awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}; FNR==NR{a[$1, $2]; next} FNR!=NR{if($3=="CDS"){b[$1, ceil($4/1000)]=($5-$4)}} END{for (wind in a){sum+=b[wind]; c++};print sum "\t" c}'  all_windows_temp ../leptidea_sinapis_rc1.gff > CDS_density_all_temp

cds_all=$(cat CDS_density_all_temp) 



#Creating DAF spectra and synthesizing information

awk '{if($8 == "REF"){if($7 == "Neutre"){print $6}} else{if($7 == "Neutre"){print$5}}}' ${population}_noExons_GC${c} |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c > temp_C
awk -f ../DAF_gap_fill.awk temp_C | awk -f ../transpose.awk | awk -v cds_all="$cds_all" -v cds_pol="$cds_pol" -v invar_L="$invar_L" -v pi_n_L_unanchor="$pi_n_L_unanchor" -v pi_n_L_pol_anchor="$pi_n_L_pol_anchor" -v var="$c" \
 '{if(FNR==1){split(pi_n_L_unanchor, a, "\t"); split(pi_n_L_pol_anchor, b, "\t"); split(invar_L, inv, "\t"); print "C_" var "\t" (a[5]+b[5]+inv[2])/(a[4]+b[4]+inv[1])  "\t" invar_L "\t" pi_n_L_unanchor "\t" pi_n_L_pol_anchor "\t" (a[2]+b[2])/(a[4]+b[4]+inv[1]) "\t" cds_pol "\t" cds_all "\t" "Neutre" "\t" $0}}' >> GC_centiles_DAF_noExons_data

awk '{if($8 == "REF"){if($7 == "SW"){print $6}} else{if($7 == "SW"){print$5}}}' ${population}_noExons_GC${c} |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c > temp_C
awk -f ../DAF_gap_fill.awk temp_C | awk -f ../transpose.awk | awk -v cds_all="$cds_all" -v cds_pol="$cds_pol" -v invar_L="$invar_L" -v pi_n_L_unanchor="$pi_n_L_unanchor" -v pi_n_L_pol_anchor="$pi_n_L_pol_anchor" -v var="$c" \
 '{if(FNR==1){split(pi_n_L_unanchor, a, "\t"); split(pi_n_L_pol_anchor, b, "\t"); split(invar_L, inv, "\t"); print "C_" var "\t" (a[5]+b[5]+inv[2])/(a[4]+b[4]+inv[1]) "\t" invar_L "\t" pi_n_L_unanchor "\t" pi_n_L_pol_anchor "\t" (a[2]+b[2])/(a[4]+b[4]+inv[1]) "\t" cds_pol "\t" cds_all "\t" "SW" "\t" $0}}' >> GC_centiles_DAF_noExons_data

awk '{if($8 == "REF"){if($7 == "WS"){print $6}} else{if($7 == "WS"){print$5}}}' ${population}_noExons_GC${c} |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n  | uniq -c > temp_C
awk -f ../DAF_gap_fill.awk temp_C | awk -f ../transpose.awk | awk -v cds_all="$cds_all" -v cds_pol="$cds_pol" -v invar_L="$invar_L" -v pi_n_L_unanchor="$pi_n_L_unanchor" -v pi_n_L_pol_anchor="$pi_n_L_pol_anchor" -v var="$c" \
 '{if(FNR==1){split(pi_n_L_unanchor, a, "\t"); split(pi_n_L_pol_anchor, b, "\t"); split(invar_L, inv, "\t"); print "C_" var "\t" (a[5]+b[5]+inv[2])/(a[4]+b[4]+inv[1]) "\t" invar_L "\t" pi_n_L_unanchor "\t" pi_n_L_pol_anchor "\t" (a[2]+b[2])/(a[4]+b[4]+inv[1]) "\t" cds_pol "\t" cds_all "\t" "WS" "\t" $0}}' >> GC_centiles_DAF_noExons_data


done