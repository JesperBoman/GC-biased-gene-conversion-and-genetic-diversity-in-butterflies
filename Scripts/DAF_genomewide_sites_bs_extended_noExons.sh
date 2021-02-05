#!/bin/bash -l


SNP_list="swe_sin_DAF_type_comb_noExons"
num_rows=$(wc -l < "$SNP_list" | cut -f1 -d ' ')



for ((c=1; c<=1000; c++))
do
shuf -r -n $num_rows $SNP_list > SNP_list_sampled
awk '{if($8 == "REF"){if($7 == "Neutre"){print $6}} else{if($7 == "Neutre"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "Neutre" "\t" $0}' >> DAF_genomewide_bs_noExons
awk '{if($8 == "REF"){if($7 == "SW"){print $6}} else{if($7 == "SW"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "SW" "\t" $0}' >> DAF_genomewide_bs_noExons
awk '{if($8 == "REF"){if($7 == "WS"){print $6}} else{if($7 == "WS"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n  | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "WS" "\t" $0}' >> DAF_genomewide_bs_noExons
done

awk '{if(NR==1){n=1; print n "\t" $0} else if($1 =="Neutre"){n=n+1; print n "\t" $0} else{print n "\t" $0}}' DAF_genomewide_bs_noExons > DAF_genomewide_bs_noExons_num
rm DAF_genomewide_bs_noExons

echo "DAF_genomewide_bs_noExons finished"




SNP_list="swe_sin_DAF_type_comb_noExons_noCpG"
num_rows=$(wc -l < "$SNP_list" | cut -f1 -d ' ')



for ((c=1; c<=1000; c++))
do
shuf -r -n $num_rows $SNP_list > SNP_list_sampled
awk '{if($8 == "REF"){if($7 == "Neutre"){print $6}} else{if($7 == "Neutre"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "Neutre" "\t" $0}' >> DAF_genomewide_bs_noExons_noCpG
awk '{if($8 == "REF"){if($7 == "SW"){print $6}} else{if($7 == "SW"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "SW" "\t" $0}' >> DAF_genomewide_bs_noExons_noCpG
awk '{if($8 == "REF"){if($7 == "WS"){print $6}} else{if($7 == "WS"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n  | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "WS" "\t" $0}' >> DAF_genomewide_bs_noExons_noCpG
done

awk '{if(NR==1){n=1; print n "\t" $0} else if($1 =="Neutre"){n=n+1; print n "\t" $0} else{print n "\t" $0}}' DAF_genomewide_bs_noExons_noCpG > DAF_genomewide_bs_noExons_noCpG_num
rm DAF_genomewide_bs_noExons_noCpG

echo "DAF_genomewide_bs_noExons_noCpG finished"


SNP_list="swe_sin_DAF_type_comb_noExons_onlyCpG"
num_rows=$(wc -l < "$SNP_list" | cut -f1 -d ' ' )



for ((c=1; c<=1000; c++))
do
shuf -r -n $num_rows $SNP_list > SNP_list_sampled
awk '{if($8 == "REF"){if($7 == "Neutre"){print $6}} else{if($7 == "Neutre"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "Neutre" "\t" $0}' >> DAF_genomewide_bs_noExons_onlyCpG
awk '{if($8 == "REF"){if($7 == "SW"){print $6}} else{if($7 == "SW"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "SW" "\t" $0}' >> DAF_genomewide_bs_noExons_onlyCpG
awk '{if($8 == "REF"){if($7 == "WS"){print $6}} else{if($7 == "WS"){print$5}}}' SNP_list_sampled |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n  | uniq -c | awk -f ../DAF_gap_fill.awk | awk -f ../transpose.awk | awk '{if(FNR==1)print "WS" "\t" $0}' >> DAF_genomewide_bs_noExons_onlyCpG
done

awk '{if(NR==1){n=1; print n "\t" $0} else if($1 =="Neutre"){n=n+1; print n "\t" $0} else{print n "\t" $0}}' DAF_genomewide_bs_noExons_onlyCpG > DAF_genomewide_bs_noExons_onlyCpG_num
rm DAF_genomewide_bs_noExons_onlyCpG

echo "DAF_genomewide_bs_noExons_onlyCpG finished"
