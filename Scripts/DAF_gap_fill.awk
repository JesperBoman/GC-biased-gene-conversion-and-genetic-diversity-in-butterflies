#!/usr/bin/awk -f

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# awk -f DAF_gap_fill.awk <DAF spectrum without 0's> > <Derived allele spectrum with 0's>
# ==========================================================================================================
# Add zeros in 
# ==========================================================================================================
# Modified from https://www.biostars.org/p/152592/
# ==========================================================================================================
# Jesper Boman                      15 nov 2020
# ==========================================================================================================
# Example:
# Use in pipe like this (start with a file of counts for each allele at a position):
# awk '{if($8 == "REF"){if($7 == "Neutre"){print $6}} else{if($7 == "Neutre"){print$5}}}' ${pop}_DAF_type_comb_noExons_noCpG |  sed -E 's/[A-Z]:([0-9]+)/\1/' | sort -n | uniq -c | awk -f ../DAF_gap_fill.awk 
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #

n= #sample size

BEGIN{X=1;}
{
 if($2==X){print $1 "\t" $2;} 
 else{ while (X<$2){print "0" "\t" X;  X++}
	print $1 "\t" $2;
	}
 X++;
}
END{
while (X<=n){print "0" "\t" X;  X++}
}
