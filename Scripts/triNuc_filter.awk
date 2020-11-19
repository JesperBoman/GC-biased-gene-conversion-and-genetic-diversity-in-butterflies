#u!/usr/bin/awk -f
#Define genome as input variable

{
cmd = "samtools faidx " genome " "$1":"$2-1"-"$2+1" | grep -v '>'"
cmd | getline triNuc 

if(triNuc~/CG/ || triNuc~/TG/ || triNuc~/CA/ || triNuc~/NG/ || triNuc~/TN/ || triNuc~/CN/ || triNuc~/NA/){
        print $0 "\t" "CpG_ans";}
else{print $0 "\t" "Other"}
close(cmd);
}
