#u!/usr/bin/awk -f
{
if(FNR==NR){
	a[$2]=$1;
} 

else if(NR!=FNR){ 
	if($2 in a){a[$2]+=$1} 
	else {a[$2]=$1}
}
}

END {
	for (key in a) print a[key] "\t"  key
}
