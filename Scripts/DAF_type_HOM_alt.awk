#!/usr/bin/awk -f

{ 
##Assign Strong -> Weak
if (($6~/C/ || $6~/G/) && ($5~/A/ || $5~/T/))
	{print $0, "	", "SW"}

##Assign Weak -> Strong
else if (($6~/A/ || $6~/T/) && ($5~/C/ || $5~/G/))
	{print $0, "	", "WS"}

else
	{print $0, "	", "Neutre"}
}
