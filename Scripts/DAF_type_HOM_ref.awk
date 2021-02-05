#!/usr/bin/awk -f

{ 
##Assign Strong -> Weak
if (($5~/C/ || $5~/G/) && ($6~/A/ || $6~/T/))
	{print $0, "	", "SW"}

##Assign Weak -> Strong
else if (($5~/A/ || $5~/T/) && ($6~/C/ || $6~/G/))
	{print $0, "	", "WS"}

##Assign Neutral
else
	{print $0, "	", "Neutre"}
}
