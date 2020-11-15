#!/usr/bin/awk -f

{

split($5, a, ":"); 
split($6, b, ":");

##REF is ancestral
if(a[1] == $7){
        ##Assign Strong -> Weak
        if (($7~/C/ || $7~/G/) && ($6~/A/ || $6~/T/))
                {print $0, "\t", "SW", "\t", "REF" }

        ##Assign Weak -> Strong
        else if (($7~/A/ || $7~/T/) && ($6~/C/ || $6~/G/))
                {print $0, "\t", "WS", "\t", "REF" }

        ##Assign Neutral
        else
                {print $0, "\t", "Neutre", "\t", "REF" }
}

##ALT is ancestal
else if(b[1] == $7){
        ##Assign Strong -> Weak
        if (($7~/C/ || $7~/G/) && ($5~/A/ || $5~/T/))
                {print $0, "\t", "SW", "\t", "ALT" }

        ##Assign Weak -> Strong
        else if (($7~/A/ || $7~/T/) && ($5~/C/ || $5~/G/))
                {print $0, "\t", "WS", "\t", "ALT" }

        ##Assign Neutral
        else
                {print $0, "\t", "Neutre", "\t", "ALT" }
}

else {next;

#OPTIONAL assign multiple polarizations - directions will be unsure if ancestral allele is not known for certain
        ##Assign multiple polarizations
        ##First allele
       # if (($7~/C/ || $7~/G/) && ($5~/A/ || $5~/T/))
       #         {first_allele="SW"}
       # else if (($7~/A/ || $7~/T/) && ($5~/C/ || $5~/G/))
       #         {first_allele="WS"}
       # else
       #        {first_allele="Neutre"}

        ##Second allele
        #if (($7~/C/ || $7~/G/) && ($6~/A/ || $6~/T/))
        #        {second_allele="SW"}
        #else if (($7~/A/ || $7~/T/) && ($6~/C/ || $6~/G/))
        #        {second_allele="WS"}
        #else
        #        {second_allele="Neutre"}
        
        #print $0, "\t", "Other:" first_allele ":" second_allele
}
}
