#!/bin/bash

while IFS=' ' read -r domain target
do
	echo ""
	echo $target
    # Long-range only
    perl ../../../scripts/coneva.pl -all -pdb ../../../data/casp13/chains/$domain.pdb -rr ./out-rr/$target.contact.rr
    # Medium- and long-range
    #perl ../../../scripts/coneva.pl -smin 12 -all -pdb ../../../data/casp13/chains/$domain.pdb -rr ./predictions-raptorx/$target.rr
done < "domains.lst"
