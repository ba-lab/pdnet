#!/bin/bash

for id in `sort -R $1 | cut -d' ' -f1`; do
	#echo $id
	fseq=`tail -n 1 ./fasta/$id.fasta`
    aseq=`head -n 1 ./aln.deepcov/${id}0.aln`
    
    if [[ -f "./features/$id.pkl" ]]; then
        #echo " already done!"
        continue
    fi

    echo $id

    if [[ $fseq == $aseq ]]; then
        cp ./trRosetta-pssm/$id.pssm.npy ./raw-features/$id/
    fi
    #Should automatically pick up existing alignment a3mfile and make new PSSMs
    python3 ../../scripts/make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id/ -p ./features/$id.pkl &> ./raw-features/$id.log
done

echo "### Finished at: $(date) ###"


