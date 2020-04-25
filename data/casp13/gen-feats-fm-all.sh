#!/bin/bash

while read id; do
    echo $id
    # Generate dmap
    python3 ../../scripts/pdb2dist.py -i ./chains/$id.pdb -b ./distance/$id-cb -a ./distance/$id-ca -p
    if [ -f ./features/$id.pkl ]; then
        echo "already done!!"
    else
        echo "running.."
        head -n 2 ./a3m_from_trRosetta/casp13_a3m/$id.a3m > ./fasta/$id.fasta
        # Generate feature
        python3 ../../scripts/make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m &> ./raw-features/$id.log
    fi
done < 'casp13-public-targets.lst'
