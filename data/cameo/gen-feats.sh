#!/bin/bash

for id in `sort -R ../cameo-hard.lst | cut -d' ' -f1`; do
    echo $id
    if [ -f ./features/$id.pkl ]; then
        echo "already done!!"
    else
        echo "running.."
        echo ">$id" > ./fasta/$id.fasta
        head -n 2 ./cameo_a3m/$id.a3m | tail -n 1 >> ./fasta/$id.fasta
        # Generate dmap
        python3 ../../scripts/pdb2dist.py -i ./chains/$id.pdb -b ./distance/$id-cb -a ./distance/$id-ca -p
        # Generate feature
        python3 ../../scripts/make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./cameo_a3m/$id.a3m &> ./raw-features/$id.log
    fi
done
