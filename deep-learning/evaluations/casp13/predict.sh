#!/bin/bash

while read id; do
    echo $id
    python3 ../../src/predict-contacts.py -w ../../experiments/contact256.hdf5 -p ../../../data/casp13/features/$id.pkl -o ./out-rr/$id.contact.rr
    python3 ../../src/predict-distances.py -w ../../experiments/distance256.hdf5 -p ../../../data/casp13/features/$id.pkl -o ./out-rr/$id.4bydistance.rr
done < '../../../data/casp13/casp13-public-targets.lst'
