#!/bin/bash

# First learn the feature to distance mapping
python3 ../src/train.py -t distance -w distance.hdf5 -n -1 -c 128 -e 64 -d 128 -f 64 -p ../../data/ -v 0 -o 1 &> distance.log &

# Run 8 epochs of distance, contact at cropsize 256
python3 ../src/train.py -t distance -w distance.hdf5 -n -1 -c 256 -e 16 -d 128 -f 64 -p ../../data/ -v 0 -o 1 &> distance256.log &
python3 ../src/train.py -t contact -w contact.hdf5 -n -1 -c 256 -e 8 -d 128 -f 64 -p ../../data/ -v 0 -o 1 &> contact256.log &

# Run 4 epochs of binned (loading distance weights)
python3 ../src/train.py -t binned -w binned256.hdf5 -n -1 -c 128 -e 8 -d 128 -f 64 -p ../../data/ -v 0 -o 1 &> binned256.log &
