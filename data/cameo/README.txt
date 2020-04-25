
wget https://yanglab.nankai.edu.cn/trRosetta/benchmark/cameo_a3m.tar.bz2 --no-check-certificate

tar -xjvf cameo_a3m.tar.bz2

# Features
./gen-feats.sh

# Distances
while read id; do     echo $id;     python3 ../../scripts/pdb2dist.py -i ./chains/$id.pdb -b ./distance/$id-cb -a ./distance/$id-ca -p; done < 'cameo-hard.lst' &> ./logs-data-prep/cameo-dist-gen.log &


