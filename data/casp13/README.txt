------------------------------------------------------------------------
# Prepare distance files from PDB
------------------------------------------------------------------------
# Download chains
cd chains
wget http://predictioncenter.org/download_area/targets/casp13.targets.T-D.4public.tar.gz
tar zxvf casp13.targets.T-D.4public.tar.gz
wget http://predictioncenter.org/download_area/targets/casp13.targets.T.4public.tar.gz
tar zxvf casp13.targets.T.4public.tar.gz

while read id; do
    echo $id
    python3 pdb2dist.py -i ./casp13_trRosetta/chains/$id.pdb -b ./casp13_trRosetta/distance/$id-cb -a ./casp13_trRosetta/distance/$id-ca -p
done < 'casp13-trrosetta.lst'

# Crosscheck fasta
while read id flag; do
    echo $id
    cat ./fasta/$id.fasta
    cat ./a3m_from_trRosetta/fasta/$id.fasta
done < 'casp13-public.lst'

------------------------------------------------------------------------
Generate features (in Lewis server)
------------------------------------------------------------------------
# Old approach (ignore)
while read id; do
    echo $id
    if [ -f ./a3m_from_trRosetta/casp13_a3m/$id.a3m ]; then
        sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
    else
        sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
    fi
done < 'casp13-public.lst'

# Generate features (using trRosetta a3m)
id="T0950"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh
id="T0953s2"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0957s1"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0957s2"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0960"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0963"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0968s1"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"
id="T0968s2"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl -m ./a3m_from_trRosetta/casp13_a3m/$id.a3m" "./raw-features/$id.log"

# Generate features (using RaptorX a3m)
id="T0953s1"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_ure3/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T0954"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_ure5/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T0958"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_uce0/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T0966"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_uce0/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1005"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_uce3/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1008"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_uce0/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1009"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-k40.sh "../make-features.py -m ./a3m_from_raptorx/CASP13MyDM-Features/${id}_contact/${id}_ure5/$id.a3m -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"

# Generate features (without aln)
id="T0951"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T0955"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1003"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1011"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
id="T1016"
sbatch /storage/hpc/group/prayog/ba-lab/sbatch-any-gpu.sh "../make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id -p ./features/$id.pkl" "./raw-features/$id.log"
