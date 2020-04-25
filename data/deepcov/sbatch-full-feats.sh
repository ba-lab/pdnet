#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition gpu3,gpu4
#SBATCH --nodes=1
#SBATCH --ntasks=1         # leave at '1' unless using a MPI code
#SBATCH --cpus-per-task=1  # cores per task
#SBATCH --mem-per-cpu=10G  # memory per core (default is 1GB/core)
#SBATCH --time 1-23:00     # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=general-gpu  # investors will replace this with their account name
#SBATCH --gres gpu:1
##SBATCH --exclusive

## labels and outputs
#SBATCH --job-name=usagpu
#SBATCH --output=_%j.out  # %j is the unique jobID

## notifications
#SBATCH --mail-user=adhikarib@umsl.edu  # email address for notifications
#SBATCH --mail-type=FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# load modules then list loaded modules
source /group/prayog/virt-env3.5/bin/activate
module load cudnn/cudnn-7.4.2-cuda-10.0.130

# Science goes here:

export HDF5_USE_FILE_LOCKING=FALSE

for id in `sort -R $1 | cut -d' ' -f1`; do
	#echo $id
	fseq=`tail -n 1 ./fasta/$id.fasta`
    aseq=`head -n 1 ./aln.original.deepcov.paper/$id.aln`
    
    if [[ -f "./features/$id.pkl" ]]; then
        #echo " already done!"
        continue
    fi

    echo $id

    if [[ $fseq == $aseq ]]; then
        echo "use existing alignments.."
        cp ./trRosetta-pssm/$id.pssm.npy ./raw-features/$id/
        python3 ../../scripts/make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id/ -p ./features/$id.pkl -a ./aln.original.deepcov.paper/$id.aln &> ./raw-features/$id.log
    else
        echo "generate new alignments.."
        #Should automatically pick up existing alignment a3mfile and make new PSSMs
        python3 ../../scripts/make-features.py -f ./fasta/$id.fasta -j ./raw-features/$id/ -p ./features/$id.pkl &> ./raw-features/$id.log
    fi
done

echo "### Finished at: $(date) ###"


