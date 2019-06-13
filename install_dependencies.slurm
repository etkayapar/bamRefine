#!/bin/bash -l
#SBATCH --partition=bonobo                    # Run a single task  
#SBATCH --ntasks=1                  # Run a single task        
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --output=/mnt/NEOGENE1/home/etka/tmp/bamrefine.out
#SBATCH  --error=/mnt/NEOGENE1/home/etka/tmp/bamrefine.err


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"


pip3 install --user cython
pip3 install --user pysam



exit

