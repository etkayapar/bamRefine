#!/bin/bash -l
#SBATCH --partition=bonobo                    # Run a single task  
#SBATCH --ntasks=8                  # Run a single task        
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --output=/mnt/NEOGENE1/home/etka/tmp/cythonize/bamrefine.out
#SBATCH  --error=/mnt/NEOGENE1/home/etka/tmp/cythonize/bamrefine.err


echo "SLURM_NODELIST $SLURM_NODELIST"
echo "NUMBER OF CORES $SLURM_NTASKS"


# py='/mnt/NEOGENE1/home/etka/miniconda3/bin/python3.7'

cd $dir

time python3 main_wrapper.py -i deneme2.bam -o deneme2_filtered.bam -p 16

exit

