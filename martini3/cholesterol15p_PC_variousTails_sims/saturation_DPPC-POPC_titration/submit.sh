#!/bin/bash

#### Jobname and output file
#SBATCH --job-name=CHOL15__NN__
#SBATCH --output=%J.out
####
#SBATCH -N 1 --ntasks-per-node=24
####
#SBATCH -t 24:00:00
#SBATCH --mail-type=FAIL  
#SBATCH --mail-user=j.melcr@rug.nl
#SBATCH --dependency=singleton

# Set flags for mdrun (-cpi is set by the script when needed, no need to to include it here)
FLAGS='-maxh 23.5 -cpi state.cpt -append -awh ../../../awh_ref___PHASE__.xvg'
# Set working directory
echo "workfolder:" `pwd`
echo "On host:" `hostname`
echo

module purge
module load 2019
module load GROMACS/2019.3-intel-2018b
module list

srun gmx_mpi mdrun $FLAGS  >> mdrun.outlog 2>&1 

