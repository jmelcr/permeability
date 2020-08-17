#!/bin/bash

#### Jobname and output file
#SBATCH --job-name=awhana
#SBATCH --output=%J.out
#SBATCH -n 1 -c 24
#SBATCH -t 1-0:00:00
#SBATCH --mail-type=FAIL  
#SBATCH --mail-user=j.melcr@rug.nl
#SBATCH --dependency=singleton


module purge
module load 2019
module load GROMACS/2019.3-foss-2018b

# path to task spooler
export PATH=$HOME/.local/bin/:$PATH

echo "workfolder:" `pwd`
echo "On host:" `hostname`
echo
echo "Env variables:"
echo $env
echo
echo "Path is:"
echo $PATH
echo
echo "OMP number threads is:"
echo ${OMP_NUM_THREADS}
echo
echo

# project folder and on scratch folder allow for fast parallel I/O ops
# Task Spooler is used for parallel processing of individual tasks
echo "setting task spooler to use ${OMP_NUM_THREADS} parallel tasks for analysis"
ts -S ${OMP_NUM_THREADS}
echo

bash analyze.sh

# this will ensure that SLURM will wait till all jobs in TaskSpooler end
# ts -f  ->  do not fall into background 
# ts -N ${OMP_NUM_THREADS}  ->  Use all threads (previously set as max_threads) 
#                               => will execute only after all previous jobs finish
ts -N ${OMP_NUM_THREADS} -f echo "Finished!"

