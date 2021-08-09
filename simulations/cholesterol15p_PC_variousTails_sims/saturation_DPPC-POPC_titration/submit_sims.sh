#!/bin/bash

rootdir=`pwd`

for phase in gel fluid
do

for d in `find -type d -name "*sim1_awh*${phase}"`
do
  dir=$d
  cd $rootdir

  if `cd ${dir}`
  then
    cd ${dir}
    cp $rootdir/submit.sh ./submit.sh
    simname=`echo ${dir} | sed -r 's/\//_/g'`
    echo $simname 
    sed `echo s/__NN__/${simname}/` submit.sh -i
    sed `echo s/__PHASE__/${phase}/` submit.sh -i

    # if there is file ".submitted .NOT. -> submit the job and if successful -> create the file
    if [ ! -s jobs.submitted ] 
    then
       sbatchoutput=`sbatch submit.sh` && touch jobs.submitted
       jobid=`echo $sbatchoutput | cut -f4 -d " "`
       # store the jobsnames locally in sim folder
       echo $jobid >> jobs.submitted
       # store the jobsnames globally in main folder
       echo "# sim "${simname} >> $rootdir/jobs.submitted
       echo $jobid >> $rootdir/jobs.submitted
    fi

    for j in `seq 4`
    do
      jobid=`sbatch --dependency=singleton submit.sh `
      # store the jobsnames locally in sim folder
      echo $jobid >> jobs.submitted
      # store the jobsnames globally in main folder
      echo $jobid >> $rootdir/jobs.submitted
    done 

    cd $rootdir

  else
    "Could not step into the directory "$dir
  fi

done   

done   
