#!/bin/bash

# this bash script invokes analysis of each individual trajectory in the folder structure
# calculating the following properties:
#   - Free energy dG profile from AWH
#   - density profile (solvent) -> bilayer thickness
#   - box-x,y fluctuation -> compressibility modulus
#   - LATER: poreformation Free energy from additional MetaD run (new Gromacs compilation req'd)
#   - MAYBE LATER: density profile (ethanol, water, solvent) -> surface excess of eth.
#   - SKIPPED HERE: bending/tilt moduli (resis analysis via lipidator)
#      - & APL will be calcualted along the way too (mentioned in lipidator.out)
#
# These properties are not recalculated, if their output files exist. 
#
# The script is supposed to be run on Cartesius in the environment of jmelcr
# assumes installed Task Spooler ("ts" or "tsp") to batch process tasks
# and Gromacs tools (v2019+)
#
# python script/library get_conc_ion_bulk.py
# 
# All simulations are assumed to have the very same rigid file and folder structure format
#
#------------------------------------------------------------
# Made by J.Melcr,  Last edit 2020/03/09
#------------------------------------------------------------
#

export PATH=$HOME/.local/bin/:$PATH

module purge
module load 2019
module load GROMACS/2019.3-foss-2018b

scriptdir=`dirname $0`
rd=$PWD

# GROMACS command definition
gmx="gmx"

# file-naming conventions used in this analysis script
traj_file_name="traj_comp.xtc"
tpr_file_name="topol.tpr"
dens_file_name="density_solvent.xvg"
#top="../system.top"


#for file in `find -name "topol.tpr"`
for d in `find -type d -name "*sim1_awh*"`
do
  dir=$d
  cd $rd

  if `cd ${dir}`
  then
    cd ${dir}
    # workaround: now $d, $dir, $wd have the same value 
    # because the code is stitched from different sources. Huh. 
    wd=${dir}  
  else
    echo "Could not step into the directory "$dir
    cd $rd
    continue
  fi

   # is the simulation finished?
   if ! [ -s confout.gro ] 
   then
       echo "This is likely an unfinished simulation and as such we will skip it."
       echo "Place:" $wd
       echo $wd >> ${rd}/skipped.sims.list
       # skip to the next iteration of the most inner iterator (here for loop)
       cd $rd
       continue
   fi

   # is it a simulation folder? (shall contain topol.tpr)
   if ! [ -s $tpr_file_name ] 
   then
       echo "We really need " $tpr_file_name " , but we can't find it!"
       echo "Place:" $wd
       echo $wd >> ${rd}/problematic.sims.list
       # skip to the next iteration of the most inner iterator (here for loop)
       cd $rd
       continue
   fi

   # Assuming x=y for sim-box size!!
   [ ! -f boxx.xvg ] && ts -L $wd bash -c "echo Box-X | $gmx energy -f ener.edr -o boxx.xvg"

   # get densities of (ethanol, water, solvent) centered around POPC
   if ! [ -s $dens_file_name ] 
   then
       ts -L $wd bash -c "echo Membrane Solvent | $gmx density -sl 400 -dens number -ng 1 -f $traj_file_name -center -symm -relative -o $dens_file_name -n ../index.ndx"
   fi
   #[ ! -s lipidator.out ] && ts -L ${wd} bash -c "${rd}/resis ${traj_file_name} ${rd}/martini.ldx ../../lipidator.ndx > lipidator.out"

    # get dG profile from AWH
    #[ ! -s awh_t2e+06.xvg ] && ts gmx awh -skip 1000 -more -kt -fric
    ts gmx awh -skip 1000 -more -kt -fric
   cd $rd
done


