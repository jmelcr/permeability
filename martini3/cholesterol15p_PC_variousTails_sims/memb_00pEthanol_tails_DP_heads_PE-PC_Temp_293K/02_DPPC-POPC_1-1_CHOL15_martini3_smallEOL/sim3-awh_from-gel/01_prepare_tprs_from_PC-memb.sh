#!/bin/bash

# in a previous simulation 
# I equilibrated a patch of a DPPC/POPC/15%chol membrane 
#
# Now I will use this equilibrated configuration to
# run a test whether using smaller dt will lead to smoother friction profile from AWH

source_dir="../../01_DPPC-CHOL15_martini3/sim1-awh_gel" 
mdpfname="../../../../martini_3_awh-md_293K.mdp"

init_tpr="${source_dir}/topol.tpr"
init_cpt="${source_dir}/state.cpt"
ln -s ${source_dir}/awh_ref_profile.xvg

# prep the sim
gmx grompp -f ${mdpfname} -p ../system_1mol.top -c ${init_tpr} -t ${init_cpt} -n ../index.ndx -v -o topol.tpr -maxwarn 1
# run a short test simulation!
tsp -L test_sim_${PWD} gmx mdrun -nsteps 1000 -awh awh_ref_profile.xvg
#tsp -L ${deffnm} gmx mdrun -deffnm ${deffnm} -cpi ${deffnm}.cpt -append -nsteps 1000
