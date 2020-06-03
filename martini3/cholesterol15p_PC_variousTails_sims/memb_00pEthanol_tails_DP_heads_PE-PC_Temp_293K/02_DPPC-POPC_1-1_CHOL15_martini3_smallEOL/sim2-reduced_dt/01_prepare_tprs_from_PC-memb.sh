#!/bin/bash

# in a previous simulation 
# I equilibrated a patch of a DPPC/POPC/15%chol membrane 
#
# Now I will use this equilibrated configuration to
# run a test whether using smaller dt will lead to smoother friction profile from AWH

source_dir="../sim1-awh"
mdpfname="martini_3_awh-md_293K_small-dt.mdp"

# change the dt in the basic mdp file
cp    ../../../../martini_3_awh-md_293K.mdp $mdpfname
sed -e 's/dt                       = 0.02/dt                       = 0.002/g' -i $mdpfname

init_tpr="${source_dir}/topol.tpr"
init_cpt="${source_dir}/state.cpt"
ln -s ${source_dir}/awh_ref_profile.xvg

# prep the sim
gmx grompp -f ${mdpfname} -p ../system_1mol.top -c ${init_tpr} -t ${init_cpt} -n ../index.ndx -v -o topol.tpr -maxwarn 1
# run a short test simulation!
tsp -L test_sim_${PWD} gmx mdrun -nsteps 1000 -awh awh_ref_profile.xvg
#tsp -L ${deffnm} gmx mdrun -deffnm ${deffnm} -cpi ${deffnm}.cpt -append -nsteps 1000
