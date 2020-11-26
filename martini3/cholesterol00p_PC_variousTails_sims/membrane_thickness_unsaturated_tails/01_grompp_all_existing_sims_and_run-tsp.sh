#! /bin/bash
# simple script, which prepares TPR files 
# from an existing set of specific simualations 
# using AWH to calculate permeability
#

rd=`pwd`

mdp_fname="martini_3_awh-md_293K.mdp"
init_gro_fname="confout.gro"

# for each existing simulations configuration
for init_gro in $(find -name $init_gro_fname)
do
wd=$(dirname $init_gro)
cd $wd

# use the awh profile as an initial guess for quick convergence
[ awh_init.xvg -s ] || ln -s awh_t2e+06.xvg awh_init.xvg

# correct the permeating particle name in index file
sed -e "s/EOL/particle/g" index.ndx -i

# correct the permeating particle name and WF→W in topology file
sed -e "s/EOL/EOLT/g" -e "s/WF/W /g" system_1mol.ndx -i

# grompp the sim, set final time to 1 us and try it to run it using TSP if succesful
gmx grompp -c ${init_gro_fname}  -f ${rd}/../../${mdp_fname} -n index.ndx -p ../system_1mol.top && gmx convert-tpr -until 1000000 -o topol.tpr &&  tsp -N 3 gmx mdrun -awh awh_init.xvg -v  # nsteps 10000 

cd ${rd}
done

