#! /bin/bash
# simple script, which prepares TPR files 
# from an existing set of specific simualations 
# using AWH to calculate permeability
#

rd=`pwd`

mdp_fname="martini_3_awh-md_293K.mdp"
init_gro_fname="confout.gro"
top_final="system.top"

ln -s ../saturation_DPPC-POPC_titration/system_template.top ./

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
sed -e "s/EOL /EOLS/g" -e "s/WF/W /g" ../system_1mol.top -i

# combine the correct header from topology template 
# with the correct system-molecular definition from the system_1mol.top
grep -e " system " -B 10 ${rd}/system_template.top > ../${top_final}
echo  "Membrane simulation of varying tail length, both tails monounsaturated" >> ../${top_final}
echo  " " >> ../${top_final}
grep -e " molecules " -A 30 ../system_1mol.top >> ../${top_final}

# grompp the sim, set final time to 1 us and try it to run it using TSP if succesful
gmx grompp -c ${init_gro_fname}  -f ${rd}/../../${mdp_fname} -n index.ndx \
    -p ../${top_final} \
    -maxwarn 1  && \
   gmx convert-tpr -until 1000000 -o topol.tpr && \
   tsp -N 3 gmx mdrun -awh awh_init.xvg -v  # nsteps 10000 

cd ${rd}
done

