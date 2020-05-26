#! /bin/bash
# simple script, which prepares a set of simulations

for d in PC_*
do
rd=`pwd`
cd $d
echo $PWD
simdirname=sim1_awh_fluid
mkdir -p $simdirname && cd $simdirname || exit 66
tsp gmx grompp -c ../../confout_fluid.gro  -f ${rd}/../../../martini_3_awh-md_293K.mdp -n ${rd}/index.ndx -p ../system_1mol.top -maxwarn 1
tsp -L $PWD gmx mdrun -v -awh ../../awh_ref_profile_fluid.xvg -cpi state.cpt -append -nsteps 1000
cd $rd
done

