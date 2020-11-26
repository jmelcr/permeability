#!/bin/bash
# simple analysis - calculate AWH profile from simulations

rd=$PWD

for file in $(find ./D* -name  "topol.tpr")
do
   wd=`dirname $file`
   cd $wd
   [ ! -f ener_all.edr ] && gmx eneconv -f ener.*edr -o ener_all.edr -noerror
   gmx awh -f ener_all.edr -fric -skip 5000 -kt -more
   #[ ! -f traj_cat_all_pbc.xtc ] && [ ! -f traj_cat_all.xtc ] && gmx trjcat -f traj*.xtc -o traj_cat_all.xtc
   #[ ! -f traj_cat_all_pbc.xtc ] && echo System | gmx trjconv -f traj_cat_all.xtc -o traj_cat_all_pbc.xtc -pbc mol && rm traj_cat_all.xtc 
   #bash ${rd}/surface_excess_calculate.sh traj_cat_all_pbc.xtc
   #[ ! -f boxx.xvg ] && echo Box-X | gmx energy -f ener_all.edr -o boxx.xvg
   #[ ! -f lipidator.out ] && tsp -L RESIS_${wd} bash -c "${rd}/resis traj_cat_all_pbc.xtc ${rd}/martini.ldx ../../lipidator.ndx > lipidator.out"
   #mv CHECK unCHECK 
   cd $rd
done

