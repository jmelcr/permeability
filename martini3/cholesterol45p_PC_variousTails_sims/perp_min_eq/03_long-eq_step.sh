#!/bin/bash

source_mdp="../../martini_3_awh-md_293K.mdp"

# Equilibrated structure exists from prev step.
#  + initialize initial configuration
dt=020
eqno="eq_dt_${dt}"
conf_init=../${eqno}.gro


simdir="long_eq"
mkdir -p ${simdir} 
cd ${simdir}


# LONG EQ with pressure coupling

eqno="eq-long_dt_${dt}"
mdpfname="martini_3_awh-${eqno}_293K.mdp"

# create equilibration mdp file
sed -e "s/awh1-user-data = yes/awh1-user-data = no/g" \
    -e 's/Pcoupl                   = parrinello-rahman/Pcoupl                   = berendsen/g' \
    -e "s/dt                       = 0.02/dt                       = 0.${dt}/g" \
    ../${source_mdp}  > ${mdpfname}

gmx grompp -f ${mdpfname} -c ${conf_init} -p ../system_1mol.top -n ../index.ndx 


gmx mdrun -v -nsteps 10000 

