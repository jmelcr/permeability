#!/bin/bash

source_mdp="../../martini_3_awh-md_293K.mdp"

# Minimized structure exists from prev step.
#  + initialize initial configuration
conf_init="min.gro"


# EQ without pressure coupling

for dt in 001 002 #003 004 005 010 015 020
do

eqno="eq_noP_dt_${dt}"
mdpfname="martini_3_awh-${eqno}_293K.mdp"

# primitive chekcpointing
if ! [ -s ${eqno}.xtc ] 
then
	# create equilibration mdp file
	#sed -e 's/Pcoupl                   = parrinello-rahman/Pcoupl                   = berendsen/g' \
	sed -e "s/awh1-user-data = yes/awh1-user-data = no/g" \
	    -e 's/Pcoupl                   = parrinello-rahman/Pcoupl                   = no/g' \
	    -e "s/dt                       = 0.02/dt                       = 0.${dt}/g" \
	    -e "s/awh1-user-data = yes/awh1-user-data = no/g" \
	    ${source_mdp}  > ${mdpfname}

	gmx grompp -f ${mdpfname} -c ${conf_init} -p system_1mol.top -n index.ndx -o ${eqno} -maxwarn 1

	gmx mdrun -v -nsteps 10000 -deffnm ${eqno}
fi

conf_init=${eqno}

done


# EQ with pressure coupling

for dt in 002 004 005 007 010 015 020
do

eqno="eq_dt_${dt}"
mdpfname="martini_3_awh-${eqno}_293K.mdp"

# primitive chekcpointing
if ! [ -s ${eqno}.xtc ] 
then
	# create equilibration mdp file
	sed -e "s/awh1-user-data = yes/awh1-user-data = no/g" \
	    -e 's/Pcoupl                   = parrinello-rahman/Pcoupl                   = berendsen/g' \
	    -e "s/dt                       = 0.02/dt                       = 0.${dt}/g" \
	    ${source_mdp}  > ${mdpfname}

	gmx grompp -f ${mdpfname} -c ${conf_init} -t ${conf_init}.cpt -p system_1mol.top -n index.ndx -o ${eqno}

	gmx mdrun -v -nsteps 100000 -deffnm ${eqno}
fi

conf_init=${eqno}

done
