#!/bin/bash

conf_init_stitched="conf_CHOL35p_stitched.pdb"

# make index file - numbers shall work withing the same version of GMX (2019)
gmx make_ndx -f ${conf_init_stitched} <<< ' 2 | 3 | 4
 5 | 6
 q 
'

# adapt names in index file
sed -e 's/EOL_W/Solvent/g' \
    -e 's/DPPC_POPC_CHOL/Membrane/g' \
    index.ndx -i

grep EOL index.ndx -A1 | sed -e 's/EOL/Particle/g' >> index.ndx

# topology adapted manually (duplicating CHOL line=]) as was the configuration (in VMD)

gmx grompp -f ../../martini_3_awh-min.mdp -c ${conf_init_stitched} -p system_1mol.top -n index.ndx -o min

gmx mdrun -v -nsteps 50000 -deffnm min
