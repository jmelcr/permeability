#!/bin/bash

# make index file - numbers shall work withing the same version of GMX (2019)
gmx make_ndx -f fluid-ish-30pchol_stitched.pdb <<< ' 2 | 3 | 4
 5 | 6
 q '

# adapt names in index file
sed -e 's/EOL_W/Solvent/g' 
    -e 's/DPPC_POPC_CHOL/Membrane/g'

grep EOL index.ndx -A1 | sed -e 's/EOL/Particle/g' >> index.ndx

# topology adapted manually (duplicating CHOL line=]) as was the configuration (in VMD)

gmx grompp -f ../../martini_3_awh-min.mdp -c fluid-ish-30pchol_stitched.pdb -p system_1mol.top -n index.ndx -o min

gmx mdrun -v -nsteps 5000 -deffnm min
