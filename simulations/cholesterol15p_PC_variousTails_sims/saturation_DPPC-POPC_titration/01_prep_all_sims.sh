#! /bin/bash
# simple script, which prepares a set of simulations
# for permeability measurement
#
# both sturation level as well as permeating particle size is varied

#only phospholipids, which are being titrated
n_lipids_tot=179

rd=`pwd`

mdp_fname="martini_3_awh-md_293K.mdp"

# for various (sizes of) permeating particles
for perm in EOLS EOLT EOLR
do
mkdir -p perm_particle_${perm} && cd perm_particle_${perm} || exit 8

# for various level of (un)saturation
for n_unsat in `seq 0 18 $n_lipids_tot` $n_lipids_tot 
do
n_sat=`echo "$n_lipids_tot - $n_unsat" | bc`
d=`echo "$n_unsat * 50 / $n_lipids_tot ; scale=1" | bc`
echo working on d= $d , n_sat = $n_sat , n_unsat = $n_unsat , total lipids = $n_lipids_tot
mkdir -p sat-unsat_${n_sat}-${n_unsat}_d_${d} && cd sat-unsat_${n_sat}-${n_unsat}_d_${d} || exit 9

# Create system topology
sed -e "s/_SYS_NAME_/sat-unsat_${n_sat}-${n_unsat}_d_${d}/g" \
    -e "s/_N_SAT_/${n_sat}/g" \
    -e "s/_N_UNSAT_/${n_unsat}/g" \
    -e "s/_PERM_/${perm}/g" \
    ${rd}/system_template.top > system.top

# First entry in the molecules list seems to require non-zero number, correcting ...
sed -e '/DPPC           0/d' system.top  -i 

# for various initial conditions - phases
for init_conf in fluid gel
do
mkdir -p sim1_awh_${init_conf} && cd sim1_awh_${init_conf} || exit 11
tsp gmx grompp -c ${rd}/../prep_eq/conf_${init_conf}.gro  -f ${rd}/../../${mdp_fname} -n ${rd}/index.ndx -p ../system.top -maxwarn 1
#tsp -d -N 3 gmx mdrun -awh ${rd}/awh_ref_${init_conf}.xvg -nsteps 10000 -v

cd ../
done

cd ../
done

cd ${rd}
done

