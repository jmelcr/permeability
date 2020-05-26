#! /bin/bash
# simple script, which chagnes top files from M3 to M2 parameters
#
#  RUN only ONCE!!

rd=`pwd`
topfname="system_1mol.top"
topfnamem3="system_1mol_m3.top"

for d in PC_*
do

cd $d || exit 66
echo $PWD

mv $topfname $topfnamem3

cp ${rd}/template_m2_topology.top $topfname
grep -e " system " -A 20 $topfnamem3 >> $topfname

sed 's/EOLR/EOL /g' $topfname -i

cd $rd
done

