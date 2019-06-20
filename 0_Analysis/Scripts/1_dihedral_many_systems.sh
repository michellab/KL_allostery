# $ script.sh *number_of_systems_in_system_list* *number_of_bins_for_histogram*

rep=0
while [ $rep -lt $1 ];
do
        python 1_dihedral.py $rep $2
	echo $rep
	rep=$((rep+1))

done


