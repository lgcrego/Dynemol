#!/bin/bash 

###################################################################
#Script Name	: panel_elhl_wp
#Description	: gathers el_wp_energy.dat and hl_wp_energy.dat in the local directory and 
#                 plot a panel of wavepacket energies for each sample simulation                                                                              
#Args           :                                                                                           
#Author       	: Luis G C Rego
#date           : 12/Feb/2024
###################################################################

#  >>>  edit paramters  <<<
#------------------------------------------------------------------
#number of files , rows and columns of the panel
nf=50    nrows=10    ncolumns=5

#gather the files in the local directory
for i in `seq -w 01 ${nf}`; do 
      cp ../inpt-"${i}"/dyn.trunk/el_wp_energy.dat el_data-"${i}" 
      cp ../inpt-"${i}"/dyn.trunk/hl_wp_energy.dat hl_data-"${i}"
done
#------------------------------------------------------------------

# Initialize an empty array
declare -a min_vals
declare -a max_vals

#find the min and max values for plotting the graphs
for i in $(seq -w 01 ${nf}); do
    # Convert padded number to integer
    index=$(echo $i | sed 's/^0*//')
    # Capture the output of the awk command and assign it to the array element
    min_vals[$index]=$(awk 'NR == 1 || $2 < min {min = $2} END {print min}' hl_data-$i)
    max_vals[$index]=$(awk 'NR == 1 || $2 > max {max = $2} END {print max}' el_data-${i}) 
done

y_min=${min_vals[1]}
y_max=${max_vals[1]}
for ((i = 2; i <= ${nf}; i++)); do
    #finding y_min 
    if (( $(echo "${min_vals[$i]} < ${y_min}" | bc -l) )); then
        y_min=${min_vals[$i]}
    fi
    #finding y_max 
    if (( $(echo "${max_vals[$i]} > ${y_max}" | bc -l) )); then
        y_max=${max_vals[$i]}
    fi
done

echo "parameters: " "${y_min}" "${y_max}" "${nf}" "${nrows}" "${ncolumns}"

gnuplot -c gnupanel.gps "${y_min}" "${y_max}" "${nf}" "${nrows}" "${ncolumns}"

