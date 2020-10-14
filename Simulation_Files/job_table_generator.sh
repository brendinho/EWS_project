#!/bin/bash

FILENAME="job_table_temp.txt"

rm job_table.txt
touch $FILENAME

declare -a NET_TYPES=('smallworld');
declare -a SIZES=(10000);
declare -a AVERAGE_NODE_DEGREES=(50)
declare -a ENSEMBLE=($(seq 1 5));
declare -a INFEC_PROB=(0.2);
declare -a VACCINE_RISK=($(seq -1 .03125 1));
declare -a REPLENISHMENT=(2.4e-4);
declare -a IMPORTATION=(2.5e-5);
declare -a SOCIAL_NORM=($(seq 0 0.125 2.5));
declare -a INIT_VACC_PROPORTION=(0.25 0.75);
declare -a RANDOM_OPINION_SWITCH=(0.0001);
declare -a NETWORK_PROPORTION=(0.6 0.8 1);

for size in "${SIZES[@]}"; do
for structure in "${NET_TYPES[@]}"; do
for neighbours in "${AVERAGE_NODE_DEGREES[@]}"; do
for init_vacc in "${INIT_VACC_PROPORTION[@]}"; do

	for switch_prob in "${RANDOM_OPINION_SWITCH[@]}"; do
	for infec in "${INFEC_PROB[@]}"; do

		for replenish in "${REPLENISHMENT[@]}"; do
		for import in "${IMPORTATION[@]}"; do

			for net_prop in "${NETWORK_PROPORTION[@]}"; do
			for instance in "${ENSEMBLE[@]}"; do
			for norm in "${SOCIAL_NORM[@]}"; do
			for risk in "${VACCINE_RISK[@]}"; do

        		echo "$structure $size $neighbours $instance $infec $risk $replenish $import $norm $init_vacc $switch_prob $net_prop" >> $FILENAME;

			done
			done
			done
			done

        done
        done

    done
	done

done
done
done
done

dos2unix $FILENAME

# g++ clean_the_parameter_series.cvpp -o PARAMETER_CLEANER
# . PARAMETER_CLEANER $FILENAME
# rm PARAMETER_CLEANER
# rm $FILENAME
#
# dos2unix job_table.txt
