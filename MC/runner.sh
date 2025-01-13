#!/bin/bash

tasks_to_be_done=28
max_tasks=8
current_task=0
running_tasks=0

ENS=NVT
DYN=VMMC

declare -a _params=(
	 "1.0000 250 1"  "1.0000 250 2"  "1.0000 250 3"  "1.0000 250 4"  "1.0000 500 1"  "1.0000 500 2"  "1.0000 500 3"  "1.0000 500 4" 
	 "1.1000 250 1"  "1.1000 250 2"  "1.1000 250 3"  "1.1000 250 4"  "1.1000 500 1"  "1.1000 500 2"  "1.1000 500 3"  "1.1000 500 4" 
	 "1.2000 250 1"  "1.2000 250 2"  "1.2000 250 3"  "1.2000 250 4"  "1.2000 500 1"  "1.2000 500 2"  "1.2000 500 3"  "1.2000 500 4" 
	 "1.3000 100 1"  "1.3000 100 2"  "1.3000 100 3"  "1.3000 100 4" 
)

full_path=/home/danimuzi/Desktop/MD_IPPs/code_for_git/MC/res/EXPO/NVT_VMMC/EE0.1_PP4.0_D0.2_K13/L10/RotMax0.1_TranMax0.05/Smax25_Tmax1.8


current_task=0
while (($current_task < $tasks_to_be_done)); do
	running_tasks=`ps -C IPC --no-headers | wc -l`
	while (($running_tasks < $max_tasks && ${current_task} < ${tasks_to_be_done})); do

		set -- ${_params[$current_task]}
		T=${1}
		N=${2}
		LastFold="T${T}_N${N}"

		n=${3}

		f1="${LastFold}/settings_num${n}.txt"

		f="${full_path}/${f1}"
		# echo "$f"

		./IPC $f &

		running_tasks=`ps -C IPC --no-headers | wc -l`
		((current_task++))

	done
done
wait
