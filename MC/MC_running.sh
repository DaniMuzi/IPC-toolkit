#!/bin/bash

export LC_ALL=en_US.UTF-8

################################################################################
################################################################################

# This script is used to construct folders and input files for MC simulations
# of the IPP model. Its usage is detailed in the associated README. However,
# here we briefly summarize the parameters that a user might want to change.
# They are comprised between line 84 and line 116 and are:

# - max_translation = maximum displacement of the seed of the VMMC cluster
# - max_rotation = maximum rotation of the seed of the VMMC
# - L = box size
# - max_cluster_size = maximum VMMC cluster size
# - max_cluster_translation = maximum translation length of a VMMC cluster
# - tot_steps = number of MC step of a simulation
# - save_steps = configurations are saved in output by the simulation every print_steps MC steps
# - print_steps = energy and acceptance rate are saved in output by the simulation every print_steps MC steps
# - restart = whether or not the simulation starts from scratch (set to 0) or if it restarts from a previously saved configuration (set to 1)
# - max_tasks = maximum number of parallel simulations to be run
# - runs_per_conf = number of parallel simulations per state point

# - _params =  the array params contains pairs. Each pair is of the form (T, N) where T is the temperature of the simulation and N is the number of particles.
#              Each pair represents a state point.

################################################################################
################################################################################

## INIT: basic variables of the script, functions, input variables

source ./funcs.sh

# Interaction: 0 = overlap of spheres, 1 = exponential
model_file=$1
interaction=$2

ensemble=0
dynamics=1

################################################################################

## Input check

# Ensemble is allowed?
if [[ $ensemble == 0 ]]; then
	ENS="NVT"
else
	echo "UNKNOWN ENSEMBLE"
	exit 1
fi

# Dynamics is allowed?
if [[ $dynamics == 1 ]]; then
	DYN="VMMC"
else
	echo "UNKNOWN DYNAMICS"
	exit 1
fi

# Interaction is allowed?
if [[ $interaction == 0 ]]; then
	INT="SPHERES"
elif [[ $interaction == 1 ]]; then
	INT="EXPO"
else
	echo "UNKNOWN INTERACTION"
	exit 1
fi

echo "Using $ENS with $DYN and $INT"

################################################################################
################################################################################
################################################################################

## Settings for the physics of the model

n_patches=2

read -r a sigma_c sigma_p delta_c delta_p e_EE e_EP e_PP k_EE k_EP k_PP rcut <<< "$(read_model $model_file)"

################################################################################

## Settings of a single MC simulation

max_translation=0.05
max_rotation=0.1
L=10

# VMMC settings

max_cluster_size=25                               													# Try to keep this < Nmax/10 (at least)
max_cluster_translation=1.8                    														  #$Try to keep this <~ L/5

# "Unphysical" details of the MC simulation

tot_steps=82000 #000                                                               # Total number of steps
save_steps=500 #00                                                                # Configuration print frequency
print_steps=100 #0                                                                # Observables print frequency
restart=0                                                                       # Wether to start the simulation from t=0 or to restarty from a configuration saved. 0 -> NO restarting, 1 YES

################################################################################

max_tasks=8
runs_per_conf=5                                                                 # How many MC simulations for each given state point

################################################################################

## State points: each is pair of the for (T, N)

declare -a _params=("1.0000 250" "1.0000 500" "1.1000 250" "1.1000 500" "1.2000 250" "1.2000 500" "1.3000 100")
num_of_confs=${#_params[@]}                                                     # Number of state points

################################################################################
#																			 #
#           				CREATION OF FOLDERS AND INPUT FILES           			       #
#																			 #
################################################################################


# Paths (part 1)

HOME=$(pwd)
model=$(basename ${model_file} .txt)

# Build folders build paths (part 2)

cd $HOME
full_path=$(build_directories "$HOME" "$ENS" "$DYN" "$INT" "$model" "$L" "$max_rotation" "$max_translation" "$max_cluster_size" "$max_cluster_translation")
cd $full_path

## NOW WE BUILD THE LAST SET OF FOLDERS AND WE BUILD THE INPUT FILES

i=0
while (($i < $num_of_confs)); do

  # Extract state point from the array _params
	set -- ${_params[$i]}
	T=${1}
	N=${2}
	LastFold="T${T}_N${N}"

  # Build directory of the state point and enter it
	build_dir $LastFold
	cd $LastFold

	j=0
	while (($j < $runs_per_conf)); do
		n=$((j+1))

    # Build directory for files strictly related to run number n: contains configurations generated during the run
		fold="confs_num${n}"
		build_dir $fold

		if [[ $restart == 0 ]]; then                                                # Generate seed of simulation and give name to initial condition file
			seed=$(shuf -i 1-2147483645 -n 1)
			init_cond_file="init_cond_num$n.txt"
		else
			setf="settings_num$n.txt"
			l=$(sed '36!d' $setf)
			seed=${l##* }
			init_cond_file="confs_num${n}/last.rrr"
		fi

    ## Now we write the settings file

		write_mandatory_settings "$n" "$tot_steps" "$print_steps" "$save_steps" "$dynamics" "$ensemble" "$interaction" "$T" "$max_translation" "$max_rotation" "$L" "$N" "$full_path" "$LastFold" "$init_cond_file"
		write_VMMC_settings "$n" "$max_cluster_translation" "$max_cluster_size"
		write_non_mandatory_settings "$n" "$restart" "$seed" "$full_path" "$LastFold"
		write_particle_system_settings "$n" "$n_patches" "$delta_c" "$sigma_c" "$delta_p" "$sigma_p" "$a" "$e_EE" "$e_EP" "$e_PP" "$k_EE" "$k_EP" "$k_PP" "$rcut"

		if [[ $restart == 0 ]]; then                                                # Generates the initial condition file if it is needed
			cd $HOME
			f="settings_num${n}.txt"
			fname="${full_path}/${LastFold}/${f}"
			./init_cond_creator.sh $fname &
			wait
			cd $full_path/$LastFold
		fi

		((j++))
	done

	cd $full_path
	((i++))
done

################################################################################
#																			 #
#           							  CREATION OF ACTUAL RUNNER  		           			     #
#																			 #
################################################################################


cd $HOME

stop=8
ntasks=0

i=0
n=1
line="declare -a _params=(\n\t"                                                 # The variable line contains the array _params to be put in the runner.sh file

while (($i < $num_of_confs)); do

  # Extract state point from the array _params
	set -- ${_params[$i]}
	T=${1}
	N0=${2}
	LastFold="T${T}_N${N0}"

  m=1
  while (($m < $runs_per_conf)); do

		line="${line} \"${T} ${N0} ${m}\" "

  	((ntasks++))

  	r=$((ntasks % stop))
  	if [[ $r = 0 ]]; then
  		line="${line}\n\t"
  	fi

  	((m++))
  done
  ((i++))
done





line="${line}\n)\n\n"

write_runner_header "$ntasks" "$max_tasks" "$HOME" "$ENS" "$DYN"

f="runner.sh"
printf "${line}" >> $f
printf "full_path=${full_path}\n" >> $f

write_runner
echo "Creating runner.sh"

((n++))
ntasks=0
line="declare -a _params=(\n\t"



################################################################################
#																			 #
#           					LAUNCHING AND MONITORING SLURM JOBS	 	 					    		 #
#																			 #
################################################################################

cd $HOME
chmod u+x runner.sh                                                          # Need to make the .slrm files executables
# ./runner.sh
#
