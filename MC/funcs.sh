#!/bin/bash

# SET OF FUNCTIONS 1: NEEDED TO BUILD FOLDERS

# This function takes in input a string representing a file name.
# The file MUST have the structure explained in the README
# Given such file, the function returns all the parameters that are relevant to the chosen IPP model

function read_model() {

	f=$1

	l=$(sed '1!d' $f)
	a=${l##* }

	l=$(sed '2!d' $f)
	sigma_c=${l##* }

	l=$(sed '3!d' $f)
	sigma_p=${l##* }

	l=$(sed '4!d' $f)
	delta_c=${l##* }

	l=$(sed '5!d' $f)
	delta_p=${l##* }

	l=$(sed '6!d' $f)
	e_EE=${l##* }

	l=$(sed '7!d' $f)
	e_EP=${l##* }

	l=$(sed '8!d' $f)
	e_PP=${l##* }

	l=$(sed '9!d' $f)
	k_EE=${l##* }

	l=$(sed '10!d' $f)
	k_EP=${l##* }

	l=$(sed '11!d' $f)
	k_PP=${l##* }

	l=$(sed '12!d' $f)
	rcut=${l##* }

	echo "$a $sigma_c $sigma_p $delta_c $delta_p $e_EE $e_EP $e_PP $k_EE $k_EP $k_PP $rcut"

}

## Functions for directories

# This function takes in input a string.
# If a folder with name given by the input string (locally) exists, it does nothing
# If the folder does not (locally) exists, the function (locally) creates it

function build_dir() {
	if [[ ! -d $1 ]]; then
		mkdir $1 &
		wait
	fi
}

# This function creates a series of folders and the return the ABSOLUTE path to it.
# The folders built by this function are dictated by the parametes that set the physics of the model and by the setting of the MC simulation

function build_directories() {

	HOME=$1
	ENS=$2
	DYN=$3
	INT=$4
	model=$5
	L=$6
	max_rotation=$7
	max_translation=$8
	max_cluster_size=$9
	max_cluster_translation=${10}

	cd $HOME
	full_path=$HOME

	fold="res"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"


	fold="${INT}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"


	fold="${ENS}_${DYN}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"


	fold="${model}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"


	fold="L${L}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"

	fold="RotMax${max_rotation}_TranMax${max_translation}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"

	# Dynamics folder
	fold="Smax${max_cluster_size}_Tmax${max_cluster_translation}"
	build_dir $fold
	cd $fold
	full_path="${full_path}/${fold}"

	echo "$full_path"                                                             # DO NOT REMOVE!!!! THIS IS THE RETURN LINE OF THE FUNCTION!

}


# SET OF FUNCTIONS 2: NEEDED TO BUILD INPUT FILES

# These functions takes in input the various parameters that have been set for the simulations
# After they have all been called, a proper input file will be generated and placed in the proper folder
# Note that all these functions will be called, despite not all parameters are needed.
# The reason behind this is to have settings file that always have the same structure
# In particular, the seed is always at line 59, which is necessary for the restart option
# However, (say) we use NVT: the function for GC and NPT will write ????

# MANDATORY: these parameters are needed

function write_mandatory_settings() {

	n=$1
	tot_steps=$2
	print_steps=$3
	save_steps=$4
	dynamics=$5
	ensemble=$6
	interaction=$7
	T=$8
	max_translation=$9
	max_rotation=${10}
	L=${11}
	N=${12}
	full_path=${13}
	LastFold=${14}
	init_cond_file=${15}

	f="settings_num${n}.txt"

	printf "###############################################\n" > $f               # The file is created: we have > $f
	printf "#                   MANDATORY                 #\n" >> $f              # The file exists and lines are appended: we have >> $f
	printf "###############################################\n" >> $f
	printf "\n" >> $f
	printf "Steps = ${tot_steps}\n" >> $f
	printf "Print_every = ${print_steps}\n" >> $f
	printf "Save_every = ${save_steps}\n" >> $f
	printf "\n" >> $f
	printf "# Dynamics: 0 = RTMC, 1 = VMMC\n# Ensamble: 0 = NVT, 1 = GC, 2 = NPT\n" >> $f
	printf "Dynamics = ${dynamics}\n" >> $f
	printf "Ensemble = ${ensemble}\n" >> $f
	printf "Interaction = ${interaction}\n" >> $f
	printf "Temperature = ${T}\n" >> $f
	printf "\n" >> $f
	printf "Disp_max = ${max_translation}\n" >> $f
	printf "Theta_max = ${max_rotation}\n" >> $f
	printf "\n" >> $f
	printf "L = ${L}\n" >> $f
	printf "N = ${N}\n" >> $f
	printf "\n" >> $f
	printf "Initial_conditions_file = ${full_path}/${LastFold}/${init_cond_file}\n" >> $f
	printf "\n" >> $f

}


function write_VMMC_settings() {

	n=$1
	max_cluster_translation=$2
	max_cluster_size=$3

	f="settings_num${n}.txt"

	printf "###############################################\n" >> $f
	printf "#               MANDATORY FOR VMMC            #\n" >> $f
	printf "###############################################\n" >> $f
	printf "\n" >> $f
	printf "vmmc_max_move = ${max_cluster_translation}\n" >> $f
	printf "vmmc_max_cluster = ${max_cluster_size}\n" >> $f
	printf "\n" >> $f

}

# Mandatory only if Lx_move is True

function write_non_mandatory_settings() {

	n=$1
	restart=$2
	seed=$3
	full_path=$4
	LastFold=$5

	f="settings_num${n}.txt"

	printf "###############################################\n" >> $f
	printf "#                 NON MANDATORY               #\n" >> $f
	printf "###############################################\n" >> $f
	printf "\n" >> $f
	printf "Restart_step_counter = ${restart}\n" >> $f
	printf "Seed = ${seed}\n" >> $f
	printf "\n" >> $f
	printf "Log_file = ${full_path}/${LastFold}/log_num${n}.txt\n" >> $f
	printf "Energy_file = ${full_path}/${LastFold}/energy_num${n}.txt\n" >> $f
	printf "Acceptance_file = ${full_path}/${LastFold}/acceptance_num${n}.txt\n" >> $f
	printf "Configuration_folder = ${full_path}/${LastFold}/confs_num${n}/\n" >> $f
	printf "Configuration_last = ${full_path}/${LastFold}/confs_num${n}/last.rrr\n" >> $f
	printf "\n" >> $f

}

# MANDATORY: will be called, these parameters are needed.
# The parameters for the exponential coarse grain are absent when we use SPHERES
# It's ok

function write_particle_system_settings() {

	n=$1
	n_patches=$2
	delta_c=$3
	sigma_c=$4
	delta_p=$5
	sigma_p=$6
	a=$7
	e_EE=$8
	e_EP=$9
	e_PP=${10}
	k_EE=${11}
	k_EP=${12}
	k_PP=${13}
	rcut=${14}

	f="settings_num${n}.txt"

	printf "###############################################\n" >> $f
	printf "#                PARTICLE SYSTEM              #\n" >> $f
	printf "###############################################\n" >> $f
	printf "\n" >> $f
	printf "n_patches = ${n_patches}" >> $f
	printf "\n" >> $f
	printf "delta_c = ${delta_c}\n" >> $f
	printf "sigma_c = ${sigma_c}\n" >> $f
	printf "delta_p = ${delta_p}\n" >> $f
	printf "sigma_p = ${sigma_p}\n" >> $f
	printf "a = ${a}\n" >> $f
	printf "\n" >> $f
	printf "e_EE = ${e_EE}\n" >> $f
	printf "e_EP = ${e_EP}\n" >> $f
	printf "e_PP = ${e_PP}\n" >> $f
	printf "\n" >> $f
	printf "k_EE = ${k_EE}\n" >> $f
	printf "k_EP = ${k_EP}\n" >> $f
	printf "k_PP = ${k_PP}\n" >> $f
	printf "rcut = ${rcut}\n" >> $f

}


# SET OF FUNCTIONS 3: NEEDED TO BUILD THE .sh FILE

function write_runner_header() {

	ntasks=$1
	max_tasks=$2
	HOME=$3
	ENS=$4
	DYN=$5

	f="runner.sh"

	printf "#!/bin/bash\n" > $f
	printf "\n" >> $f
	printf "tasks_to_be_done=${ntasks}\n" >> $f
	printf "max_tasks=${max_tasks}\n" >> $f
	printf "current_task=0\n" >> $f
	printf "running_tasks=0\n" >> $f
	printf "\n" >> $f
	printf "ENS=${ENS}\n" >> $f
	printf "DYN=${DYN}\n" >> $f
	printf "\n" >> $f

}

# Second part of the .slrm file

function write_runner() {

	f="runner.sh"

	IFS=''
	while read Fin; do
		printf "${Fin}\n" >> $f
	done < runner_bottomfile.sh

}





















#
