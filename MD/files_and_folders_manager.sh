#!/bin/bash


################################################################################
# This script contains functions that are called from the script lammps_running.sh
# These functions are used for:
#    - Build folders
#    - Build initial conditions and interaction potentials from Silvano's scripts
#    - Write input files of LAMMPS simulations
#    - Write the .slrm files needed to run the simuations with SLURM
################################################################################

# SET OF FUNCTIONS 1: NEEDED TO BUILD FOLDERS

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
# The folders built by this function are dictated by the parametes that set the physics of the model
function build_first_set_of_folders() {

  # Args
  HOME=$1
  mapping=$2
  symmetry=$3
  ipc_model=$4

  model_name=$5

  A=$6
  n=$7


	cd $HOME

  # res folder
	next_fold="res"
	build_dir $next_fold
	cd $next_fold

  # system folder
  if [[ $ipc_model == 1 ]]; then
    next_fold="${symmetry}_${mapping}"
  else
    next_fold="${symmetry}_${mapping}_IPC-OFF"
  fi

	build_dir $next_fold
	cd $next_fold

  # model folder
	next_fold=${model_name}
	build_dir $next_fold
	cd $next_fold

  # LJ-params folder
	next_fold="A${A}_n${n}"
	build_dir $next_fold
	cd $next_fold

  path=$(pwd)
	echo "$path"                                                                  # DO NOT REMOVE!!!! THIS IS THE RETURN LINE OF THE FUNCTION!

}

# This function computes the total number of particles given Nside and the lattice model
# This number is needed by the next function
function compute_Ntot() {

  lattice=$1
  Nside=$2

  if [[ $lattice == SC ]]; then
    Ntot=$((${Nside}*${Nside}*${Nside}))
  elif [[ $lattice == FCC ]]; then
    Ntot=$((4*${Nside}*${Nside}*${Nside}))
  else
    echo "lattice (line 58) can be: SC (Simple Subic) or FCC (Face Centered Cubic)"
    echo "Current selection of lattice is ${lattice}"
    echo "This is not admitted. EXIT FAILURE"
    exit 1
  fi

  echo $Ntot

}


# This function creates a series of folders and the return the ABSOLUTE path to it.
# The folders built by this function are dictated by the setting of the simulations
function build_second_set_of_folders() {

  HOME=$1
  path0=$2

  type=$3
  k=$4
  dump=$5
  T=$6
  Ntot=$7
  rho=$8
  m=$9

  cd $HOME
  cd $path0

  # Thermostat folder
  next_fold="${type}"
  build_dir $next_fold
	cd $next_fold

  # Bond model folder
  if [[ $k == NONE ]]; then
    next_fold="rigid"
  else
    next_fold="k_ang${k}"
  fi

  build_dir $next_fold
	cd $next_fold

  # Dumping factor folder
	if [[ $type == NH ]]; then
  	next_fold="Td${dump}"
	elif [[ $type == LANG ]]; then
  	next_fold="gamma${dump}"
	else
		echo "THERMOSTAT IS NOT KNOWN!"
		echo "The thermostat is ${type}"
		echo "Accepted options are NH (for Nosé-Hoover) and LANG (for Langevin)"
		echo "Interrupting the scritp"
		exit 1
	fi
	build_dir $next_fold
	cd $next_fold

  # Ntot folder
  next_fold="N${Ntot}"
  build_dir $next_fold
	cd $next_fold

  # State point folder
  next_fold="T${T}_rho${rho}"
  build_dir $next_fold
	cd $next_fold

  # run folder
  next_fold="run_num${m}"
  build_dir $next_fold

  path=$(pwd)
	cd $HOME
	echo "$path"                                                                  # DO NOT REMOVE!!!! THIS IS THE RETURN LINE OF THE FUNCTION!

}


################################################################################
################################################################################
################################################################################

# SET OF FUNCTIONS 2: NEEDED TO BUILD POTENTIAL TABLES AND INITIAL CONDITION


# This function creates the tables for the interpartile potentials
# It basically replaces Silivano's script "from_contact_values_to_potentials.sh"
# It returns nothing
function build_potential() {

  cutoff=0.0001

  cd OSPC_LAMMPS/1-make_potentials/

	model_name=$1
	delta=$2
	ecc1=$3
	rad1=$4
	vEE=$5
	vEP1=$6
	vP1P1=$7

	A=$8
	n=$9

	ecc2=${10}
	rad2=${11}
	vEP2=${12}
	vP1P2=${13}
	vP2P2=${14}

  symmetry=${15}
  mapping=${16}

	pushd sources
	  # generate inputfile
	  echo $model_name > inputfile
	  echo $delta >> inputfile
	  echo $cutoff >> inputfile
	  echo $ecc1 >> inputfile
	  echo $rad1 >> inputfile
	  echo $vEE >> inputfile
	  echo $vEP1 >> inputfile
	  echo $vP1P1 >> inputfile

	  echo $A >> inputfile
	  echo $n >> inputfile


	  echo $ecc2 >> inputfile
	  echo $rad2 >> inputfile
	  echo $vEP2 >> inputfile
	  echo $vP1P2 >> inputfile
	  echo $vP2P2 >> inputfile


	  ./build.sh

	  target="../target_${model_name}_${symmetry}_${mapping}_contact"
	  [ -d $target ] && rm -rf $target
	  mkdir -p $target

	  [ $ipc_model -eq 1 ] && is_ipc="-p"
	  ./bld/lammps_pot_generator \
	    -c $is_ipc -s $symmetry -m $mapping \
	    -i inputfile -o $target
	  rm inputfile
	popd

}


# This function creates the initial condition file and return its name to the main script
function build_initial_condition() {

	symmetry=$1
	Nside=$2
	rho=$3
	ecc1=$4
	ecc2=$5
	Ntot=$6

	cd OSPC_LAMMPS/2-startingstate_creators

	if [[ $symmetry == janus ]]; then
		python3 face_centered_cubic_janus.py $Nside $rho $ecc1 > mute
	  init_file="startingstate_FCC_1p_${Ntot}_${rho}_${ecc1}.txt"
	else
		if [[ $lattice = SC ]]; then
			python3 simple_cubic_two_patch.py $Nside $rho $ecc1 $ecc2 > mute
			init_file="startingstate_SC_2p_${Ntot}_${rho}_${ecc1}_${ecc2}.txt"
		elif [[ $lattice == FCC ]]; then
			python3 face_centered_cubic_two_patch.py $Nside $rho $ecc1 $ecc2 > mute
			init_file="startingstate_FCC_2p_${Ntot}_${rho}_${ecc1}_${ecc2}.txt"
		else
	    echo "lattice (line 58) can be: SC (Simple Subic) or FCC (Face Centered Cubic)"
	    echo "Current selection of lattice is ${lattice}"
	    echo "This is not admitted. EXIT FAILURE"
	    exit 1
	  fi

	fi

	echo $init_file

}


################################################################################
################################################################################
################################################################################

# SET OF FUNCTIONS 3: NEEDED TO BUILD INPUT FILES

# This function calls the input files generators, depending on the model (NH/LANG & rigid/harmonic)
# It returns nothing
function write_input_file() {

	m=$1
	full_path=$2
	init_cond_file=$3
	dump=$4
	T=$5
	k=$6
	ecc1=$7

	if [ $type == NH ] && [ $k == NONE ]; then
		write_NH_rigid_input "$m" "$full_path" "$init_cond_file" "$dump" "$T"
	elif [ $type == NH ] && [ $k != "NONE" ]; then
		write_NH_harmonic_input "$m" "$full_path" "$init_cond_file" "$dump" "$T" "$k" "$ecc1"
	elif [ $type == LANG ] && [ $k == NONE ]; then
		write_LANG_rigid_input "$m" "$full_path" "$init_cond_file" "$dump" "$T"
	elif [ $type == LANG ] && [ $k != "NONE" ]; then
		write_LANG_harmonic_input "$m" "$full_path" "$init_cond_file" "$dump" "$T" "$k" "$ecc1"
	fi



}

# This function build the input file for a NH simulation with rigid bonds.
# Once the input is written, it is copied in the folder where the output of the simulation
# will be written. After that, the function terminates without returning anything
function write_NH_rigid_input() {

	m=$1
	full_path=$2
	init_cond_file=$3
	Td=$4
	T=$5

	input_file="input_num${m}.txt"
	output_file="output_num${m}.txt"
	equilibration_file="equil_num${m}.lammpstrj"
	trajectory_file="traj_num${m}.lammpstrj"

	cp templates/NH_rigid_2patches_template $input_file

	sed -i "15 c\read_data ${full_path}/${init_cond_file}" ${input_file}

	sed -i "21 c\pair_coeff 1 1 ${full_path}/lammpspot/BB.table BB" ${input_file}
	sed -i "22 c\pair_coeff 1 2 ${full_path}/lammpspot/Bs1.table Bs1" ${input_file}
	sed -i "23 c\pair_coeff 1 3 ${full_path}/lammpspot/Bs2.table Bs2" ${input_file}
	sed -i "24 c\pair_coeff 2 2 ${full_path}/lammpspot/s1s1.table s1s1" ${input_file}
	sed -i "25 c\pair_coeff 2 3 ${full_path}/lammpspot/s1s2.table s1s2" ${input_file}
	sed -i "26 c\pair_coeff 3 3 ${full_path}/lammpspot/s2s2.table s2s2" ${input_file}

	sed -i "35 c\log ${full_path}/run_num${m}/${output_file}" ${input_file}

	sed -i "39 c\dump           dumpy   all atom 10000 ${full_path}/run_num${m}/${equilibration_file}" ${input_file}
	sed -i "52 c\dump           dumpy   all atom 10000 ${full_path}/run_num${m}/${trajectory_file}" ${input_file}

	seed=$(shuf -i 1-800000000 -n 1)
	sed -i "41 c\velocity       all create 1.0 ${seed} dist gaussian rot yes" ${input_file}

	Td_command=$(printf '$(%.1f*dt)' "$Td")

	sed -i "42 c\fix    1       all rigid/nvt/small molecule temp 1.0 1.0 ${Td_command}" ${input_file}
	sed -i "46 c\fix    1       all rigid/nvt/small molecule temp 1.0 ${T} ${Td_command}" ${input_file}
	sed -i "54 c\fix    1       all rigid/nvt/small molecule temp ${T} ${T} ${Td_command}" ${input_file}

	cp $input_file ${full_path}/run_num${m}
	rm $input_file

}

# This function build the input file for a LANG simulation with rigid bonds.
# Once the input is written, it is copied in the folder where the output of the simulation
# will be written. After that, the function terminates without returning anything
function write_LANG_rigid_input() {

	m=$1
	full_path=$2
	init_cond_file=$3
	gamma=$4
	T=$5

	input_file="input_num${m}.txt"
	output_file="output_num${m}.txt"
	equilibration_file="equil_num${m}.lammpstrj"
	trajectory_file="traj_num${m}.lammpstrj"

	cp templates/LANG_rigid_2patches_template $input_file

	sed -i "15 c\read_data ${full_path}/${init_cond_file}" ${input_file}

	sed -i "21 c\pair_coeff 1 1 ${full_path}/lammpspot/BB.table BB" ${input_file}
	sed -i "22 c\pair_coeff 1 2 ${full_path}/lammpspot/Bs1.table Bs1" ${input_file}
	sed -i "23 c\pair_coeff 1 3 ${full_path}/lammpspot/Bs2.table Bs2" ${input_file}
	sed -i "24 c\pair_coeff 2 2 ${full_path}/lammpspot/s1s1.table s1s1" ${input_file}
	sed -i "25 c\pair_coeff 2 3 ${full_path}/lammpspot/s1s2.table s1s2" ${input_file}
	sed -i "26 c\pair_coeff 3 3 ${full_path}/lammpspot/s2s2.table s2s2" ${input_file}

	sed -i "35 c\log ${full_path}/run_num${m}/${output_file}" ${input_file}

	sed -i "39 c\dump           dumpy   all atom 10000 ${full_path}/run_num${m}/${equilibration_file}" ${input_file}
	sed -i "52 c\dump           dumpy   all atom 10000 ${full_path}/run_num${m}/${trajectory_file}" ${input_file}

	seed1=$(shuf -i 1-800000000 -n 1)
	seed2=$(shuf -i 1-800000000 -n 1)
	seed3=$(shuf -i 1-800000000 -n 1)

	sed -i "42 c\fix    1       all rigid/small molecule langevin 1.0 1.0 ${gamma} ${seed1}" ${input_file}
	sed -i "46 c\fix    1       all rigid/small molecule langevin 1.0 ${T} ${gamma} ${seed2}" ${input_file}
	sed -i "54 c\fix    1       all rigid/small molecule langevin ${T} ${T} ${gamma} ${seed3}" ${input_file}

	cp $input_file ${full_path}/run_num${m}
	rm $input_file

}



# This function build the input file for a NH simulation with harmonic bonds.
# Once the input is written, it is copied in the folder where the output of the simulation
# will be written. After that, the function terminates without returning anything
function write_NH_harmonic_input() {

	m=$1
	full_path=$2
	init_cond_file=$3
	Td=$4
	T=$5
	k=$6
	ecc1=$7

	# kappa=$(($k*1000))

	input_file="input_num${m}.txt"
	output_file="output_num${m}.txt"
	equilibration_file="equil_num${m}.lammpstrj"
	trajectory_file="traj_num${m}.lammpstrj"

	cp templates/NH_hb_2patches_template $input_file

	sed -i "15 c\read_data ${full_path}/${init_cond_file}" ${input_file}

	sed -i "23 c\pair_coeff 1 1 ${full_path}/lammpspot/BB.table BB" ${input_file}
	sed -i "24 c\pair_coeff 1 2 ${full_path}/lammpspot/Bs1.table Bs1" ${input_file}
	sed -i "25 c\pair_coeff 1 3 ${full_path}/lammpspot/Bs2.table Bs2" ${input_file}
	sed -i "26 c\pair_coeff 2 2 ${full_path}/lammpspot/s1s1.table s1s1" ${input_file}
	sed -i "27 c\pair_coeff 2 3 ${full_path}/lammpspot/s1s2.table s1s2" ${input_file}
	sed -i "28 c\pair_coeff 3 3 ${full_path}/lammpspot/s2s2.table s2s2" ${input_file}
	# sed -i "29 c\bond_coeff 1 ${kappa} ${ecc1}" ${input_file}
	# sed -i "30 c\bond_coeff 2 ${kappa} ${ecc1}" ${input_file}
	sed -i "29 c\bond_coeff 1 10000 ${ecc1}" ${input_file}
	sed -i "30 c\bond_coeff 2 10000 ${ecc1}" ${input_file}
	sed -i "31 c\angle_coeff 1 ${k} 180.0" ${input_file}

	sed -i "40 c\log ${full_path}/run_num${m}/${output_file}" ${input_file}

	sed -i "44 c\dump           dumpy   all atom 20000 ${full_path}/run_num${m}/${equilibration_file}" ${input_file}
	sed -i "57 c\dump           dumpy   all atom 20000 ${full_path}/run_num${m}/${trajectory_file}" ${input_file}

	seed=$(shuf -i 1-800000000 -n 1)
	sed -i "46 c\velocity       all create 1.0 ${seed} dist gaussian rot yes" ${input_file}


	Td_command=$(printf '$(%.1f*dt)' "$Td")

	sed -i "47 c\fix    1       all nvt temp 1.0 1.0 ${Td_command}" ${input_file}
	sed -i "51 c\fix    1       all nvt temp 1.0 ${T} ${Td_command}" ${input_file}
	sed -i "59 c\fix    1       all nvt temp ${T} ${T} ${Td_command}" ${input_file}

	cp $input_file ${full_path}/run_num${m}
	rm $input_file

}


# This function build the input file for a LANG simulation with harmonic bonds.
# Once the input is written, it is copied in the folder where the output of the simulation
# will be written. After that, the function terminates without returning anything
function write_LANG_harmonic_input() {

	m=$1
	full_path=$2
	init_cond_file=$3
	gamma=$4
	T=$5
	k=$6
	ecc1=$7

	input_file="input_num${m}.txt"
	output_file="output_num${m}.txt"
	equilibration_file="equil_num${m}.lammpstrj"
	trajectory_file="traj_num${m}.lammpstrj"

	cp templates/LANG_hb_2patches_template $input_file

	sed -i "15 c\read_data ${full_path}/${init_cond_file}" ${input_file}

	sed -i "23 c\pair_coeff 1 1 ${full_path}/lammpspot/BB.table BB" ${input_file}
	sed -i "24 c\pair_coeff 1 2 ${full_path}/lammpspot/Bs1.table Bs1" ${input_file}
	sed -i "25 c\pair_coeff 1 3 ${full_path}/lammpspot/Bs2.table Bs2" ${input_file}
	sed -i "26 c\pair_coeff 2 2 ${full_path}/lammpspot/s1s1.table s1s1" ${input_file}
	sed -i "27 c\pair_coeff 2 3 ${full_path}/lammpspot/s1s2.table s1s2" ${input_file}
	sed -i "28 c\pair_coeff 3 3 ${full_path}/lammpspot/s2s2.table s2s2" ${input_file}
	# sed -i "29 c\bond_coeff 1 ${kappa} ${ecc1}" ${input_file}
	# sed -i "30 c\bond_coeff 2 ${kappa} ${ecc1}" ${input_file}
	sed -i "29 c\bond_coeff 1 10000 ${ecc1}" ${input_file}
	sed -i "30 c\bond_coeff 2 10000 ${ecc1}" ${input_file}
	sed -i "31 c\angle_coeff 1 ${k} 180.0" ${input_file}

	sed -i "40 c\log ${full_path}/run_num${m}/${output_file}" ${input_file}

	sed -i "44 c\dump           dumpy   all atom 20000 ${full_path}/run_num${m}/${equilibration_file}" ${input_file}
	sed -i "57 c\dump           dumpy   all atom 20000 ${full_path}/run_num${m}/${trajectory_file}" ${input_file}

	seed1=$(shuf -i 1-800000000 -n 1)
	seed2=$(shuf -i 1-800000000 -n 1)
	seed3=$(shuf -i 1-800000000 -n 1)

	sed -i "47 c\fix    1       all langevin 1.0 1.0 ${gamma} ${seed1}" ${input_file}
	sed -i "51 c\fix    1       all langevin 1.0 ${T} ${gamma} ${seed2}" ${input_file}
	sed -i "59 c\fix    1       all langevin ${T} ${T} ${gamma} ${seed3}" ${input_file}

	cp $input_file ${full_path}/run_num${m}
	rm $input_file

}

################################################################################
################################################################################
################################################################################

# SET OF FUNCTIONS 4: NEEDED TO BUILD SLURM FILES


function init_new_exec() {

	max_tasks=$1

	exec_file="runner.sh"
	cp templates/runner_template.sh ${exec_file}

	sed -i "20 c\max_tasks=${max_tasks}" ${exec_file}

}


function update_exec() {

	m=$1
	full_path=$2
	Ncores=$3

	exec_file="runner.sh"

	printf '\nin_file="%s/run_num%d/input_num%d.txt"\n' "${full_path}" "$m" "$m" >> $exec_file
	printf 'log_file="%s/run_num%d/log_num%d.lammps"\n' "${full_path}" "$m" "$m" >> $exec_file
	# printf 'srun --mem=$mem_per_task --cpus-per-task=1 --ntasks=%d ${LAMMPSpath} -screen none -in ${in_file} -log ${log_file} &\n' "${Ncores}" >> $exec_file   # FOR SLURM
	printf 'mpirun -np=%d --bind-to none ${LAMMPSpath} -screen none -in ${in_file} -log ${log_file} &\n' "${Ncores}" >> $exec_file                               # STD
	printf 'sleep 10\n' >> $exec_file
	printf 'check_and_sleep "$max_tasks"\n' >> $exec_file

}



function close_exec() {


	exec_file="runner.sh"

	printf '\n\nwait\n' >> $exec_file

}


# function update_SLRM_Emanuele() {
#
# 	Nslrm=$1
# 	m=$2
# 	full_path=$3
# 	Ncores=$4
# 	TheStr=$5
#
# 	exec_file="SLRM_runner_${Nslrm}.slrm"
#
# 	printf '\nin_file="%s/run_num%d/input_num%d.txt"\n' "${full_path}" "$m" "$m" >> $exec_file
# 	printf 'log_file="%s/run_num%d/log_num%d.lammps"\n' "${full_path}" "$m" "$m" >> $exec_file
#
# 	# printf 'TheStr=%s' "${TheStr}"  >> $exec_file
# 	# printf 'srun -n 4 --overlap --cpu_bind=mask_cpu:%s lmp_intel_cpu_intelmpi -screen none -in ${in_file} -log ${log_file} &'  "${TheStr}" >> $exec_file
# 	printf 'srun -n 4 --overlap --cpu_bind=maks lmp_intel_cpu_intelmpi -screen none -in ${in_file} -log ${log_file} &\n'  "${TheStr}" >> $exec_file
#
# 	# printf 'mpirun -n %d ${LAMMPSpath} -screen none -in ${in_file} -log ${log_file} &\n' "${Ncores}" >> $exec_file
# 	printf 'sleep 10\n' >> $exec_file
# 	printf 'check_and_sleep_Emanuele "$max_tasks"\n' >> $exec_file
#
# }


















#
