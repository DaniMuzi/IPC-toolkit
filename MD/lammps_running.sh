#!/bin/bash

export LC_ALL=en_US.UTF-8

################################################################################

HOME=$(pwd)
source ./files_and_folders_manager.sh                                           # A script containing shell functions. It's like "import" in python or <#include> in C/C++

################################################################################

## Settings for the physics of the model

mapping=g                                                                       # g = geometrical, e = exponential
symmetry="symm"                                                                 # symm / asymm
ipc_model=1

delta=0.2                                                                       # delta = TWICE the distance from the colloid surface (radius set to 0.5) where the potential goes to zero
ecc1=0.22                                                                       # NOTE: delta / 2 + 0.5 = ecc1 + rad1 MUST be true
rad1=0.38

# mapping=e                                                                       # g = geometrical, e = exponential
# symmetry="symm"                                                                 # symm / asymm
# ipc_model=1

# parameters for the first patch
# delta=0.077
# ecc1=0.22
# rad1=0.3185


# contact values -> the naming is referring to two patch systems.
vEE=0.1
vEP1=-1.0
vP1P1=4.0

# parameters for regualble softness
A=500
n=15

# parameters for the second patch -- **IGNORED** for tsIPC
# written as it is, the following settings assume that the IPC constraint is ON: we set the repulsions to be the same as for the first patch
ecc2=$ecc1
rad2=$rad1
vEP2=$vEP1
vP1P2=$vP1P1
vP2P2=$vP1P1

# This line might require changes if the IPC constraint is off
model_name="EE${vEE}_PP${vP1P1}_d${delta}_a${ecc1}"

path0=$(build_first_set_of_folders "$HOME" "$mapping" "$symmetry" "$ipc_model" "$model_name" "$A" "$n")
cd $HOME

################################################################################
## Settings for the simulations

lattice='SC'                                                                    # SC / FCC

# One entry of the array _param is a triplet of the form (Temperature, Nside, rho).
declare -a _params=("0.15 10 0.25" "0.15 10 0.5" "0.15 10 0.75")

# One entry of the array _thermostat is a triplet. Each triplet is of the form (A B C), where:
# A = NH / LANG = Nosé–Hoover / Langevin
# B = spring constant for harmonic bonds. NONE stands for rigid bonds
# C = Value for dumping effect = T_d for NH and gamma for LANG.
declare -a _thermostat=("NH 5200 100.0" "NH NONE 100.0" "LANG 5200 1.0" "LANG NONE 1.0")

NstatePoints=${#_params[@]}
Nalgorithms=${#_thermostat[@]}

Nparal=8                                                                        # How many parallel runs for each simulation
Ncores=4                                                                        # How many cores per run
max_tasks=12

################################################################################

## PART 1: creation of potentials

potential_table_folder="target_${model_name}_${symmetry}_${mapping}_contact"

cd $HOME
build_potential "$model_name" "$delta" "$ecc1" "$rad1" "$vEE" "$vEP1" "$vP1P1" "$A" "$n" "$ecc2" "$rad2" "$vEP2" "$vP1P2" "$vP2P2" "$symmetry" "$mapping"
potential_table_path="OSPC_LAMMPS/1-make_potentials/${potential_table_folder}/lammpspot"
cd $HOME

## PART 2-4:
##    - Initialize exec file                                                           A
##    - For each setting of the thermostat
##        - For each state point
##            - Build initial condition                                                       B
##            - For each paraller run
##                - Build folders; put potential and initial condition in folder;             C, D, E
##                - Write input file of the run                                               F
##                - Update exec file with new run                                            G
##    - Close the current exec file                                                          I


init_new_exec "$max_tasks"                                                                   # A

runs_in_this_slrm_file=0
for i in $(seq 0 $((Nalgorithms-1))); do


  set -- ${_thermostat[$i]}

  type=$1
  k=$2
  dump=$3                                                                         # dump=Td if type == NH, dump=gamma if type == LANG

  for j in $(seq 0 $((NstatePoints-1))); do

    set -- ${_params[$j]}

    T=$1
    Nside=$2
    rho=$3

    Ntot=$(compute_Ntot "$lattice" "$Nside")                                    # 4 Nside^3 for FCC, 4 Nside for SC, exit otherwise

    cd $HOME
    init_cond_file=$(build_initial_condition "$symmetry" "$Nside" "$rho" "$ecc1" "$ecc2"	"$Ntot")                 # B
    cd $HOME


    for m in $(seq 1 $Nparal); do

      full_path=$(build_second_set_of_folders "$HOME" "$path0" "$type" "$k" "$dump" "$T" "$Ntot" "$rho" "$m")      # C

      cp -r $potential_table_path $full_path                                                                       # D
      cp "OSPC_LAMMPS/2-startingstate_creators/${init_cond_file}" $full_path                                       # E

      write_input_file "$m" "$full_path" "$init_cond_file" "$dump" "$T" "$k" "$ecc1"                               # F

      update_exec "$m" "$full_path" "$Ncores"                                                             # G

    done
    rm -R "OSPC_LAMMPS/2-startingstate_creators/${init_cond_file}"
  done
done

close_exec

################################################################################
## Launch runs

# chmod u+x runner.sh
# ./runner.sh &
# wait


























#
