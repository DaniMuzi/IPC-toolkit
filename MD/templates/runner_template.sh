#!/bin/bash

################################################################################

function check_and_sleep() {

  max_tasks=$1
  process_name="lmp_mpi"                                                             # use lmp if mpi is not available

  running_tasks=`ps -C ${process_name} --no-headers | wc -l`
  while (( $running_tasks >= $max_tasks )); do
    sleep 60
    running_tasks=`ps -C ${process_name} --no-headers | wc -l`
  done

}

################################################################################

#MAX TASKS
mem_per_task=6GB
LAMMPSpath="your/path/to/lammps"
