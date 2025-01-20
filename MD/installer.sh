#! /bin/bash

ipc_lammps=$( ls IPC/ )

if [ -z "$ipc_lammps" ]; then
  git submodule init
  git submodule update
fi
