#!/usr/bin/env python

## LIBRARIES

import sys
import numpy as np


def generate_fcc_lattice(N, L):
    # FCC lattice unit cell points
    unit_cell = np.array([
        [0, 0, 0],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5]
    ])

    if L / ((N / 4) ** (1/3)) <= 1:
        print('Error: The smallest distance between particles in the lattice will be <= 1.')
        raise SystemExit

    num_cells_per_dim = int(np.ceil((N / 4) ** (1/3)))

    lattice = []
    for i in range(num_cells_per_dim):
        for j in range(num_cells_per_dim):
            for k in range(num_cells_per_dim):
                for point in unit_cell:
                    lattice_point = (np.array([i, j, k]) + point) * (L / num_cells_per_dim)
                    lattice.append(lattice_point)

                    if len(lattice) == N:
                        return np.array(lattice)

    return np.array(lattice)


### INPUT

in_file = sys.argv[1]
output_file = in_file.replace('settings_num', 'init_cond_num')

with open(in_file) as F:
    for line in F:

        if line.startswith('L ='):
            L = float(line.split(' ')[2])

        if line.startswith('N ='):
            N = float(line.split(' ')[2])

        if line.startswith('a ='):
            a = float(line.split(' ')[2])


fcc_lattice = generate_fcc_lattice(N, L)

### OUTPUT

F = open(output_file, 'w')

F.write('0 %d %f %f %f\n' %(N, L, L, L))
for site in fcc_lattice:
    rx, ry, rz = site
    F.write('%f %f %f\n' %(rx, ry, rz))
    F.write('%f %f %f\n' %(rx, ry, rz+a))
    F.write('%f %f %f\n' %(rx, ry, rz-a))

F.close()








































#
