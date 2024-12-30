#ifndef VMMC_H_
#define VMMC_H_

#define VMMC_ROTATION (1)
#define VMMC_TRANSLATION (2)

#include "defs.h"

typedef struct input_file input_file;

typedef struct vmmc_d {
	double max_move;
	int max_cluster;
	PatchyParticle ** possible_links;
	int n_possible_links;
	PatchyParticle ** prelinked_particles;
	int n_prelinked_particles;
	PatchyParticle ** clust;
	int n_clust;
	int * is_in_cluster;
	int which_move;
	matrix rotation;
} vmmc_d;


void VMMC_init(input_file *input, System *syst, Output *IO);

void VMMC_dynamics(System *syst, Output *IO);

void _populate_possible_links(System *syst, Output *output_files, PatchyParticle *p);
void _move_particle(System * syst, PatchyParticle * p, vector move, double t);
double _compute_cluster_energy(System *syst);
double _pair_energy(System * syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K);

int _mycomp(const void *p, const void * q, void * s);
void _store_dof(PatchyParticle * p);
void _restore_dof(PatchyParticle * p);
void normalized_diff_vector(System *syst, vector a, vector b, vector c);
// void VMMC_reset_counters();
// void print_VMMC_counters(Output *output_files);
void VMMC_free();

#endif
