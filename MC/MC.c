#include "MC.h"
#include "vmmc.h"
#include "utils.h"
#include "basic_functions.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


void MC_init(input_file *input, System *syst, Output *IO) {

	// Here we set the pointer to the function that will be used to make a Monte Carlo step
	// according to the ensemble specified in the input file
	switch(syst->ensemble) {

		case NVT:
				syst->do_ensemble = &do_NVT;
				break;
		default:
			output_exit(IO, "Ensemble %d not supported\n", syst->ensemble);
			break;

	}


	// Here we set the pointer to the function that will be used to perform the type of dynamics
	// specified in the input file
	switch(syst->dynamics) {
		case VMMC:
			VMMC_init(input, syst, IO);
			syst->do_dynamics = &VMMC_dynamics;
			break;
		default:
			output_exit(IO, "Dynamics %d not supported\n", syst->dynamics);
			break;
	}

	// Here we set the pointer to the function that computes the energy weights
	// as specified in the input file
	switch (syst->energy_weights) {
		case SPHERES:
		syst->energetic_weights = &overlap_volume;
		break;
	case EXPO:
		syst->energetic_weights = &exponential_weight;
		break;
	default:
		output_exit(IO, "Energy weights %d not supported\n", syst->dynamics);
		break;
	}

}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Ensembles


void do_NVT(System *syst, Output *output_files) {
	int i;
	for(i = 0; i < syst->N; i++) {
		syst->do_dynamics(syst, output_files);
	}
}


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Energy calculation

double MC_energy(System *syst, PatchyParticle *p) {

	vector O;
	double E = 0.;
	int j, k, l, val;
	int ind[3], loop_ind[3], loop_index;

	syst->overlap = 0;

	cells_fill_and_get_idx_from_particle(syst, p, ind);

	for(j = -1; j < 2; j++) {                                                                          // Iterate over first dimension of the box
		loop_ind[0] = (ind[0] + j + syst->cells->N_side[0]) % syst->cells->N_side[0];
		for(k = -1; k < 2; k++) {                                                                        // second dimension
			loop_ind[1] = (ind[1] + k + syst->cells->N_side[1]) % syst->cells->N_side[1];
			for(l = -1; l < 2; l++) {                                                                      // Third dimension
				loop_ind[2] = (ind[2] + l + syst->cells->N_side[2]) % syst->cells->N_side[2];
				loop_index = (loop_ind[0] * syst->cells->N_side[1] + loop_ind[1]) * syst->cells->N_side[2] + loop_ind[2];         // Index of a cell relevant to particle p

				PatchyParticle *q = syst->cells->heads[loop_index];                     // One particle that is in the cell with index loop_index (108)
				while(q != NULL) {                                                      // At second stage q = 32
					if(q->index != p->index) {
						val = MC_interact(syst, p, q, O, syst->K);

						if(val == BOND) {
							E += SCALAR(syst->e, O);
						}
						else if(val == OVERLAP) {
							syst->overlap = 1;
							return 1e8;
						}
					}
					q = syst->cells->next[q->index];                                      // next particle in the cell with index loop_index (ex: next[108] = 32)
				}                                                                       // At the end of the second stage, next[32] = NULL -> break the loop
			}
		}
	}

	return E;
}


inline int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K) {
	return (*syst->energetic_weights)(syst, p, q, O, K);
}


int overlap_volume(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K){

	int i, j;
	double R, r, d;

	R = syst->sigma_c + syst->delta_c;
	r = syst->sigma_p;

	vector rab = {q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);
	d = sqrt(SCALAR(rab, rab));

	if (d - 2.0*syst->sigma_c < -1e-12) return OVERLAP;
	if (2.0*R - d < -1e-12) return NO_BOND;                                              // This MUST be changed if the IPC constraint is OFF!!!


	double O_EE, O_EP, O_PP, norm;

	// EE
	O_EE = 0;
	O_EE += bonding_volume(syst, R, R, p->r, q->r);

	// EP
	O_EP = 0;
	for(i=0; i<syst->n_patches; i++) {
		O_EP += bonding_volume(syst, R, r, p->r, q->patches[i]);                           // E1, P_i
		O_EP += bonding_volume(syst, r, R, p->patches[i], q->r);                           // E2, P_i
	}


	// PP
	O_PP = 0;
	for (i=0; i<syst->n_patches; i++) {
		for (j=0; j<syst->n_patches; j++) {
			O_PP += bonding_volume(syst, r, r, p->patches[i], q->patches[j]);                 // P_i, P_j
		}
	}

	norm = (4.0/3.0) * M_PI * (syst->sigma_c*syst->sigma_c*syst->sigma_c);

	O[0] = O_EE / norm;
	O[1] = O_EP / norm;
	O[2] = O_PP / norm;

	return BOND;

}

double bonding_volume(System *syst, double Ra, double Rb, vector ra, vector rb) {

	double Rm, rmax, rmin, d, f;

	vector rab = {ra[0] - rb[0], ra[1] - rb[1], ra[2] - rb[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);

	d = sqrt(SCALAR(rab, rab));

	Rm = (Ra <= Rb) ? Ra : Rb;
	rmax = Ra + Rb;
	rmin = fabs(Ra - Rb);

	if (d >= rmax) { f = 0; }
	else if (d <= rmin) { f = (4.0/3.0)*M_PI*Rm*Rm*Rm; }
	else {                                                                        	// THAT IS: else if (rmin <= d && d <= rmax){

		double f1, f2, f3 ,f4;

		f1 = 2*Ra + (SQR(Ra) - SQR(Rb) + SQR(d)) / (2*d);
		f2 = Ra - (SQR(Ra) - SQR(Rb) + SQR(d)) / (2*d);

		f3 = 2*Rb - (SQR(Ra) - SQR(Rb) - SQR(d)) / (2*d);
		f4 = Rb + (SQR(Ra) - SQR(Rb) - SQR(d)) / (2*d);

		f = (f1*f2*f2 + f3*f4*f4) * M_PI / 3.0;

	}

	return f;


}

int exponential_weight(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K){

	int i, j;
	double d;

	vector rab = {q->r[0] - p->r[0], q->r[1] - p->r[1], q->r[2] - p->r[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);
	d = sqrt(SCALAR(rab, rab));

	if (d - 2.0*syst->sigma_c < -1e-12) return OVERLAP;
	if (syst->r_cut - d < -1e-12) return NO_BOND;                                              // This MUST be changed if the IPC constraint is OFF!!!


	double ree = 2.0 * syst->sigma_c;
	double rep = 2.0 * syst->sigma_c - syst->a;
	double rpp = 2.0 * (syst->sigma_c - syst->a);

	double O_EE, O_EP, O_PP;


	// EE
	O_EE = 0;
	O_EE += single_expo_weight(syst, p->r, q->r, ree, syst->K[0]);

	// EP
	O_EP = 0;
	for(i=0; i<syst->n_patches; i++) {
		O_EP += single_expo_weight(syst, p->patches[i], q->r, rep, syst->K[1]);                           // E1, P_i
		O_EP += single_expo_weight(syst, p->r, q->patches[i], rep, syst->K[1]);                           // E2, P_i
	}

	// PP
	O_PP = 0;
	for (i=0; i<syst->n_patches; i++) {
		for (j=0; j<syst->n_patches; j++) {
			O_PP += single_expo_weight(syst, p->patches[i], q->patches[j], rpp, syst->K[2]);                 // P_i, P_j
		}
	}


	O[0] = O_EE;
	O[1] = O_EP;
	O[2] = O_PP;


	return BOND;

}

double single_expo_weight(System *syst, vector ra, vector rb, double rc, double k) {

	double d;

	vector rab = {ra[0] - rb[0], ra[1] - rb[1], ra[2] - rb[2]};

	rab[0] -= syst->box[0] * rint(rab[0] / syst->box[0]);
	rab[1] -= syst->box[1] * rint(rab[1] / syst->box[1]);
	rab[2] -= syst->box[2] * rint(rab[2] / syst->box[2]);

	d = sqrt(SCALAR(rab, rab));

	return exp(-k * (d-rc));
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// MC moves


void MC_change_cell(System *syst, PatchyParticle *p) {

	int ind[3];
	int cell_index = cells_fill_and_get_idx_from_particle(syst, p, ind);
	if(cell_index == p->cell) {
		p->cell_old = p->cell;
		return;
	}

	// Remove the particle from the old cell
	Cells *cells = syst->cells;
	PatchyParticle *previous = NULL;
	PatchyParticle *current = cells->heads[p->cell];
	assert(current != NULL);
	while(current->index != p->index) {
		previous = current;
		current = cells->next[current->index];
		assert(cells->next[previous->index]->index == current->index);
	}
	if(previous == NULL) cells->heads[p->cell] = cells->next[p->index];
	else cells->next[previous->index] = cells->next[p->index];

	// Add the particle to the new cell
	cells->next[p->index] = cells->heads[cell_index];
	cells->heads[cell_index] = p;
	p->cell_old = p->cell;
	p->cell = cell_index;
}



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


void MC_check_energy(System *syst, Output *IO, llint t) {

	int i;
	double E = 0.;
	syst->overlap = 0;
	for(i = 0; i < syst->N; i++) {
		PatchyParticle *p = syst->particles + i;
		E += MC_energy(syst, p);
	}
	E *= 0.5;

	if(syst->overlap) output_exit(IO, "Computing energy from scratch resulted in an overlap, aborting\n");
	if(fabs(syst->energy) > 1e-5 && fabs((E - syst->energy) / syst->energy) > 1e-5) {
		// output_exit(IO, "\nEnergy check failed, old energy = %lf, new energy = %lf. This happens when N = %d and t=%lld.\n", syst->energy, E, syst->N, t);
		output_log_msg(IO, "Energy check failed, old energy = %lf, new energy = %lf. This happens when N = %d and t = %lld.\n", syst->energy, E, syst->N, t);
	}

	syst->energy = E;
}



void MC_free(System *syst) {
	switch(syst->dynamics) {
	case VMMC:
		VMMC_free();
		break;
	}
}
