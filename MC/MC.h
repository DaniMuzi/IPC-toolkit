#ifndef MC_H_
#define MC_H_

#include "defs.h"

typedef struct input_file input_file;


void MC_init(input_file *input, System *syst, Output *IO);                 // This exists in the MC.c file
void do_NVT(System *syst, Output *output_files);

double MC_energy(System *syst, PatchyParticle *p);
int MC_interact(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K);

int overlap_volume(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K);
double bonding_volume(System *syst, double Ra, double Rb, vector ra, vector rb);
int exponential_weight(System *syst, PatchyParticle *p, PatchyParticle *q, vector O, vector K);
double single_expo_weight(System *syst, vector ra, vector rb, double rc, double k);

void MC_change_cell(System *syst, PatchyParticle *p);


void MC_check_energy(System *syst, Output *output_files, llint t);
void MC_free(System *syst);


#endif /* MC_H_ */
