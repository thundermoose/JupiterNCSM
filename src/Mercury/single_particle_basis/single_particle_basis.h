#ifndef __SINGLE_PARTICLE_BASIS__
#define __SINGLE_PARTICLE_BASIS__

#include <stdlib.h>

struct _single_particle_basis_;
typedef struct _single_particle_basis_ *single_particle_basis_t;

typedef struct
{
	int n,l,j,m,tz;
	size_t neptune_index;
} single_particle_state_t;

single_particle_basis_t new_single_particle_basis(int energy_max);

single_particle_state_t get_state(single_particle_basis_t basis,
				  size_t index);

void free_single_particle_basis(single_particle_basis_t basis);

int get_m(single_particle_state_t state);

#endif
