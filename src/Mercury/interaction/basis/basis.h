#ifndef __BASIS__
#define __BASIS__

#include <stdlib.h>

struct _basis_;
typedef struct _basis_ *basis_t;


basis_t read_basis(const char *interaction_path,
		   size_t num_particles);

size_t find_basis_state(basis_t basis,
		  int *state,
		  size_t num_particles);

void free_basis(basis_t basis);


#endif
