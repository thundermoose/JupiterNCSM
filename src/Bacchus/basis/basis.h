#ifndef __BASIS__
#define __BASIS__

#include <vector/vector.h>

struct _basis_;
typedef struct _basis_ *basis_t;

// Warning, for the moment the vector pointer is stored
// inside the basis_t structure, this means that it is 
// not safe to free the pointer vectors. Further more
// the free_basis function, at the end, will free this pointer
// this will however be changed in the future
basis_t new_basis_from_vectors(
		vector_t *vectors,
		size_t num_vectors);

void construct_vector(vector_t result,
		      basis_t basis,
		      double *amplitudes,
		      size_t num_amplitudes);

void free_basis(basis_t basis);

#endif
