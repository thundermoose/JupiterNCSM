#ifndef __BASIS__
#define __BASIS__

#include <vector/vector.h>
#include <stdlib.h>

struct _basis_;
typedef struct _basis_ *basis_t;

basis_t new_basis_empty(vector_settings_t vector_settings,
			char *basis_directory,
			size_t max_num_vectors);

basis_t new_basis_copy(basis_t origin);

void basis_append_vector(basis_t basis);

void basis_remove_last(basis_t basis);

vector_t basis_get_vector(basis_t basis,
			  size_t vector_index);

vector_t *basis_get_all_vectors(basis_t basis);

size_t basis_get_dimension(basis_t basis);

void basis_construct_vector(vector_t result,
			    basis_t basis,
			    double *amplitudes,
			    size_t num_amplitudes);

void free_basis(basis_t basis);

#endif
