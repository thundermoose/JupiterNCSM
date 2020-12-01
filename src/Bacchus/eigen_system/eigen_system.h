#ifndef __EIGEN_SYSTEM__
#define __EIGEN_SYSTEM__

#include <stdlib.h>
#include <basis/basis.h>

struct _eigen_system_;
typedef struct _eigen_system_ *eigen_system_t;

eigen_system_t new_empty_eigensystem(const size_t num_eigen_values);

void set_basis(eigen_system_t eigen_system,
	       basis_t basis);

void set_eigen_values(eigen_system_t eigen_system,
	double *eigen_values);

void set_raw_eigen_vectors(eigen_system_t eigen_system,
			   double *eigen_vectors);

size_t get_num_eigen_values(eigen_system_t eigen_system);

double get_eigen_value(eigen_system_t eigen_system,
		size_t eigen_value_index);

double* get_eigen_values(eigen_system_t eigen_system);

void  get_eigen_vector(vector_t result,
		       eigen_system_t eigen_system,
		       size_t eigen_vector_index);

void print_eigen_system(const eigen_system_t eigen_system);

void free_eigen_system(eigen_system_t eigen_system);
#endif
