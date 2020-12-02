#ifndef __EIGEN_SYSTEM__
#define __EIGEN_SYSTEM__

#include <stdlib.h>
#include <basis/basis.h>

struct _eigensystem_;
typedef struct _eigensystem_ *eigensystem_t;

eigensystem_t new_empty_eigensystem(const size_t num_eigenvalues);

void set_basis(eigensystem_t eigensystem,
	       basis_t basis);

void set_eigenvalues(eigensystem_t eigensystem,
	double *eigenvalues);

void set_raw_eigenvectors(eigensystem_t eigensystem,
			   double *eigenvectors);

size_t get_num_eigenvalues(eigensystem_t eigensystem);

double get_eigenvalue(eigensystem_t eigensystem,
		size_t eigenvalue_index);

double *get_eigenvector_amplitudes(eigensystem_t eigensystem,
				    size_t eigenvector_index);

double* get_eigenvalues(eigensystem_t eigensystem);

void  get_eigenvector(vector_t result,
		       eigensystem_t eigensystem,
		       size_t eigenvector_index);

void print_eigensystem(const eigensystem_t eigensystem);

void free_eigensystem(eigensystem_t eigensystem);
#endif
