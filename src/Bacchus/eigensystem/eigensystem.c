#include <eigensystem/eigensystem.h>
#include <stdio.h>
#include <basis/basis.h>
#include <string.h>
#include <assert.h>

struct _eigensystem_
{
	size_t num_eigenvalues;
	double *eigenvalues;
	double *eigenvector_amplitudes;
	basis_t vector_space_basis;	
};

eigensystem_t new_empty_eigensystem(const size_t num_eigenvalues)
{
	eigensystem_t eigensystem = 
		(eigensystem_t)malloc(sizeof(struct _eigensystem_));
	eigensystem->num_eigenvalues = num_eigenvalues;
	eigensystem->eigenvalues = 
		(double*)calloc(num_eigenvalues,sizeof(double));
	eigensystem->eigenvector_amplitudes =
		(double*)calloc(num_eigenvalues*num_eigenvalues,
				sizeof(double));
	return eigensystem;	
}

void set_basis(eigensystem_t eigensystem,
	       basis_t basis)
{
	eigensystem->vector_space_basis = basis;
}

void set_eigenvalues(eigensystem_t eigensystem,
	double *eigenvalues)
{
	assert(eigensystem);
	assert(eigensystem->eigenvalues);
	memcpy(eigensystem->eigenvalues,
		eigenvalues,
		eigensystem->num_eigenvalues*sizeof(double));
}

void set_raw_eigenvectors(eigensystem_t eigensystem,
			   double *eigenvectors)
{
	memcpy(eigensystem->eigenvector_amplitudes,
       		eigenvectors,
 		eigensystem->num_eigenvalues*eigensystem->num_eigenvalues*
		sizeof(double));		
}

size_t get_num_eigenvalues(eigensystem_t eigensystem)
{
	assert(eigensystem != NULL);
	return eigensystem->num_eigenvalues;
}

double get_eigenvalue(eigensystem_t eigensystem,
		size_t eigenvalue_index)
{
	assert(eigensystem != NULL);
	assert(eigensystem->eigenvalues != NULL);
	assert(eigensystem->num_eigenvalues>eigenvalue_index);
	return eigensystem->eigenvalues[eigenvalue_index];
}

double *get_eigenvector_amplitudes(eigensystem_t eigensystem,
				   size_t eigenvector_index)
{
	return eigensystem->eigenvector_amplitudes + 
		eigenvector_index*eigensystem->num_eigenvalues;
}

double* get_eigenvalues(eigensystem_t eigensystem)
{
	assert(eigensystem != NULL);
	double *eigenvalues = 
		(double*)malloc(eigensystem->num_eigenvalues*sizeof(double));
	memcpy(eigenvalues,
		eigensystem->eigenvalues,
		eigensystem->num_eigenvalues*sizeof(double));
	return eigenvalues;
}

void get_eigenvector(vector_t result,
		      eigensystem_t eigensystem,
		      size_t eigenvector_index)
{
	double *amplitudes = 
		eigensystem->eigenvector_amplitudes +
		eigenvector_index*eigensystem->num_eigenvalues;
	basis_construct_vector(result,
			       eigensystem->vector_space_basis,
			       amplitudes,
			       eigensystem->num_eigenvalues);
}

void print_eigensystem(const eigensystem_t eigensystem)
{
	printf("eigensystem:\n"
	       "{\n");
	printf("\tnum_eigenvalues = %lu;\n",
	       eigensystem->num_eigenvalues);
	printf("\teigenvalues = ( ");
	for (size_t i = 0; i<eigensystem->num_eigenvalues; i++)
		printf("%.17lg%s ",
		       eigensystem->eigenvalues[i],
		       i == eigensystem->num_eigenvalues-1 ? "" : ","); 
	printf(");\n"
	       "}\n");

}

void free_eigensystem(eigensystem_t eigensystem)
{
	assert(eigensystem != NULL);
	free(eigensystem->eigenvalues);
	free(eigensystem);
}
