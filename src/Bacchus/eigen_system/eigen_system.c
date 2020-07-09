#include <eigen_system/eigen_system.h>
#include <basis/basis.h>
#include <string.h>
#include <assert.h>

struct _eigen_system_
{
	size_t num_eigen_values;
	double *eigen_values;
	double *eigen_vector_amplitudes;
	basis_t vector_space_basis;	
};

eigen_system_t new_empty_eigensystem(const size_t num_eigen_values)
{
	eigen_system_t eigen_system = 
		(eigen_system_t)malloc(sizeof(struct _eigen_system_));
	eigen_system->num_eigen_values = num_eigen_values;
	eigen_system->eigen_values = 
		(double*)calloc(num_eigen_values,sizeof(double));
	return eigen_system;	
}

size_t get_num_eigen_values(eigen_system_t eigen_system)
{
	assert(eigen_system != NULL);
	return eigen_system->num_eigen_values;
}

double get_eigen_value(eigen_system_t eigen_system,
		size_t eigen_value_index)
{
	assert(eigen_system != NULL);
	assert(eigen_system->eigen_values != NULL);
	assert(eigen_system->num_eigen_values>eigen_value_index);
	return eigen_system->eigen_values[eigen_value_index];
}

double* get_eigen_values(eigen_system_t eigen_system)
{
	assert(eigen_system != NULL);
	double *eigen_values = 
		(double*)malloc(eigen_system->num_eigen_values*sizeof(double));
	memcpy(eigen_values,
		eigen_system->eigen_values,
		eigen_system->num_eigen_values*sizeof(double));
	return eigen_values;
}

void set_eigen_values(eigen_system_t eigen_system,
	double *eigen_values)
{
	assert(eigen_system);
	assert(eigen_system->eigen_values);
	memcpy(eigen_system->eigen_values,
		eigen_values,
		eigen_system->num_eigen_values*sizeof(double));
}

void free_eigen_system(eigen_system_t eigen_system)
{
	assert(eigen_system != NULL);
	free(eigen_system->eigen_values);
	free(eigen_system);
}
