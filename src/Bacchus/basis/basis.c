#include <basis/basis.h>
#include <assert.h>
#include <string_tools/string_tools.h>

struct _basis_
{
	vector_settings_t vector_settings;
	char *basis_directory;
	vector_t *vectors;
	size_t *num_vectors;
	size_t max_num_vectors;
	int *num_instances;
};

basis_t new_basis_empty(vector_settings_t vector_settings,
			char *basis_directory,
			size_t max_num_vectors)
{
	basis_t basis = (basis_t)malloc(sizeof(struct _basis_));
	basis->vector_settings = vector_settings;
	basis->basis_directory = copy_string(basis_directory);
	basis->vectors = (vector_t*)calloc(max_num_vectors,
					   sizeof(vector_t));
	basis->num_vectors = (size_t*)calloc(1,sizeof(size_t));
	basis->max_num_vectors = max_num_vectors;
	basis->num_instances = (int*)malloc(sizeof(int));
	*basis->num_instances = 1;
	return basis;
}

basis_t new_basis_copy(basis_t origin)
{
	basis_t basis = (basis_t)malloc(sizeof(struct _basis_));
	*basis = *origin;
	(*basis->num_instances)++;
	return basis;
}

void basis_append_vector(basis_t basis)
{
	assert(*basis->num_vectors < basis->max_num_vectors);
	char vector_directory_name[2048] = {0};
	sprintf(vector_directory_name,
		"%s/vector_%lu",
		basis->basis_directory,
		*basis->num_vectors);
	vector_settings_t vector_settings =
		basis->vector_settings;
	vector_settings.directory_name = copy_string(vector_directory_name);
	basis->vectors[*basis->num_vectors] = 
		new_zero_vector(vector_settings);
	(*basis->num_vectors)++;
}

void basis_remove_last(basis_t basis)
{
	(*basis->num_vectors)--;
	free_vector(basis->vectors[*basis->num_vectors]);
}

vector_t basis_get_vector(basis_t basis,
			  size_t vector_index)
{
	return basis->vectors[vector_index];
}

vector_t *basis_get_all_vectors(basis_t basis)
{
	return basis->vectors;
}

size_t basis_get_dimension(basis_t basis)
{
	return *basis->num_vectors;
}

void basis_construct_vector(vector_t result,
			    basis_t basis,
			    double *amplitudes,
			    size_t num_amplitudes)
{
	assert(basis != NULL);
	assert(amplitudes != NULL);
	assert(basis->vectors != NULL);
	assert(num_amplitudes < *basis->num_vectors);	
	for (size_t i = 0; i<num_amplitudes; i++)
		vector_add_scaled(result,
				  amplitudes[i],
				  basis->vectors[i]);
}

void free_basis(basis_t basis)
{
	if (*basis->num_instances == 1)
	{
		for (size_t i = 0; i<*basis->num_vectors; i++)
			free_vector(basis->vectors[i]);
		free(basis->vectors);
		free(basis->basis_directory);
		free(basis->num_vectors);
		free(basis->num_instances);
	}
	else
		(*basis->num_instances)--;
	free(basis);
}
