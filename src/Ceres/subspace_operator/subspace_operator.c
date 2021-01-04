#include <subspace_operator/subspace_operator.h>
#include <matrix/matrix.h>


void create_subspace_operator(const char *subspace_operator_path,
			      const char *operator_path,
			      const char *workspace_path,
			      const char *index_list_base_directory,
			      combination_table_t combination_table,
			      evaluation_order_t evaluation_order,
			      vector_settings_t vector_setting,
			      vector_t *training_vectors,
			      size_t num_training_vectors,
			      size_t maximum_loaded_memory)
{
	matrix_t operator = new_generative_matrix(combination_table,
						  evaluation_order,
						  index_list_base_directory,
						  operator_path,
						  maximum_loaded_memory);
	vector_t *intermediate_vectors = 
		(vector_t*)malloc(num_training_vectors*sizeof(vector_t));
	char *vector_path_buffer =
	       	calloc(strlen(workspace_path)+128,sizeof(char));
	for (size_t i = 0; i < num_training_vectors; i++)
	{
		sprintf(vector_path_buffer,
			"%s/intermediate_vector_%lu\n",
			workspace_path,
			i);
		vector_setting.directory_name = vector_path_buffer;
		intermediate_vectors[i] = new_zero_vector(vector_setting);
		matrix_vector_multiplication(intermediate_vectors[i],
					     operator,
					     training_vectors[i]);
	}
	free(vector_path_buffer);
	double *matrix_elements = (double*)malloc(num_training_vectors*
						  num_training_vectors*
						  sizeof(double));	
	for (size_t i = 0; i<num_training_vectors; i++)
	{
		for (size_t j = i; j<num_training_vectors; i++)
		{
			double element =
				scalar_multiplication(training_vectors[j],
						      intermediate_vectors[i]);	      
			matrix_elements[i*num_training_vectors+j] = element;
			matrix_elements[j*num_training_vectors+i] = element;
		}
		free_vector(intermediate_vectors[i]);
	}
	free_matrix(operator);
	save_as_numpy_matrix(subspace_operator_path,
			     matrix_elements,
			     num_training_vectors,
			     num_training_vectors);
	free(elements);
}
