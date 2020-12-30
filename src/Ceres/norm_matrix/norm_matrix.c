#include <norm_matrix/norm_matrix.h>
#include <numpy_matrix/numpy_matrix.h>

void create_norm_matrix(const char *norm_matrix_path,
			vector_t *training_vectors,
			size_t num_training_vectors)
{
	double *matrix_elements = (double*)malloc(num_training_vectors*
						  num_training_vectors*
						  sizeof(double));	

	for (size_t i = 0; i<num_training_vectors; i++)
	{
		for (size_t j = i; j<num_training_vectors; j++)
		{
			double element = 
				scalar_multiplication(training_vectors[i],
						      training_vectors[j]);
			matrix_elements[i*num_training_vectors+j] = element;
			matrix_elements[j*num_training_vectors+i] = element;
		}
	}
	save_as_numpy_matrix(norm_matrix_path,
			     matrix_elements,
			     num_training_vectors,
			     num_training_vectors);
	free(matrix_elements);
}
