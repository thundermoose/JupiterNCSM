#include <matrix/matrix.h>
#include <explicit_matrix/explicit_matrix.h>
#include <log/log.h>
#include <error/error.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <thundertester/test.h>
#include <scheduler/scheduler.h>
#include <math.h>

typedef enum
{
	EXPLICIT_MATRIX,
	GENERATIV_MATRIX
} matrix_type_t;
struct _matrix_
{
	matrix_type_t type;
	explicit_matrix_t explicit_matrix;
	scheduler_t scheduler;
};

matrix_t new_zero_matrix(size_t num_rows,
			 size_t num_columns)
{
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->type = EXPLICIT_MATRIX;
	matrix->explicit_matrix = new_zero_explicit_matrix(num_rows,
							   num_columns);
	return matrix;
}

matrix_t new_random_symmetric_matrix(size_t side_length)
{
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->type = EXPLICIT_MATRIX;
	matrix->explicit_matrix = 
		new_random_symmetric_explicit_matrix(side_length);
	return matrix;
}

matrix_t new_matrix_from_numpy(FILE* matrix_file)
{
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->type = EXPLICIT_MATRIX;
	matrix->explicit_matrix = new_explicit_matrix_from_numpy(matrix_file);
	return matrix;
}

matrix_t new_generative_matrix(execution_order_t execution_order,
			       combination_table_t combination_table,
			       const char *index_lists_base_directory,
			       const char *matrix_file_base_directory)
{
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->type = GENERATIV_MATRIX;
	matrix->scheduler = new_scheduler(execution_order,
					  combination_table,
					  index_lists_base_directory,
					  matrix_file_base_directory);
	return matrix;
}

void matrix_vector_multiplication(vector_t result_vector,
				  const matrix_t matrix,
				  const vector_t vector)
{
	switch (matrix->type)
	{
		case EXPLICIT_MATRIX:
			explicit_matrix_vector_multiplication
				(result_vector,
				 matrix->explicit_matrix,
				 vector);
			break;
		case GENERATIV_MATRIX:
			run_matrix_vector_multiplication
				(get_vector_path(result_vector),
				 get_vector_path(vector),
				 matrix->scheduler);
			break;
	}
}

size_t get_num_rows(matrix_t matrix)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	return get_explicit_matrix_num_rows(matrix->explicit_matrix);
}

size_t get_num_columns(matrix_t matrix)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	return get_explicit_matrix_num_columns(matrix->explicit_matrix);
}

void set_matrix_elements(matrix_t matrix,
			 const double *elements)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	set_explicit_matrix_elements(matrix->explicit_matrix,
				     elements);
}

void set_matrix_element(matrix_t matrix,
			const size_t i,
			const size_t j,
			const double element)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	set_explicit_matrix_element(matrix->explicit_matrix,
				    i,j,element);
}

double *get_matrix_elements(matrix_t matrix)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	return get_explicit_matrix_elements(matrix->explicit_matrix);
}

void save_numpy_matrix(FILE *file,
		       const matrix_t matrix)
{
	assert(matrix->type == EXPLICIT_MATRIX);
	save_numpy_explicit_matrix(file,matrix->explicit_matrix);	
}

equality_status_t compare_matrices(const matrix_t first_matrix,
				   const matrix_t second_matrix)
{
	if (first_matrix->type != second_matrix->type)
		return NOT_EQUAL;
	switch (first_matrix->type)
	{
		case EXPLICIT_MATRIX:
			return compare_explicit_matrices
				(first_matrix->explicit_matrix,
				 second_matrix->explicit_matrix);
		case GENERATIV_MATRIX:
			if (first_matrix->scheduler == 
				second_matrix->scheduler)
				return EQUAL;
			else
				return NOT_EQUAL;
	}
	return EQUAL;
}

void free_matrix(matrix_t matrix)
{
	switch(matrix->type)
	{
		case EXPLICIT_MATRIX:
			free_explicit_matrix(matrix->explicit_matrix);
			break;
		case GENERATIV_MATRIX:
			free_scheduler(matrix->scheduler);
			break;
	}
	free(matrix);
}
new_test(save_to_numpy_file,
	 const char *filename = get_test_file_path("matrix.npy");
	 matrix_t matrix = new_random_symmetric_matrix(10);
	 FILE *outfile = fopen(filename,"w");
	 save_numpy_matrix(outfile,matrix);
	 fclose(outfile);
	 FILE *infile = fopen(filename,"r");
	 matrix_t loaded_matrix = new_matrix_from_numpy(infile);
	 fclose(infile);
	 assert_that(compare_matrices(matrix,loaded_matrix) == EQUAL);
	 free_matrix(matrix);
	 free_matrix(loaded_matrix);
	);
