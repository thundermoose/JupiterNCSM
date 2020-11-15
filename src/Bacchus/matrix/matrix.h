#ifndef __MATRIX__
#define __MATRIX__

#include <vector/vector.h>
#include <equality_status/equality_status.h>
#include <execution_order/execution_order.h>
#include <combination_table/combination_table.h>
#include <stdio.h>
#include <stdlib.h>

struct _matrix_;
typedef struct _matrix_ *matrix_t;

matrix_t new_zero_matrix(size_t num_rows,
			 size_t num_columns);

matrix_t new_random_symmetric_matrix(size_t side_length);

matrix_t new_matrix_from_numpy(FILE* matrix_file);

matrix_t new_generative_matrix(execution_order_t execution_order,
			       combination_table_t combination_table,
			       const char *index_lists_base_directory,
			       const char *matrix_file_base_directory,
			       size_t maximum_loaded_memory);

void matrix_vector_multiplication(vector_t result_vector,
				  const matrix_t matrix,
				  const vector_t vector);

size_t get_num_rows(matrix_t matrix);

size_t get_num_columns(matrix_t matrix);

void set_matrix_elements(matrix_t matrix,
			 const double *elements);

void set_matrix_element(matrix_t matrix,
			size_t i,
			size_t j,
			double element);

double* get_matrix_elements(matrix_t matrix);

void save_numpy_matrix(FILE *file,
		       const matrix_t matrix);

equality_status_t compare_matrices(const matrix_t first_matrix,
				   const matrix_t second_matrix);

void free_matrix(matrix_t matrix);
#endif
