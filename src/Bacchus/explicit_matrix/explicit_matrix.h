#ifndef __EXPLICIT_MATRIX__
#define __EXPLICIT_MATRIX__

struct _explicit_matrix_;
typedef struct _explicit_matrix_ *explicit_matrix_t;

explicit_matrix_t new_zero_explicit_matrix(size_t num_rows,
			 size_t num_columns);

explicit_matrix_t new_random_symmetric_explicit_matrix(size_t side_length);

explicit_matrix_t new_explicit_matrix_from_numpy(FILE* explicit_matrix_file);

void explicit_matrix_vector_multiplication(vector_t result_vector,
				  const explicit_matrix_t explicit_matrix,
				  const vector_t vector);

size_t get_explicit_matrix_num_rows(explicit_matrix_t explicit_matrix);

size_t get_explicit_matrix_num_columns(explicit_matrix_t explicit_matrix);

void set_explicit_matrix_elements(explicit_matrix_t explicit_matrix,
			 const double *elements);

void set_explicit_matrix_element(explicit_matrix_t explicit_matrix,
			size_t i,
			size_t j,
			double element);

double* get_explicit_matrix_elements(explicit_matrix_t explicit_matrix);

void save_numpy_explicit_matrix(FILE *file,
		       const explicit_matrix_t explicit_matrix);

equality_status_t 
compare_matrices(const explicit_matrix_t first_explicit_matrix,
		 const explicit_matrix_t second_explicit_matrix);

void free_explicit_matrix(explicit_matrix_t explicit_matrix);
#endif
