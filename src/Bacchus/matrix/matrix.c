#include <matrix/matrix.h>
#include <numpy_array/numpy_array.h>
#include <log/log.h>
#include <error/error.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <thundertester/test.h>
#include <math.h>

struct _matrix_
{
	size_t num_rows;
	size_t num_columns;
	double *elements;
};

static
void write_numpy_header(FILE* file,
			size_t num_rows,
			size_t num_columns);

static
size_t append_numpy_padding(char *buffer,
			    size_t message_length);

matrix_t new_zero_matrix(size_t num_rows,
			 size_t num_columns)
{
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->num_rows = num_rows;
	matrix->num_columns = num_columns;
	matrix->elements = (double*)calloc(num_rows*num_columns,
					   sizeof(double));
	return matrix;
}

matrix_t new_random_symmetric_matrix(size_t side_length)
{
	matrix_t matrix = new_zero_matrix(
					  side_length,
					  side_length);
	for (size_t i = 0; i < side_length; i++)
		for (size_t j = i; j < side_length; j++)
			matrix->elements[j*side_length+i] = 
				matrix->elements[i*side_length+j] =
				2*drand48()-1;
	return matrix;
}

matrix_t new_matrix_from_numpy(FILE* matrix_file)
{
	numpy_array_t numpy_array = read_numpy_file(matrix_file);
	const numpy_header_t header = get_numpy_header(numpy_array);
	double *elements = get_numpy_array_elements(numpy_array);
	free_numpy_array(numpy_array);
	matrix_t matrix = (matrix_t)malloc(sizeof(struct _matrix_));
	matrix->num_rows = header.num_rows;
	matrix->num_columns = header.num_columns;
	matrix->elements = elements;
	return matrix;
}

void matrix_vector_multiplication(vector_t result_vector,
				  const matrix_t matrix,
				  const vector_t vector)
{
	assert(matrix != NULL);
	assert(matrix->elements != NULL);
	assert(vector != NULL);
	assert(matrix->num_columns == vector_dimension(vector));
	assert(matrix->num_rows == vector_dimension(result_vector));
	log_entry("matrix_vector_multiplication");
	for (size_t i = 0; i < matrix->num_rows; i++)
	{
		double accumulator = 0.0;
		for (size_t j = 0; j < matrix->num_columns; j++)
		{
			size_t I = i*matrix->num_columns + j;
			double element = matrix->elements[I];
			double v_element = get_element(vector,j);
			log_entry("element = %lg",element);
			log_entry("v_element = %lg",v_element);
			accumulator += v_element*element;
		}
		log_entry("accumulator = %lg",accumulator);
		set_element(result_vector,i,accumulator);
	}
	save_vector(result_vector);
}

size_t get_num_rows(matrix_t matrix)
{
	return matrix->num_rows;
}

size_t get_num_columns(matrix_t matrix)
{
	return matrix->num_columns;
}

void set_matrix_elements(matrix_t matrix,
			 const double *elements)
{
	assert(matrix != NULL);
	assert(matrix->elements != NULL);
	assert(elements != NULL);
	memcpy(matrix->elements,
	       elements,
	       sizeof(double)*matrix->num_rows*matrix->num_columns);
}

void set_matrix_element(matrix_t matrix,
			const size_t i,
			const size_t j,
			const double element)
{
	size_t I = i*matrix->num_columns + j;
	matrix->elements[I] = element;
}

double *get_matrix_elements(matrix_t matrix)
{
	assert(matrix != NULL);
	assert(matrix->elements != NULL);
	double *elements = 
		(double*)malloc(matrix->num_rows*matrix->num_columns*
				sizeof(double));
	memcpy(elements,
	       matrix->elements,
	       matrix->num_rows*matrix->num_columns*
	       sizeof(double));
	return elements;
}

void save_numpy_matrix(FILE *file,
		       const matrix_t matrix)
{
	write_numpy_header(file,matrix->num_rows,matrix->num_columns);
	const size_t num_elements = matrix->num_rows*matrix->num_columns;
	if (fwrite(matrix->elements,
		   sizeof(double),
		   num_elements,
		   file) < num_elements)
		error("Could not write elements to numpy file\n");
}

equality_status_t compare_matrices(const matrix_t first_matrix,
				   const matrix_t second_matrix)
{
	if (first_matrix->num_rows != second_matrix->num_rows ||
	    first_matrix->num_columns != second_matrix->num_columns)
		return NOT_EQUAL;
	const size_t num_elements =
	       	first_matrix->num_rows*first_matrix->num_columns;
	const double tolerance = 1e-10;
	for (size_t i = 0; i < num_elements; i++)
		if (fabs(first_matrix->elements[i] -
			 second_matrix->elements[i])>tolerance)
			return NOT_EQUAL;
	return EQUAL;
}

void free_matrix(matrix_t matrix)
{
	free(matrix->elements);
	free(matrix);
}

static
void write_numpy_header(FILE* file,
			size_t num_rows,
			size_t num_columns)
{
	char buffer[2048];
	size_t message_length = sprintf(buffer,
		"\x93NUMPY\x01%cF%c{'descr': '<f8','fortran_order': False, 'shape': (%ld,%ld), }",
		0,0,num_rows,num_columns);
	size_t padding_length = append_numpy_padding(buffer,message_length);
	uint64_t header_end = 0x0a20202020202020;
	memcpy(buffer+message_length+padding_length,
		&header_end,
		sizeof(uint64_t));
	size_t total_length =
	       	message_length + padding_length + sizeof(uint64_t);
	if (fwrite(buffer,
		   sizeof(char),
		   total_length,
		   file) < total_length)
		error("Could not write header to numpy file\n");
}

static
size_t append_numpy_padding(char *buffer,
			    size_t message_length)
{
	const size_t padding_length = 8 - message_length % 8;
	for (size_t i = 0; i<padding_length; i++)
		buffer[i+message_length] = 0x20;
	return padding_length;
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
