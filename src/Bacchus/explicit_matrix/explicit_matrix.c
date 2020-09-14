#include <explicit_matrix/explicit_matrix.h>
#include <numpy_array/numpy_array.h>
#include <log/log.h>
#include <error/error.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <unit_testing/test.h>
#include <math.h>

struct _explicit_matrix_
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

explicit_matrix_t new_zero_explicit_matrix(size_t num_rows,
			 size_t num_columns)
{
	explicit_matrix_t explicit_matrix =
	       	(explicit_matrix_t)malloc(sizeof(struct _explicit_matrix_));
	explicit_matrix->num_rows = num_rows;
	explicit_matrix->num_columns = num_columns;
	explicit_matrix->elements = (double*)calloc(num_rows*num_columns,
					   sizeof(double));
	return explicit_matrix;
}

explicit_matrix_t new_random_symmetric_explicit_matrix(size_t side_length)
{
	explicit_matrix_t explicit_matrix = new_zero_explicit_matrix(
					  side_length,
					  side_length);
	for (size_t i = 0; i < side_length; i++)
		for (size_t j = i; j < side_length; j++)
			explicit_matrix->elements[j*side_length+i] = 
				explicit_matrix->elements[i*side_length+j] =
				2*drand48()-1;
	return explicit_matrix;
}

explicit_matrix_t new_explicit_matrix_from_numpy(FILE* explicit_matrix_file)
{
	numpy_array_t numpy_array = read_numpy_file(explicit_matrix_file);
	const numpy_header_t header = get_numpy_header(numpy_array);
	double *elements = get_numpy_array_elements(numpy_array);
	free_numpy_array(numpy_array);
	explicit_matrix_t explicit_matrix =
	       	(explicit_matrix_t)malloc(sizeof(struct _explicit_matrix_));
	explicit_matrix->num_rows = header.num_rows;
	explicit_matrix->num_columns = header.num_columns;
	explicit_matrix->elements = elements;
	return explicit_matrix;
}

void explicit_matrix_vector_multiplication(vector_t result_vector,
				  const explicit_matrix_t explicit_matrix,
				  const vector_t vector)
{
	assert(explicit_matrix != NULL);
	assert(explicit_matrix->elements != NULL);
	assert(vector != NULL);
	assert(explicit_matrix->num_columns == vector_dimension(vector));
	assert(explicit_matrix->num_rows == vector_dimension(result_vector));
	log_entry("explicit_matrix_vector_multiplication");
	for (size_t i = 0; i < explicit_matrix->num_rows; i++)
	{
		double accumulator = 0.0;
		for (size_t j = 0; j < explicit_matrix->num_columns; j++)
		{
			size_t I = i*explicit_matrix->num_columns + j;
			double element = explicit_matrix->elements[I];
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

size_t get_explicit_matrix_num_rows(explicit_matrix_t explicit_matrix)
{
	return explicit_matrix->num_rows;
}

size_t get_explicit_matrix_num_columns(explicit_matrix_t explicit_matrix)
{
	return explicit_matrix->num_columns;
}

void set_explicit_matrix_elements(explicit_matrix_t explicit_matrix,
			 const double *elements)
{
	assert(explicit_matrix != NULL);
	assert(explicit_matrix->elements != NULL);
	assert(elements != NULL);
	memcpy(explicit_matrix->elements,
	       elements,
	       sizeof(double)*explicit_matrix->num_rows*
	       explicit_matrix->num_columns);
}

void set_explicit_matrix_element(explicit_matrix_t explicit_matrix,
			const size_t i,
			const size_t j,
			const double element)
{
	size_t I = i*explicit_matrix->num_columns + j;
	explicit_matrix->elements[I] = element;
}

double *get_explicit_matrix_elements(explicit_matrix_t explicit_matrix)
{
	assert(explicit_matrix != NULL);
	assert(explicit_matrix->elements != NULL);
	double *elements = 
		(double*)malloc(explicit_matrix->num_rows*explicit_matrix->num_columns*
				sizeof(double));
	memcpy(elements,
	       explicit_matrix->elements,
	       explicit_matrix->num_rows*explicit_matrix->num_columns*
	       sizeof(double));
	return elements;
}

void save_numpy_explicit_matrix(FILE *file,
		       const explicit_matrix_t explicit_matrix)
{
	write_numpy_header(file,explicit_matrix->num_rows,explicit_matrix->num_columns);
	const size_t num_elements = explicit_matrix->num_rows*explicit_matrix->num_columns;
	if (fwrite(explicit_matrix->elements,
		   sizeof(double),
		   num_elements,
		   file) < num_elements)
		error("Could not write elements to numpy file\n");
}

equality_status_t 
compare_explicit_matrices(const explicit_matrix_t first_explicit_matrix,
			  const explicit_matrix_t second_explicit_matrix)
{
	if (first_explicit_matrix->num_rows != second_explicit_matrix->num_rows ||
	    first_explicit_matrix->num_columns != second_explicit_matrix->num_columns)
		return NOT_EQUAL;
	const size_t num_elements =
	       	first_explicit_matrix->num_rows*first_explicit_matrix->num_columns;
	const double tolerance = 1e-10;
	for (size_t i = 0; i < num_elements; i++)
		if (fabs(first_explicit_matrix->elements[i] -
			 second_explicit_matrix->elements[i])>tolerance)
			return NOT_EQUAL;
	return EQUAL;
}

void free_explicit_matrix(explicit_matrix_t explicit_matrix)
{
	free(explicit_matrix->elements);
	free(explicit_matrix);
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

new_test(save_explicit_matrix_to_numpy_file,
	 const char *filename = get_test_file_path("explicit_matrix.npy");
	 explicit_matrix_t explicit_matrix = 
	 new_random_symmetric_explicit_matrix(10);
	 FILE *outfile = fopen(filename,"w");
	 save_numpy_explicit_matrix(outfile,explicit_matrix);
	 fclose(outfile);
	 FILE *infile = fopen(filename,"r");
	 explicit_matrix_t loaded_explicit_matrix = 
	 new_explicit_matrix_from_numpy(infile);
	 fclose(infile);
	 assert_that(compare_explicit_matrices(explicit_matrix,
					       loaded_explicit_matrix) == EQUAL);
	 free_explicit_matrix(explicit_matrix);
	 free_explicit_matrix(loaded_explicit_matrix);
	);
