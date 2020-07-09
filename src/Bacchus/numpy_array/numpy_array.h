#ifndef __NUMPY_ARRAY__
#define __NUMPY_ARRAY__

#include <stdlib.h>
#include <stdio.h>

struct _numpy_array_;
typedef struct _numpy_array_ *numpy_array_t;
typedef struct
{
	size_t num_rows;
	size_t num_columns;
} numpy_header_t;

numpy_array_t read_numpy_file(FILE *numpy_file);

numpy_header_t get_numpy_header(numpy_array_t array);

double *get_numpy_array_elements(numpy_array_t array);

void free_numpy_array(numpy_array_t array);

#endif
