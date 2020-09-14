#include <numpy_array/numpy_array.h>
#include <log/log.h>
#include <error/error.h>
#include <stdint.h>
#include <string.h>
#include <unit_testing/test.h>

struct _numpy_array_
{
	numpy_header_t header;
	double *data;
};

static
numpy_header_t read_numpy_header(FILE *numpy_file);

static
double *read_numpy_elements(FILE *numpy_file);

static 
char *find_parenthesis(char *string);

numpy_array_t read_numpy_file(FILE *numpy_file)
{
	numpy_array_t array =
		(numpy_array_t)calloc(1,sizeof(struct _numpy_array_));
	array->header = read_numpy_header(numpy_file);
	array->data = read_numpy_elements(numpy_file);
	return array;
}

numpy_header_t get_numpy_header(numpy_array_t array)
{
	return array->header;
}

double *get_numpy_array_elements(numpy_array_t array)
{
	const size_t num_elements =
		array->header.num_rows*array->header.num_columns;
	double *elements = (double*)malloc(num_elements*sizeof(double));
	memcpy(elements,array->data,
	       num_elements*sizeof(double));
	return elements;
}

void free_numpy_array(numpy_array_t array)
{
	free(array->data);
	free(array);
}

static
numpy_header_t read_numpy_header(FILE *numpy_file)
{
	uint64_t watermark = 0;
	if (fread(&watermark,sizeof(uint64_t),1,numpy_file) != 1)
		error("Could not read numpy watermark\n");
	log_entry("Numpy watermark is %lX (%s)\n",
		  watermark,(char*)(&watermark));
	if (watermark != 0x000159504d554e93)
		error("Numpy watermark is missing\n");
	unsigned short kind = 0;
	if (fread(&kind,sizeof(unsigned short),1,numpy_file) != 1)
		error("Could not read numpy kind\n");
	char *header_dictionary_string = NULL;
	size_t header_dictionary_string_length = 0;
	if (getline(&header_dictionary_string,
		    &header_dictionary_string_length,
		    numpy_file)<0)
		error("Could not read numpy dictionary string\n");
	numpy_header_t header;
	int scan_result = 0;
	char *tuple = find_parenthesis(header_dictionary_string);
	if ((scan_result = sscanf(tuple,
	      "(%lu, %lu)",
	      &header.num_rows,&header.num_columns)) != 2)
		error("Could not parse header tupple(%d): \"%s\"\n",
		      scan_result,
		      header_dictionary_string);
	free(header_dictionary_string);
	return header;
}

static
double *read_numpy_elements(FILE *numpy_file)
{
	long current_position = ftell(numpy_file);
	fseek(numpy_file,0,SEEK_END);
	long end_position = ftell(numpy_file);
	fseek(numpy_file,current_position,SEEK_SET);
	size_t num_bytes = end_position-current_position;
	if (num_bytes % 8 != 0)
		error("The array part of the numpy file"
		      " do not contain doubles\n");
	size_t length = num_bytes/8;
	double *array = (double*)malloc(num_bytes);
	if (fread(array,sizeof(double),length,numpy_file) != length)
		error("Could not read ther array part of the numpy file\n");
	return array;
}

static
char *find_parenthesis(char *string)
{
	while (*string != '(' && *string != 0) string++;
	return string;
}

new_test(read_random_10x10_matrix,
	 FILE *file = fopen(TEST_DATA "random_matrix.npy","r");
	numpy_array_t array = read_numpy_file(file);
	assert_that(array->header.num_rows == 10);
	assert_that(array->header.num_columns == 10);
	printf("Matrix:\n");
	const size_t num_elements =
	array->header.num_rows*array->header.num_columns;
	for (size_t i = 0; i<num_elements; i++)
	{
		printf("%5.2g ",array->data[i]);
		if ((i+1) % array->header.num_columns == 0)
			printf("\n");
	}
	free_numpy_array(array);
	fclose(file);
	);

