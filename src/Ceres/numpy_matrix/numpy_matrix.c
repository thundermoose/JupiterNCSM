#include <numpy_matrix/numpy_matrix.h>
#include <error/error.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

void save_as_numpy_matrix(const char *file_path,
			  double *elements,
			  size_t num_rows,
			  size_t num_columns)
{
	const size_t num_elements = num_rows*num_columns;
	FILE *output = fopen(file_path,"w");
	if (output == 0)
		error("Could not open file \"%s\" for writing. %s\n",
		      file_path,strerror(errno));
	char header[512];
	size_t header_length = 
		sprintf(header,
			"\x93NUMPY\x01%cv%c{'descr': '<f8', "
			"'fortran_order': False, "
			"'shape': (%lu, %lu), }",
			0,0,
			num_rows,
			num_columns);
	for (size_t pad = 0; pad < 64 - header_length % 8; pad++)
		header[header_length+pad]=' ';
	header_length += 64 - header_length % 8;
	header[header_length-1] = 0x0a;
	if (fwrite(header,sizeof(char),header_length,output) < header_length)
		error("Could not write header to %s\n",file_path);
	if (fwrite(elements,sizeof(double),num_elements,output) < num_elements)
		error("Could not write elements to %s\n",file_path);
	fclose(output);
}
