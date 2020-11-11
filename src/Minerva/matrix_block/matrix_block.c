#include <matrix_block/matrix_block.h>
#include <string_tools/string_tools.h>
#include <unused/unused.h>
#include <error/error.h>
#include <string.h>
#include <errno.h>
#include <debug_mode/debug_mode.h>

struct _matrix_block_
{
	double *matrix_elements;
	size_t num_elements;
	size_t neutron_matrix_dimension;
	size_t proton_matrix_dimension;
	char *base_directory;
	size_t block_id;
};

static
FILE *open_matrix_block_file(matrix_block_t matrix_block,const char *file_mode);

matrix_block_t new_matrix_block(size_t block_id,
				const char *base_directory)
{
	matrix_block_t matrix_block =
	       	(matrix_block_t)malloc(sizeof(struct _matrix_block_));
	matrix_block->block_id = block_id;
	matrix_block->base_directory = copy_string(base_directory);
	FILE *matrix_block_file = open_matrix_block_file(matrix_block,"r");
	if (fread(&matrix_block->neutron_matrix_dimension,
		  sizeof(size_t),1,
		  matrix_block_file) != 1)
		error("Could not read neutron_matrix_dimension "
		      "from matrix file %lu\n",
		      block_id);
	if (fread(&matrix_block->proton_matrix_dimension,
		  sizeof(size_t),1,
		  matrix_block_file) != 1)
		error("Could not read proton_matrix_dimension "
		      "from matrix file %lu\n",
		      block_id);
	if (matrix_block->neutron_matrix_dimension == 0)
		matrix_block->num_elements =
		       	matrix_block->proton_matrix_dimension;
	else if (matrix_block->proton_matrix_dimension == 0)
		matrix_block->num_elements =
			matrix_block->neutron_matrix_dimension;
	else
		matrix_block->num_elements =
			matrix_block->neutron_matrix_dimension*
			matrix_block->proton_matrix_dimension;
	matrix_block->matrix_elements =
		(double*)malloc(matrix_block->num_elements*sizeof(double));
	if (fread(matrix_block->matrix_elements,
		  sizeof(double),
		  matrix_block->num_elements,
		  matrix_block_file) < matrix_block->num_elements)
		error("Could not read matrix elements form matrix file %lu\n",
		      block_id);
	fclose(matrix_block_file);
	return matrix_block;
}

double *get_matrix_block_elements(const matrix_block_t matrix_block)
{
	return matrix_block->matrix_elements;
}

size_t get_neutron_matrix_dimension(const matrix_block_t matrix_block)
{
	return matrix_block->neutron_matrix_dimension;
}

size_t get_proton_matrix_dimension(const matrix_block_t matrix_block)
{
	return matrix_block->proton_matrix_dimension;
}

void free_matrix_block(matrix_block_t matrix_block)
{
	free(matrix_block->matrix_elements);
	free(matrix_block->base_directory);
	free(matrix_block);
}

	static
FILE *open_matrix_block_file(matrix_block_t matrix_block,const char *file_mode)
{
	char filename[2048];
	sprintf(filename,
		"%s/%lu_matrix_elements",
		matrix_block->base_directory,
		matrix_block->block_id);
	FILE *block_file = fopen(filename,file_mode);
	if (block_file == NULL)
		error("Could not open block file %s. %s\n",
		      filename,
		      strerror(errno));
	return block_file;
}
