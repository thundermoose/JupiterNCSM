#include <vector_block/vector_block.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>
#include <string.h>
#include <errno.h>

struct _vector_block_
{
	size_t neutron_dimension;
	size_t proton_dimension;	
	double *elements;
	int block_id;
	char *base_directory;
};

static
FILE *open_vector_block_file(vector_block_t vector_block,
			     const char *open_mode);

vector_block_t new_vector_block(const char *base_directory,
				const basis_block_t basis_block)
{
	vector_block_t vector_block =
		(vector_block_t)malloc(sizeof(struct _vector_block_));	
	vector_block->neutron_dimension = basis_block.num_neutron_states;
	vector_block->proton_dimension = basis_block.num_proton_states;
	vector_block->block_id = basis_block.block_id;
	vector_block->base_directory = copy_string(base_directory);
	const size_t num_elements =
		vector_block->neutron_dimension*vector_block->proton_dimension;
	vector_block->elements = (double*)malloc(num_elements*sizeof(double));
	load_vector_block_elements(vector_block);
	return vector_block;
}

void load_vector_block_elements(vector_block_t vector_block)
{
	FILE *file = open_vector_block_file(vector_block,"r");
	const size_t num_elements =
		vector_block->neutron_dimension*vector_block->proton_dimension;
	log_entry("num_elements = %lu",num_elements);
	if (fread(vector_block->elements,
		  sizeof(double),
		  num_elements,
		  file) != num_elements)
		error("Could not read elements of block %d in %s\n",
		      vector_block->block_id,
		      vector_block->base_directory);
	fclose(file);
}

void save_vector_block_elements(vector_block_t vector_block)
{
	FILE *file = open_vector_block_file(vector_block,"w");
	const size_t num_elements =
		vector_block->neutron_dimension*vector_block->proton_dimension;
	if (fwrite(vector_block->elements,
		   sizeof(double),
		   num_elements,
		   file) != num_elements)
		error("Could not write elements to block %d in %s\n",
		      vector_block->block_id,
		      vector_block->base_directory);
	fclose(file);
}

size_t get_neutron_dimension(const vector_block_t vector_block)
{
	return vector_block->neutron_dimension;
}

size_t get_proton_dimension(const vector_block_t vector_block)
{
	return vector_block->proton_dimension;
}

double *get_vector_block_elements(const vector_block_t vector_block)
{
	return vector_block->elements;
}

void free_vector_block(vector_block_t vector_block)
{
	free(vector_block->base_directory);
	free(vector_block->elements);
	free(vector_block);
}

static
FILE *open_vector_block_file(vector_block_t vector_block,
			     const char *open_mode)
{
	char filename[2048];
	sprintf(filename,
		"%s/vec_%d",
		vector_block->base_directory,
		vector_block->block_id);
	FILE *file = fopen(filename,open_mode);
	if (file == NULL)
		error("Could not open vector block file %s. %s\n",
		      filename,
		      strerror(errno));
	return file;
}
