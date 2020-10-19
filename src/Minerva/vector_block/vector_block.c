#include <vector_block/vector_block.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>
#include <string.h>
#include <errno.h>
#include <omp.h>

struct _vector_block_
{
	size_t neutron_dimension;
	size_t proton_dimension;	
	size_t num_instances;
	double **elements;
	int block_id;
	char *base_directory;
};

static
FILE *open_vector_block_file(vector_block_t vector_block,
			     const char *open_mode);

static
void reduce_vector(vector_block_t vector_block);

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
	vector_block->num_instances = 1;
	vector_block->elements = (double**)malloc(sizeof(double*));
	*vector_block->elements = (double*)calloc(num_elements,sizeof(double));
	load_vector_block_elements(vector_block);
	return vector_block;
}

vector_block_t new_output_vector_block(const char *base_directory,
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
	vector_block->num_instances = omp_get_num_threads();
	vector_block->elements = 
		(double**)malloc(vector_block->num_instances*sizeof(double*));
	for (size_t i = 0; i<vector_block->num_instances; i++)
		vector_block->elements[i] = (double*)calloc(num_elements,
							    sizeof(double));
	load_vector_block_elements(vector_block);
	return vector_block;
}

void load_vector_block_elements(vector_block_t vector_block)
{
	FILE *file = open_vector_block_file(vector_block,"r");
	const size_t num_elements =
		vector_block->neutron_dimension*vector_block->proton_dimension;
	log_entry("num_elements = %lu",num_elements);
	if (fread(*vector_block->elements,
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
	reduce_vector(vector_block);
	FILE *file = open_vector_block_file(vector_block,"w");
	const size_t num_elements =
		vector_block->neutron_dimension*vector_block->proton_dimension;	
	if (fwrite(*vector_block->elements,
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
	if (vector_block->num_instances == 1)
		return *vector_block->elements;
	else
		return vector_block->elements[omp_get_thread_num()];
}

void free_vector_block(vector_block_t vector_block)
{
	free(vector_block->base_directory);
	for (size_t i = 0; i<vector_block->num_instances; i++)
		free(vector_block->elements[i]);
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

static
void reduce_vector(vector_block_t vector_block)
{
	const size_t num_elements =  
		vector_block->neutron_dimension*vector_block->proton_dimension;
#pragma omp critical(reduce_vector)
	for (size_t i = 1; i < vector_block->num_instances; i++)
		for (size_t j = 0; j<num_elements; j++)
			vector_block->elements[0][j] +=
				vector_block->elements[i][j];
}
