#include "write_to_uniform.h"
#include <directory_tools/directory_tools.h>
#include <string_tools/string_tools.h>
#include <debug_mode/debug_mode.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <errno.h>

#define BUFFER_SIZE 2048

const double acceptance_tolerance = 1e-10;

const char *generate_file_name(const Uniform_Out_File *out_file,
			       const char *name,...)
{
	va_list argument_list;
	va_start(argument_list,name);
	char name_buffer[BUFFER_SIZE/2];
	vsprintf(name_buffer,name,argument_list);
	static char buffer[BUFFER_SIZE];
	snprintf(buffer,
		 BUFFER_SIZE,
		 "%s/%s",
		 out_file->directory_name,
		 name_buffer);
	return buffer;
}

void create_header_file(Uniform_Out_File *out_file,
			const size_t num_particles,
			const quantum_number e_max,
			const size_t dimension)
{
	out_file->header_file = 
		fopen(generate_file_name(out_file, "header"),"w");
	if (out_file->header_file == NULL)
	{
		fprintf(stderr,"Could not create header file %s, %s\n",
			generate_file_name(out_file,"header"),
			strerror(errno));
	}
	fprintf(out_file->header_file,
		"=====Header=====\n"
		"num particles: %lu\n"
		"dir: %s\n"
		"Nmax: %d\n"
		"dim: %lu\n",
		num_particles,
		out_file->directory_name,
		e_max,
		dimension);
}

void create_2p_basis_file(Uniform_Out_File *out_file,
			  m_scheme_2p_basis_t ms_basis)
{
	out_file->basis_file =
		fopen(generate_file_name(out_file,"basis.bin"),"w");
	if (out_file->basis_file == NULL)
	{
		fprintf(stderr,"Could not create basis file %s, %s\n",
			generate_file_name(out_file,"basis.bin"),
			strerror(errno));
	}
	m_scheme_2p_state_t *states =
		get_m_scheme_2p_states(ms_basis);
	const size_t dimension = get_m_scheme_2p_dimension(ms_basis);
	out_file->phase = (double*)malloc(dimension*sizeof(double));
	for (size_t i = 0; i < dimension; i++)
		if (states[i].a % 2 < states[i].b % 2)
		{
			sp_state_index t = states[i].a;
			states[i].a = states[i].b;
			states[i].b = t;
			out_file->phase[i] = -1;
		}
		else
			out_file->phase[i] = 1;
	if (fwrite(states,
		   sizeof(m_scheme_2p_state_t),
		   dimension,
		   out_file->basis_file) != dimension)
	{
		fprintf(stderr,"Could not write basis file\n");
	}
	free(states);
}

double permute_3p_state(M_Scheme_3p_State *state)
{
	double phase = 1;
#define permute(a,b) if ((a)%2 < (b)%2)\
	{\
		sp_state_index t = a;\
		a = b;\
		b = t;\
		phase*=-1;\
	}
	permute(state->a,state->b);
	permute(state->b,state->c);
	permute(state->a,state->b);
	return phase;
}

void create_3p_basis_file(Uniform_Out_File *out_file,
			  M_Scheme_3p_Basis *ms_basis)
{
	out_file->basis_file =
		fopen(generate_file_name(out_file,"basis.bin"),"w");
	if (out_file->basis_file == NULL)
	{
		fprintf(stderr,"Could not create basis file %s, %s\n",
			generate_file_name(out_file,"basis.bin"),
			strerror(errno));
	}
	M_Scheme_3p_State *states = get_m_scheme_3p_states(ms_basis);
	out_file->phase = (double*)malloc(ms_basis->dimension*sizeof(double));
	for (size_t i = 0; i < ms_basis->dimension; i++)
	{
		out_file->phase[i] = permute_3p_state(&states[i]);
	}
	if (fwrite(ms_basis->states,
		   sizeof(M_Scheme_3p_State),
		   ms_basis->dimension,
		   out_file->basis_file) != ms_basis->dimension)
	{
		fprintf(stderr,"Could not write basis file\n");
	}
	free(states);
}

Uniform_Out_File *create_new_uniform_2p_out_file(const char *directory_name,
						 m_scheme_2p_basis_t ms_basis)
{
	if (!directory_exists(directory_name))
		create_directory(directory_name);
	Uniform_Out_File *out_file = 
		(Uniform_Out_File*)malloc(sizeof(Uniform_Out_File));
	out_file->directory_name = copy_string(directory_name);  
	out_file->data_blocks = NULL;
	out_file->block_status = NULL;
	out_file->current_number_blocks = 0;
	out_file->max_number_blocks = 0;
	create_header_file(out_file,
			   2, // number of particles	
			   get_m_scheme_2p_e_max1(ms_basis),
			   get_m_scheme_2p_dimension(ms_basis));
	create_2p_basis_file(out_file,
			     ms_basis);
	return out_file;
}

Uniform_Out_File *create_new_uniform_3p_out_file(const char *directory_name,
						 M_Scheme_3p_Basis *ms_basis)
{
	if (!directory_exists(directory_name))
		create_directory(directory_name);
	Uniform_Out_File *out_file = 
		(Uniform_Out_File*)malloc(sizeof(Uniform_Out_File));
	out_file->directory_name = copy_string(directory_name);  
	out_file->data_blocks = NULL;
	out_file->block_status = NULL;
	out_file->current_number_blocks = 0;
	out_file->max_number_blocks = 0;
	create_header_file(out_file,
			   3, // number of particles	
			   ms_basis->e_max,
			   ms_basis->dimension);
	create_3p_basis_file(out_file,
			     ms_basis);
	return out_file;
}

void expand_blocks_if_needed(Uniform_Out_File *out_file)
{
	if (out_file->current_number_blocks < out_file->max_number_blocks)
		return;
	out_file->max_number_blocks+=out_file->max_number_blocks+1;
	out_file->data_blocks = 
		(Data_Block*)realloc(
				     out_file->data_blocks,
				     sizeof(Data_Block)*
				     out_file->max_number_blocks);
	out_file->block_status =
		(block_status_t*)realloc(
					 out_file->block_status,
					 sizeof(block_status_t)*
					 out_file->max_number_blocks);
}

typedef struct
{
	size_t bra_index;
	size_t ket_index;
} config_t;

size_t add_uniform_block(Uniform_Out_File *out_file,
			 Data_Block new_block)
{
	size_t block_number = 0;
#pragma omp critical(uniform_block)
	{
		expand_blocks_if_needed(out_file);
		out_file->data_blocks[out_file->current_number_blocks] =
			new_block;	
		out_file->block_status[out_file->current_number_blocks] =
			unused_block;
		block_number = out_file->current_number_blocks++;
	}
	return block_number;
}

size_t setup_diagonal_block(double **elements_ptr,
			    config_t **configurations_ptr,
			    Dens_Matrix block_matrix,
			    size_t *m_indices,
			    size_t *n_indices,
			    double *phases)
{
	// We are going to set these two pointer later
	// so if they already posize_t to content, that content
	// might get lost
	assert(*elements_ptr == NULL);
	assert(*configurations_ptr == NULL);

	// Diagonal blocks must be square matricies
	assert(block_matrix.m == block_matrix.n);
	const size_t matrix_side = block_matrix.m;

	// Diagonal blocks are symmetric, so we only need to store
	// the upper triangular part and the diagonal
	const size_t max_number_elements = ((matrix_side+1)*matrix_side)/2;
	double *elements =
		(double*)malloc(max_number_elements*sizeof(double));
	config_t *configurations =
		(config_t*)malloc(max_number_elements*sizeof(config_t));

	// The final size of elements and configurations arrays may be 
	// smaller than max_number_elements since we are going to filter out
	// small values	
	size_t next_element_index = 0;

	// Loop over the upper half and diagonal of the block_matrix
	for (size_t i = 0; i<matrix_side; i++)
	{
		for (size_t j = i; j<matrix_side; j++)
		{
			config_t current_config =
			{
				// The configurations refer to
				// the indices in the full Hamiltonian
				// matrix and not just the block matrix
				// that is why I add the offsets
				.bra_index = m_indices[i],
				.ket_index = n_indices[j]
			};	
			double phase = 
				phases[current_config.bra_index]*
				phases[current_config.ket_index];
			double current_element =
				block_matrix.elements[i*matrix_side+j];
			if (fabs(current_element)<acceptance_tolerance)
				continue;
			elements[next_element_index] = current_element*phase;
			configurations[next_element_index] = 
				current_config;
			next_element_index++;
		}
	}
	*elements_ptr = elements;
	*configurations_ptr = configurations;
	return next_element_index;
}

size_t setup_off_diagonal_block(double **elements_ptr,
				config_t **configurations_ptr,
				Dens_Matrix block_matrix,
				size_t *m_indices,
				size_t *n_indices,
				double *phases)
{
	// We are going to set these two pointer later
	// so if they already posize_t to content, that content
	// might get lost
	assert(*elements_ptr == NULL);
	assert(*configurations_ptr == NULL);

	// The maximum_number of elements that we might need to store
	// are all of the block_matrix elements
	const size_t max_number_elements = block_matrix.m*block_matrix.n;
	double *elements =
		(double*)malloc(max_number_elements*sizeof(double));
	config_t *configurations =
		(config_t*)malloc(max_number_elements*sizeof(config_t));

	// The final size of elements and configurations arrays may be 
	// smaller than max_number_elements since we are going to filter out
	// small values	
	size_t next_element_index = 0;

	for (size_t i = 0; i<max_number_elements; i++)
	{
		config_t current_config =
		{
			// The configurations refer to
			// the indices in the full Hamiltonian
			// matrix and not just the block matrix
			// that is why I add the offsets
			.bra_index = m_indices[i / block_matrix.n],
			.ket_index = n_indices[i % block_matrix.n]
		};	
		double phase = 
			phases[current_config.bra_index]*
			phases[current_config.ket_index];
		double current_element =
			block_matrix.elements[i];
		if (fabs(current_element)<acceptance_tolerance)
			continue;
		elements[next_element_index] = current_element*phase;
		configurations[next_element_index] = 
			current_config;
		next_element_index++;

	}
	*elements_ptr = elements;
	*configurations_ptr = configurations;
	return next_element_index;
}

void write_config_file(Uniform_Out_File *out_file,
		       size_t block_number,
		       size_t number_configurations,
		       config_t *configurations)
{
	FILE* configuration_file = 
		fopen(generate_file_name(out_file,
					 "conf_%lu",
					 block_number+1),
		      "w");
	if (configuration_file == NULL)
	{
		fprintf(stderr,"Could not create conf_%lu, %s\n",
			block_number+1,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (fwrite(configurations,
		   sizeof(config_t),
		   number_configurations,
		   configuration_file) != number_configurations)
	{
		fprintf(stderr,"Could not write to conf_%lu, %s\n",
			block_number+1,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	fclose(configuration_file);
}

void write_element_file(Uniform_Out_File *out_file,
			size_t block_number,
			size_t number_elements,
			double *elements)
{
	FILE* element_file = 
		fopen(generate_file_name(out_file,
					 "elem_%lu",
					 block_number+1),
		      "w");
	if (element_file == NULL)
	{
		fprintf(stderr,"Could not create elem_%lu, %s\n",
			block_number+1,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	if (fwrite(elements,
		   sizeof(double),
		   number_elements,
		   element_file) != number_elements)
	{
		fprintf(stderr,"Could not write to elem_%lu, %s\n",
			block_number+1,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	fclose(element_file);
}

int is_diagonal(const Data_Block data_block)
{
	return data_block.E1 == data_block.E2;
}

void write_to_uniform_block(Uniform_Out_File *out_file,
			    size_t block_number,
			    Dens_Matrix block_matrix,
			    size_t *m_indices,
			    size_t *n_indices)
{
	if (block_matrix.m == 0 || block_matrix.n == 0)
		return;
#pragma omp critical(uniform_block)
	{
		Data_Block current_block = out_file->data_blocks[block_number];
		config_t *configurations = NULL;
		double *elements = NULL;
		size_t block_length = 0;
		if (is_diagonal(current_block))
			block_length =
				setup_diagonal_block(&elements,
						     &configurations,
						     block_matrix,
						     m_indices,
						     n_indices,
						     out_file->phase);
		else
			block_length =
				setup_off_diagonal_block(&elements,
							 &configurations,
							 block_matrix,
							 m_indices,
							 n_indices,
							 out_file->phase);
		if (block_length != 0)
		{	
			out_file->block_status[block_number] = used_block;
			write_config_file(out_file,
					  block_number,
					  block_length,
					  configurations);
			write_element_file(out_file,
					   block_number,
					   block_length,
					   elements);
		}
		if (configurations != NULL)
			free(configurations);
		if (elements != NULL)
			free(elements);
	}
}

void move_data_file(Uniform_Out_File* file,
		    const char *prefix,
		    const size_t old_block_number,
		    const size_t new_block_number)
{
	char *old_file_name = 
		copy_string(generate_file_name(file,
					       "%s%lu",
					       prefix,
					       old_block_number));
	char *new_file_name =
		copy_string(generate_file_name(file,
					       "%s%lu",
					       prefix,
					       new_block_number));

	if (rename(old_file_name,new_file_name))
	{
		fprintf(stderr,
			"Could not rename %s to %s, %s\n",
			old_file_name,
			new_file_name,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	free(old_file_name);
	free(new_file_name);
}

static void remove_unused_channels(Uniform_Out_File* file)
{

	size_t new_block_number = 0;
	for (size_t i = 0; i<file->current_number_blocks; i++)
	{
		if (file->block_status[i] == unused_block)
			continue;
		printf("Moving %lu to %lu\n",i+1,new_block_number+1);
		move_data_file(file,"conf_",i+1,new_block_number+1);
		move_data_file(file,"elem_",i+1,new_block_number+1);
		new_block_number++;
	}
}

static void save_channels_to_header(Uniform_Out_File* file)
{
	fprintf(file->header_file,
		"\tE1\tE2\tTz\tM\tconf_file\telem_file\n");
	size_t new_block_number = 0;
	for (size_t i = 0; i<file->current_number_blocks; i++)
	{
		if (file->block_status[i] == unused_block)
			continue;
		Data_Block current_data_block = file->data_blocks[i];
		fprintf(file->header_file,
			"\t%d\t%d\t%d\t%d\tconf_%lu\t\telem_%lu\n",
			current_data_block.E1,
			current_data_block.E2,
			current_data_block.Tz,
			current_data_block.M,
			new_block_number+1,
			new_block_number+1);
		new_block_number++;
	}
}

void close_uniform_out_file(Uniform_Out_File* file)
{
	remove_unused_channels(file);
	save_channels_to_header(file);
	fclose(file->header_file);
	fclose(file->basis_file);
	free(file->directory_name);
	if (file->data_blocks)
		free(file->data_blocks);
	if (file->block_status)
		free(file->block_status);
	free(file);
}

