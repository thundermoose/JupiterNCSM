#include <interaction/interaction.h>
#include <interaction/basis/basis.h>
#include <interaction/energy_block/energy_block.h>
#include <log/log.h>
#include <energy_block_info/energy_block_info.h>
#include <debug_mode/debug_mode.h>
#include <error/error.h>
#include <string_tools/string_tools.h>
#include <array_builder/array_builder.h>
#include <global_constants/global_constants.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

struct _interaction_
{
	char *interaction_path;	
	header_t header;
	basis_t basis;
	energy_block_t current_block;
	energy_block_info_t current_block_info;
};

interaction_t new_interaction(const char *interaction_path)
{
	interaction_t interaction =
		(interaction_t)malloc(sizeof(struct _interaction_));
	interaction->interaction_path = copy_string(interaction_path);
	interaction->header = read_header(interaction_path);
	interaction->basis = read_basis(interaction_path,
					get_num_particles(interaction->header));
	interaction->current_block = NULL;
	energy_block_info_t empty_block = {0};
	interaction->current_block_info = empty_block;
	return interaction;
}

static
void sort_state(int *state);

const header_t get_header(const interaction_t interaction)
{
	return interaction->header;
}

double get_matrix_element(interaction_t interaction,
			  int *bra_state,
			  int *ket_state,
			  size_t num_particles,
			  int E1,int E2, int Tz, int M)
{
	if (num_particles != get_num_particles(interaction->header))
		return 0.0;
	switch (num_particles)
	{
		case 1:
			log_entry("bra_state = %d",
				  bra_state[0]);
			log_entry("ket_state = %d",
				  ket_state[0]);
			break;
		case 2:
			log_entry("bra_state = %d %d",
				  bra_state[0],
				  bra_state[1]);
			log_entry("ket_state = %d %d",
				  ket_state[0],
				  ket_state[1]);
			break;
		case 3:
			log_entry("bra_state = %d %d %d",
				  bra_state[0],
				  bra_state[1],
				  bra_state[2]);
			log_entry("ket_state = %d %d %d",
				  ket_state[0],
				  ket_state[1],
				  ket_state[2]);
			sort_state(bra_state);
			sort_state(ket_state);
			break;
		default:
			log_entry("To many particles to print states");
	}
	size_t bra_index = find_basis_state(interaction->basis,
				      bra_state,
				      num_particles);
	if (bra_index == no_index)
		return 0.0;
	size_t ket_index = find_basis_state(interaction->basis,
				      ket_state,
				      num_particles);
	if (ket_index == no_index)
		return 0.0;
	size_t block_index = find_block_index(interaction->header,
					      E1,E2,Tz,M);
	log_entry("%lu = find_block_index(interaction, %d, %d, %d, %d);",
		  block_index,E1,E2,Tz,M);
	if (block_index == no_index)
		return 0.0;
	energy_block_info_t needed_block =
	       	get_energy_block_info(interaction->header,
				      block_index);
	energy_block_t energy_block = interaction->current_block;
	if (energy_block == NULL ||
	    compare_energy_block_info(&needed_block,
				      &interaction->current_block_info) != 0)
	{
		if (energy_block != NULL)
			free_energy_block(energy_block);
		interaction->current_block = energy_block =
			open_energy_block(interaction->interaction_path,
					  get_energy_block_info(interaction->header,
								block_index));
		interaction->current_block_info = needed_block;
	}
	double element =  
		get_energy_block_element(energy_block,bra_index,ket_index);
	log_entry("%lg = get_energy_block_element(%lu, %lu, %lu);",
		  element,block_index+1,bra_index,ket_index);
	return element;
}


void free_interaction(interaction_t interaction)
{
	free(interaction->interaction_path);
	free_header(interaction->header);
	free_basis(interaction->basis);
	if (interaction->current_block != NULL)
		free_energy_block(interaction->current_block);
	free(interaction);
}

static
void sort_state(int *state)
{
	int tmp = 0;
#define swap(a,b) \
	{\
		tmp = a;\
		a = b;\
		b = tmp;\
	}
	if (state[0] > state[1])
		swap(state[0],state[1]);
	if (state[1] > state[2])
		swap(state[1],state[2]);
	if (state[0] > state[1])
		swap(state[0],state[1]);
}
