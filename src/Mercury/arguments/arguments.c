#include <arguments/arguments.h>
#include <string_tools/string_tools.h>
#include <array_builder/array_builder.h>
#include <error/error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _arguments_
{
	int to_few_arguments;
	int single_block;
	char *program_name;
	char *interaction_path_2nf;
	char *interaction_path_3nf;
	char *combination_file_path;
	char *index_list_path;
	char *output_path;
	size_t block_id;
	size_t num_protons;
	size_t num_neutrons;
	int single_particle_energy;
	int two_particle_energy;
};

#define INTEGER_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	if (!is_integer(argument_list[++i]))\
	error(name " takes integer values only\n");\
	arguments->field = atoi(argument_list[i]);\
	continue;\
}

#define MODE_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{ \
	arguments->field = 1; \
	continue;\
}

arguments_t parse_argument_list(int num_arguments,
				char **argument_list)
{
	arguments_t arguments =
		(arguments_t)calloc(1,sizeof(struct _arguments_));
	arguments->program_name = *argument_list;
	if (num_arguments >= 4)
	{
		arguments->single_block = 0;
		arguments->combination_file_path = argument_list[1];	
		arguments->index_list_path = argument_list[2];
		arguments->output_path = argument_list[3];
		arguments->num_protons = 0; 
		arguments->num_neutrons = 0; 
		arguments->single_particle_energy = 0;
		arguments->two_particle_energy = 0;
		arguments->to_few_arguments = 1;
		for (size_t i = 4; i<num_arguments; i++)
		{
			INTEGER_ARGUMENT("--num-protons",num_protons);
			INTEGER_ARGUMENT("--num-neutrons",num_neutrons);
			INTEGER_ARGUMENT("--single-particle-energy",
					 single_particle_energy);
			INTEGER_ARGUMENT("--two-particle-energy",
					 two_particle_energy);
			INTEGER_ARGUMENT("--block-id",block_id);
			MODE_ARGUMENT("--single-block",single_block);
			if (arguments->interaction_path_2nf == NULL)
			{
				arguments->interaction_path_2nf = 
					argument_list[i];
				arguments->to_few_arguments = 0;
			}
			else if (arguments->interaction_path_3nf == NULL)
			{
				arguments->interaction_path_3nf =
					argument_list[i];
			}
			else
			{
				arguments->to_few_arguments = 1;
			}
				
		}
	}
	else
	{
		arguments->to_few_arguments = 1;
	}
	return arguments;
}

int to_few_arguments(const arguments_t arguments)
{
	return arguments->to_few_arguments;
}

void show_usage(const arguments_t arguments)
{
	printf("Usage: %s <combination file path> "
	       "<index list base directory> <output path> "
	       "[--num-protons <integer>] [--num-neutrons <integer>] "
	       "[--single-particle-energy <integer>] "
	       "[--two-particle-energy <integer>] [--single-block] "
	       "[--block-id <integer>] <interaction path 2nf> "
	       "[interaction path 3nf]\n",
	       arguments->program_name);
}

int single_block_mode(const arguments_t arguments)
{
	return arguments->single_block;	
}

const char *get_interaction_path_2nf_argument(const arguments_t arguments)
{
	return arguments->interaction_path_2nf;
}

const char *get_interaction_path_3nf_argument(const arguments_t arguments)
{
	return arguments->interaction_path_3nf;
}

const char *get_combination_file_path_argument(const arguments_t arguments)
{
	return arguments->combination_file_path;
}

const char *get_index_list_path_argument(const arguments_t arguments)
{
	return arguments->index_list_path;
}

const char *get_output_path_argument(const arguments_t arguments)
{
	return arguments->output_path;
}

size_t get_block_id(const arguments_t arguments)
{
	return arguments->block_id;
}

size_t get_num_protons_argument(const arguments_t arguments)
{
	return arguments->num_protons;
}

size_t get_num_neutrons_argument(const arguments_t arguments)
{
	return arguments->num_neutrons;
}

size_t get_num_particles_argument(const arguments_t arguments)
{
	return arguments->num_neutrons + arguments->num_protons;
}

int get_single_particle_energy_argument(const arguments_t arguments)
{
	return arguments->single_particle_energy;
}

int get_two_particle_energy_argument(const arguments_t arguments)
{
	return arguments->two_particle_energy;
}

void free_arguments(arguments_t arguments)
{
	free(arguments);
}
