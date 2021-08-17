#include <arguments/arguments.h>
#include <string_tools/string_tools.h>
#include <array_builder/array_builder.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>
#include <error/error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _arguments_
{
	int to_few_arguments;
	int single_block;
	int no_2nf;
	int lec_set;
	int exclude_kinetic_energy;
	char *program_name;
	char *interaction_path_2nf;
	char *interaction_path_3nf;
	char *combination_file_path;
	char *index_list_path;
	char *output_path;
	char *finished_energy_blocks;
	size_t block_id;
	size_t num_protons;
	size_t num_neutrons;
	size_t max_loaded_memory;
	int single_particle_energy;
	int two_particle_energy;
	double lec_CE;
	double lec_CD;
	double lec_C1;
	double lec_C3;
	double lec_C4;
};

#define INTEGER_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	if (!is_integer(argument_list[++i]))\
	error(name " takes integer values only\n");\
	arguments->field = atoi(argument_list[i]);\
	continue;\
}

#define INTEGER64_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	if (!is_integer(argument_list[++i]))\
	error(name " takes integer values only\n");\
	arguments->field = atoll(argument_list[i]);\
	continue;\
}

#define MEMORY_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	if (!is_memory_string(argument_list[++i]))\
	error(name " takes memory values only\n");\
	arguments->field = parse_memory_string(argument_list[i]);\
	continue;\
}

#define MODE_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{ \
	arguments->field = 1; \
	continue;\
}

#define STRING_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	arguments->field = argument_list[++i];\
	continue;\
}

#define LEC_ARGUMENT(name,field) \
	if (strcmp(argument_list[i],name) == 0) \
{\
	if (!is_double(argument_list[++i]))\
	error(name " takes floating point argument\n");\
	arguments->field = atof(argument_list[i]);\
	arguments->lec_set = 1;\
	continue;\
}

char default_finished_block_file[] = "finished_blocks";
arguments_t parse_argument_list(int num_arguments,
				char **argument_list)
{
	arguments_t arguments =
		(arguments_t)calloc(1,sizeof(struct _arguments_));
	arguments->program_name = *argument_list;
	if (num_arguments >= 4)
	{
		arguments->single_block = 0;
		arguments->no_2nf = 0;
		arguments->exclude_kinetic_energy = 0;
		arguments->combination_file_path = argument_list[1];	
		arguments->index_list_path = argument_list[2];
		arguments->output_path = argument_list[3];
		arguments->num_protons = 0; 
		arguments->num_neutrons = 0; 
		arguments->single_particle_energy = 0;
		arguments->two_particle_energy = 0;
		arguments->to_few_arguments = 1;
		arguments->max_loaded_memory = (size_t)(1)<<32;
		arguments->finished_energy_blocks = default_finished_block_file;
		arguments->lec_C1 = 0.0;
		arguments->lec_C3 = 0.0;
		arguments->lec_C4 = 0.0;
		arguments->lec_CE = 0.0;
		arguments->lec_CD = 0.0;
		for (size_t i = 4; i<num_arguments; i++)
		{
			log_entry("argument_list[%lu] = %s",i,
				  argument_list[i]);
			INTEGER64_ARGUMENT("--num-protons",num_protons);
			INTEGER64_ARGUMENT("--num-neutrons",num_neutrons);
			INTEGER_ARGUMENT("--single-particle-energy",
					 single_particle_energy);
			INTEGER_ARGUMENT("--two-particle-energy",
					 two_particle_energy);
			INTEGER_ARGUMENT("--block-id",block_id);
			MEMORY_ARGUMENT("--max-loaded-memory",
					 max_loaded_memory);
			MODE_ARGUMENT("--single-block",single_block);
			MODE_ARGUMENT("--no-2nf",no_2nf);
			MODE_ARGUMENT("--exclude-kinetic-energy",
				      exclude_kinetic_energy);
			STRING_ARGUMENT("--finished-blocks-file",
					finished_energy_blocks);
			LEC_ARGUMENT("--LEC-CE",lec_CE);
			LEC_ARGUMENT("--LEC-CD",lec_CD);
			LEC_ARGUMENT("--LEC-C1",lec_C1);
			LEC_ARGUMENT("--LEC-C3",lec_C3);
			LEC_ARGUMENT("--LEC-C4",lec_C4);
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
		if (arguments->no_2nf)
			arguments->interaction_path_3nf = 
				arguments->interaction_path_2nf;
		log_entry("arguments->num_protons = %lu",
			  arguments->num_protons);
		log_entry("arguments->num_neutrons = %lu",
			  arguments->num_neutrons);
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
	       "[--two-particle-energy <integer>] "
	       "[--max-loaded-memory <integer>] "
	       "[--single-block] "
	       "[--block-id <integer>] "
	       "[--finished-blocks-file <filepath>] "
	       "[--no-2nf] "
	       "[--LEC-CE <float>] "
	       "[--LEC-CD <float>] "
	       "[--LEC-C1 <float>] "
	       "[--LEC-C3 <float>] "
	       "[--LEC-C4 <float>] "
	       "<interaction path 2nf> "
	       "[interaction path 3nf]\n",
	       arguments->program_name);
}

int single_block_mode(const arguments_t arguments)
{
	return arguments->single_block;	
}

int no_2nf_argument(const arguments_t arguments)
{
	return arguments->no_2nf;
}

int lec_arguments_set(const arguments_t arguments)
{
	return arguments->lec_set;
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

const char *get_finished_energy_blocks_argument(const arguments_t arguments)
{
	return arguments->finished_energy_blocks;
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

size_t get_max_loaded_memory_argument(const arguments_t arguments)
{
	return arguments->max_loaded_memory;
}

double get_CE_lec_argument(const arguments_t arguments)
{
	return arguments->lec_CE;
}

double get_CD_lec_argument(const arguments_t arguments)
{
	return arguments->lec_CD;
}

double get_C1_lec_argument(const arguments_t arguments)
{
	return arguments->lec_C1;
}

double get_C3_lec_argument(const arguments_t arguments)
{
	return arguments->lec_C3;
}

double get_C4_lec_argument(const arguments_t arguments)
{
	return arguments->lec_C4;
}

int get_exclude_kinetic_energy_argument(const arguments_t arguments)
{
	return arguments->exclude_kinetic_energy;
}

void free_arguments(arguments_t arguments)
{
	free(arguments);
}
