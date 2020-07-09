#include <arguments/arguments.h>
#include <string_tools/string_tools.h>
#include <error/error.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _arguments_
{
	int to_few_arguments;
	int single_block;
	char *program_name;
	char *interaction_path;
	char *combination_file_path;
	char *index_list_path;
	char *output_path;
	size_t block_id;
	size_t num_protons;
	size_t num_neutrons;
	int energy_max;
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
}

arguments_t parse_argument_list(int num_arguments,
				char **argument_list)
{
	arguments_t arguments =
		(arguments_t)calloc(1,sizeof(struct _arguments_));
	arguments->program_name = *argument_list;
	if (num_arguments >= 5)
	{
		arguments->single_block = 0;
		arguments->interaction_path = argument_list[1];	
		arguments->combination_file_path = argument_list[2];	
		arguments->index_list_path = argument_list[3];
		arguments->output_path = argument_list[4];
		arguments->num_protons = 0; 
		arguments->num_neutrons = 0; 
		arguments->energy_max = 0;
		for (size_t i = 5; i<num_arguments; i++)
		{
			INTEGER_ARGUMENT("--num-protons",num_protons);
			INTEGER_ARGUMENT("--num-neutrons",num_neutrons);
			INTEGER_ARGUMENT("--energy-max",energy_max);
			INTEGER_ARGUMENT("--block-id",block_id);
			MODE_ARGUMENT("--single-block",single_block);
		}
		arguments->to_few_arguments = 0;
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
	printf("Usage: %s <interaction path> <combination file path> "
	       "<index list base directory> <output path> "
	       "[--num-protons <integer>] [--num-neutrons <integer>] "
	       "[--energy-max <integer>] [--single-block] "
	       "[--block-id <integer>]\n",
	       arguments->program_name);
}

int single_block_mode(const arguments_t arguments)
{
	return arguments->single_block;	
}

const char *get_interaction_path(const arguments_t arguments)
{
	return arguments->interaction_path;
}

const char *get_combination_file_path(const arguments_t arguments)
{
	return arguments->combination_file_path;
}

const char *get_index_list_path(const arguments_t arguments)
{
	return arguments->index_list_path;
}

const char *get_output_path(const arguments_t arguments)
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

int get_energy_max(const arguments_t arguments)
{
	return arguments->energy_max;
}

void free_arguments(arguments_t arguments)
{
	free(arguments);
}
