#include <arguments/arguments.h>
#include <array_builder/array_builder.h>
#include <stdio.h>
#include <string.h>

struct _arguments_
{
	int usage;
	char *program;
	char *combination_table_path;
	size_t num_protons;
	size_t num_neutrons;
	size_t num_operators;
	double *coefficients;
	char **operator_paths;
	char *output_path;	
};

double one = 1.0;

arguments_t parse_arguments(size_t num_arguments,
			    char **argument_list)
{
	arguments_t arguments =
	       	(arguments_t)calloc(1,sizeof(struct _arguments_));
	arguments->program = *argument_list;
	if (num_arguments < 5)
	{
		arguments->usage = 1;
		return arguments;
	}
	else
	{
		arguments->usage = 0;
		arguments->output_path = argument_list[1];
		arguments->combination_table_path = argument_list[2];
		arguments->num_protons = atoll(argument_list[3]);
		arguments->num_neutrons = atoll(argument_list[4]);
		array_builder_t operator_paths_builder = 
			new_array_builder((void**)&arguments->operator_paths,
					  &arguments->num_operators,
					  sizeof(char*));
		size_t num_coefficients = 0;
		array_builder_t coefficient_builder =
			new_array_builder((void**)&arguments->coefficients,
					  &num_coefficients,
					  sizeof(double));
		for (size_t i = 5; i<num_arguments; i++)
		{
			if (strcmp(argument_list[i],
				   "-c") == 0 ||
			    strcmp(argument_list[i],
				   "--coefficient") == 0)
			{
				arguments->coefficients[num_coefficients-1] =
					atof(argument_list[++i]);
			}
			else
			{
				append_array_element(operator_paths_builder,
						     argument_list[i]);
				append_array_element(coefficient_builder,
						     &one);
			}
		}
		free_array_builder(operator_paths_builder);
		free_array_builder(coefficient_builder);
		return arguments;
	}
}

int should_show_usage(arguments_t arguments)
{
	return arguments->usage;
}

void show_usage(arguments_t arguments)
{
	printf("Usage: %s <output_path> <comb.txt> "
	       "<num protons> <num neutrons> "
	       "[operator_path [-c/--coefficient <value>]]...\n",
	       arguments->program);
}

const char *get_combination_table_argument(arguments_t arguments)
{
	return arguments->combination_table_path;
}

size_t get_num_protons_argument(arguments_t arguments)
{
	return arguments->num_protons;
}

size_t get_num_neutron_argument(arguments_t arguments)
{
	return arguments->num_neutrons;
}

size_t num_operator_arguments(arguments_t arguments)
{
	return arguments->num_operators;
}

const double *get_coefficient_arguments(arguments_t arguments)
{
	return arguments->coefficients;
}

const char **get_operator_path_arguments(arguments_t arguments)
{
	return (const char**)arguments->operator_paths;
}

const char *get_output_path_argument(arguments_t arguments)
{
	return arguments->output_path;
}

void free_arguments(arguments_t arguments)
{
	free(arguments->coefficients);
	free(arguments->operator_paths);
	free(arguments);
}
