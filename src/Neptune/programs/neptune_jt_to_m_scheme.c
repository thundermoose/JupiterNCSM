#include <stdlib.h>
#include <stdio.h>
#include <bases/m_scheme_2p_basis.h>
#include <input/read_2nf_antoine_format.h>
#include <output/out_file.h>
#include <transform_scheduller/transform_scheduller_2p.h>
#include <string.h>

typedef struct
{
	char *program_name;
	char *input_file_name;
	char *output_file_name;
	quantum_number e_max1;
	quantum_number e_max2;
	size_t num_particles;
} arguments_t;

static
arguments_t parse_arguments(size_t num_arguments,
			    char **argument_list);

static
int not_enough_input_arguments(arguments_t arguments);

static
quantum_number get_energy_max_one_particle(arguments_t arguments);

static
quantum_number get_energy_max_two_particle(arguments_t arguments);

static
size_t num_particles(arguments_t arguments);

static
const char *input_file_name(arguments_t arguments);

static
const char *output_file_name(arguments_t arguments);

static
void display_usage(arguments_t arguments);

static
int string_to_int(char *string);

static
int is_number(char *string);

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments = parse_arguments(num_arguments,argument_list);
	if (not_enough_input_arguments(arguments))
	{
		display_usage(arguments);
		return EXIT_SUCCESS;
	}
	quantum_number e_max1 = get_energy_max_one_particle(arguments);
	quantum_number e_max2 = get_energy_max_two_particle(arguments);
	antoine_2nf_file_t data_file =
		open_antoine_2nf_file(input_file_name(arguments),
				      num_particles(arguments),
				      e_max1,
				      e_max2);
	m_scheme_2p_basis_t basis = new_m_scheme_2p_basis(e_max1,e_max2);
	print_m_scheme_2p_basis(basis);
	out_file_t output_file =
		create_new_out_file(output_file_name(arguments),
				    basis,
				    unified_2p_file);
	transform_2p_data(data_file,output_file,basis);
	close_out_file(output_file);
	free_antoine_2nf_file(data_file);
	free_m_scheme_2p_basis(basis);
	return EXIT_SUCCESS;
}

	static
arguments_t parse_arguments(size_t num_arguments,
			    char **argument_list)
{
	arguments_t arguments =
	{
		.program_name = argument_list[0],
		.input_file_name = NULL,
		.output_file_name = NULL,
		.num_particles = 2,
		.e_max1 = 0,
		.e_max2 = 0
	};
	if (num_arguments>=3)
	{
		arguments.input_file_name = argument_list[1];
		arguments.output_file_name = argument_list[2];
		for (size_t i = 3; i<num_arguments; i++)
		{
			if (strcmp(argument_list[i],"--energy-max") == 0)
			{
				quantum_number e_max =
					string_to_int(argument_list[++i]);
				arguments.e_max1 = e_max;
				arguments.e_max2 = e_max;
			}
			else if (strcmp(argument_list[i],
					"--energy-max-1") == 0)
			{
				arguments.e_max1 =
					string_to_int(argument_list[++i]);
			}
			else if (strcmp(argument_list[i],
					"--energy-max-2") == 0)
			{
				arguments.e_max2 =
					string_to_int(argument_list[++i]);
			}
			else if (strcmp(argument_list[i],
					"--num-particles") == 0)
			{
				arguments.num_particles =
					string_to_int(argument_list[++i]);
			}
			else
			{
				fprintf(stderr,"Unknown argument \"%s\"\n",
					argument_list[i]);
				exit(EXIT_FAILURE);
			}
		}
	}
	return arguments;
}

	static
int not_enough_input_arguments(arguments_t arguments)
{
	return arguments.input_file_name == NULL;
}

	static
quantum_number get_energy_max_one_particle(arguments_t arguments)
{
	return arguments.e_max1;
}

	static
quantum_number get_energy_max_two_particle(arguments_t arguments)
{
	return arguments.e_max2;
}

	static
size_t num_particles(arguments_t arguments)
{
	return arguments.num_particles;
}

	static
const char *input_file_name(arguments_t arguments)
{
	return arguments.input_file_name;
}

	static
const char *output_file_name(arguments_t arguments)
{
	return arguments.output_file_name;
}

	static
const char *program_name(arguments_t arguments)
{
	return arguments.program_name;
}

	static
void display_usage(arguments_t arguments)
{
	printf("Usage: %s <input_file> <output_file> [optional arguments]\n",
	       program_name(arguments));
	printf("The input_file should be the path to"
	       " a antoine style two-nucleon force file.\n"
	       "The output_file should be the name of the"
	       "generated interaction file (which technically is a"
	       "directory)\n"
	       "There are four optional arguments:\n"
	       "--num-particles <int>: The number of particles of the nucleus,"
	       " default 2\n"
	       "--energy-max-1 <int>: The single particle max harmonic"
	       " oscillator energy, default 0\n"
	       "--energy-max-2 <int>: The two particle max harmonic oscillator"
	       " energy, default 0\n"
	       "--energy-max <int>: Sets the single particle and two particle"
	       "max harmonic oscillator energy to the same value\n");
}

	static
int string_to_int(char *string)
{
	if (!is_number(string))
	{
		fprintf(stderr,"\"%s\" is not a number\n",string);
		exit(EXIT_FAILURE);
	}
	return atoi(string);
}

	static
int is_number(char *string)
{
	while (*string != 0)
	{
		if (*string < '0' || *string > '9')
			return 0;
		string++;
	}
	return 1;
}
