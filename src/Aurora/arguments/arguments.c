#include <arguments/arguments.h>
#include <stdio.h>

struct _arguments_
{
	char *program_name;
	char *input_list;
	char *output_list;
};

arguments_t parse_argument_list(size_t num_arguments,
				char **arguments_list)
{
	arguments_t arguments =
	       	(arguments_t)calloc(1,sizeof(struct _arguments_));
	arguments->program_name = arguments_list[0];
	if (num_arguments == 3)
	{
		arguments->input_list = arguments_list[1];
		arguments->output_list = arguments_list[2];
	}
	return arguments;
}

int should_display_usage(arguments_t arguments)
{
	return arguments->input_list == NULL || arguments->output_list == NULL;
}

void display_usage(arguments_t arguments)
{
	printf("Usage: %s <input_list> <output_list>\n",
	       arguments->program_name);
}

char *get_input_list(arguments_t arguments)
{
	return arguments->input_list;
}

char *get_output_list(arguments_t arguments)
{
	return arguments->output_list;
}

void free_arguments(arguments_t arguments)
{
	free(arguments);
}
