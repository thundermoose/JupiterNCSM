#include <stdlib.h>
#include <stdio.h>
#include <arguments/arguments.h>
#include <index_list/index_list.h>

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments =
	       	parse_argument_list(num_arguments,argument_list);
	if (should_display_usage(arguments))
	{
		display_usage(arguments);
		free_arguments(arguments);
		return EXIT_FAILURE;
	}
	index_list_t index_list = 
		parse_human_readable_index_list(get_input_list(arguments));
	save_index_list(index_list,
			get_output_list(arguments));
	free_index_list(index_list);
	free_arguments(arguments);
}
