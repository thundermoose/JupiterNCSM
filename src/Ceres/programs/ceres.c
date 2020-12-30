#include <stdlib.h>
#include <stdio.h>
#include <settings/settings.h>
#include <norm_matrix/norm_matrix.h>
#include <subspace_operator/subspace_operator.h>

int main(int num_arguments, char **argument_list)
{
	settings_t settings = parse_settings(num_arguments,
					     argument_list);
	if (should_show_help_test(settings))
	{
		show_help_text(settings);
		return EXIT_FAILURE;
	}
	create_norm_matrix(settings);
	for (size_t i = 0; i<settings_num_operators(settings); i++)
		create_subspace_operator(settings,i);		
}
