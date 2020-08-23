#include <stdlib.h>
#include <stdio.h>
#include <arguments/arguments.h>
#include <combination_table/combination_table.h>
#include <log/log.h>

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MERCURY_LOGFILE",
			 "mercury_J_scheme_to_internal.log");	 
}

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments =
	       	parse_argument_list(num_arguments, argument_list);	    
	if (to_few_arguments(arguments))
	{
		show_usage(arguments);
		return EXIT_SUCCESS;
	}
	combination_table_t combination_table =
		new_combination_table(get_combination_file_path(arguments),
				      get_num_protons_argument(arguments),
				      get_num_neutrons_argument(arguments));
	while (has_next_matrix_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_matrix_block(combination_table);
		printf("%s p: %d %d %d n: %d %d %d,"
		       "#PC: %lu #NC: %lu id: %lu\n",
		       block_type_to_string(current_matrix_block.type),
		       current_matrix_block.difference_energy_protons,
		       current_matrix_block.difference_M_protons,
		       current_matrix_block.depth_protons,
		       current_matrix_block.difference_energy_neutrons,
		       current_matrix_block.difference_M_neutrons,
		       current_matrix_block.depth_neutrons,
		       current_matrix_block.num_proton_combinations,
		       current_matrix_block.num_neutron_combinations,
		       current_matrix_block.matrix_block_id);
	}
	free_combination_table(combination_table);
	free_arguments(arguments);
}
