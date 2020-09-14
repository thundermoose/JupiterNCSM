#include <stdlib.h>
#include <stdio.h>
#include <unit_testing/test.h>
#include <log/log.h>
#include <settings/settings.h>
#include <matrix/matrix.h>
#include <combination_table/combination_table.h>
#include <execution_order/execution_order.h>
#include <lanczos/lanczos.h>
#include <eigen_system/eigen_system.h>
#include <string_tools/string_tools.h>


__attribute__((constructor(101)))
void setting_up_log()
{
	initiate_logging("BACCHUS_LOG_FILE",
			 "bacchus.log");
}

int main(int num_arguments, char **argument_list)
{
	settings_t settings = parse_settings(num_arguments,
					     argument_list);
	if (should_show_help_text(settings))
	{
		show_help_text(settings);
		return EXIT_FAILURE;
	}
	combination_table_t combination_table =
		new_combination_table
		(get_combination_table_path_setting(settings),
		 get_num_protons_setting(settings),
		 get_num_neutrons_setting(settings));	      
	execution_order_t execution_order =
		read_execution_order(get_execution_order_path_setting(settings),
				    combination_table);
	lanczos_settings_t lanczos_settings =
	{
		.dimension = get_full_dimension(combination_table),
		.vector_settings = setup_vector_settings(combination_table),
		.krylow_vectors_directory_name = 
			copy_string
			(get_krylow_vector_directory_setting(settings)),
		.max_num_iterations = 
			get_max_num_lanczos_iterations_setting(settings),
		.eigenvalue_tollerance = get_tollerance_setting(settings),
		.matrix = new_generative_matrix
			(execution_order,
			 combination_table,
			 get_index_lists_base_directory_setting(settings),
			 get_matrix_file_base_directory_setting(settings))
	};
	lanczos_environment_t lanczos_environment =
		new_lanczos_environment(lanczos_settings);
	diagonalize(lanczos_environment);
	eigen_system_t eigen_system = get_eigensystem(lanczos_environment);
	print_eigen_system(eigen_system);
	free_eigen_system(eigen_system);
	free_lanczos_environment(lanczos_environment);
	free_matrix(lanczos_settings.matrix);
	free_execution_order(execution_order);
	free_combination_table(combination_table);
	free_settings(settings);
	free(lanczos_settings.krylow_vectors_directory_name);
	return EXIT_SUCCESS;
}

__attribute__((destructor(900)))
void finalize()
{
	finalize_logging();
}

new_test(a_simple_test,
	 printf("The test works\n");
	 assert_that(1);
	);
