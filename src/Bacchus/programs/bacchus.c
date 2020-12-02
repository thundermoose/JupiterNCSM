#include <stdlib.h>
#include <stdio.h>
#include <unit_testing/test.h>
#include <log/log.h>
#include <settings/settings.h>
#include <matrix/matrix.h>
#include <combination_table/combination_table.h>
#include <execution_order/execution_order.h>
#include <lanczos/lanczos.h>
#include <eigensystem/eigensystem.h>
#include <string_tools/string_tools.h>
#include <string.h>
#include <time.h>


__attribute__((constructor(101)))
void setting_up_log()
{
	initiate_logging("BACCHUS_LOG_FILE",
			 "bacchus.log");
}

int main(int num_arguments, char **argument_list)
{
	struct timespec t_start,t_end;
	printf("Bacchus starts:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
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
		.convergence_critera = 
			get_convergece_criteria_setting(settings),
		.target_eigenvalue = 0,
		.matrix = new_generative_matrix
			(execution_order,
			 combination_table,
			 get_index_lists_base_directory_setting(settings),
			 get_matrix_file_base_directory_setting(settings),
			 get_maximum_loaded_memory_setting(settings))
	};
	lanczos_environment_t lanczos_environment =
		new_lanczos_environment(lanczos_settings);
	diagonalize(lanczos_environment);
	eigensystem_t eigensystem = get_eigensystem(lanczos_environment);
	print_eigensystem(eigensystem);
	const char *eigenvector_directory =
	       	get_eigenvector_directory_setting(settings);
	for (size_t i = 0; i<get_target_eigenvector_setting(settings); i++)
	{
		vector_settings_t vector_setting =
		       	lanczos_settings.vector_settings;
		vector_setting.directory_name =
		       	(char*)calloc(strlen(eigenvector_directory)+256,
				      sizeof(char));
		sprintf(vector_setting.directory_name,
			"%s/eigenvector_%lu",
			eigenvector_directory,
			i+1);
		vector_t eigenvector = new_zero_vector(vector_setting);
		get_eigenvector(eigenvector,
				 eigensystem,i);
		save_vector(eigenvector);
		free_vector(eigenvector);
		free(vector_setting.directory_name);
	}
	free_eigensystem(eigensystem);
	free_lanczos_environment(lanczos_environment);
	free_matrix(lanczos_settings.matrix);
	free_execution_order(execution_order);
	free_combination_table(combination_table);
	free_settings(settings);
	free(lanczos_settings.krylow_vectors_directory_name);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double bacchus_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Bacchus ends after %lg µs\n",
	       bacchus_time);
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
