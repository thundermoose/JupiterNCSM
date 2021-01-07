#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <settings/settings.h>
#include <norm_matrix/norm_matrix.h>
#include <subspace_operator/subspace_operator.h>
#include <string_tools/string_tools.h>
#include <vector/vector.h>

int main(int num_arguments, char **argument_list)
{
	settings_t settings = parse_settings(num_arguments, argument_list);
	if (settings_should_show_help_text(settings))
	{
		settings_show_help_text(settings);
		free_settings(settings);
		return EXIT_FAILURE;
	}
	combination_table_t combination_table =
		new_combination_table(settings_combination_table_path(settings),
				      settings_num_protons(settings),
				      settings_num_neutrons(settings));
	size_t num_training_vectors = settings_num_training_vectors(settings);
	vector_t *training_vectors =
	       	(vector_t*)malloc(num_training_vectors*sizeof(vector_t));
	vector_settings_t vector_settings =
	       	setup_vector_settings(combination_table);
	for (size_t i = 0; i<num_training_vectors; i++)
	{
		const char *training_vector_path =
			settings_training_vector_path(settings,i);
		vector_settings_t current_vector_settings = vector_settings;
		current_vector_settings.directory_name = 
			copy_string(training_vector_path);
		training_vectors[i] =
		       	new_existing_vector(current_vector_settings);
		free(current_vector_settings.directory_name);
	}
	const char *norm_matrix_path = settings_norm_matrix_path(settings);
	create_norm_matrix(norm_matrix_path,
			   training_vectors,
			   num_training_vectors);
	const char *workspace_path = settings_workspace_path(settings);
	const char *indexlist_path = settings_index_list_path(settings);
	evaluation_order_t evaluation_order = 
		read_evaluation_order(settings_evaluation_order_path(settings),
				     combination_table);
	size_t max_loaded_memory = settings_max_loaded_memory(settings);
	for (size_t i = 0; i<settings_num_operators(settings); i++)
	{
		const char *operator_path = 
			settings_operator_path(settings,i);
		const char *subspace_operator_path =
			settings_subspace_operator_path(settings,i);
		create_subspace_operator(subspace_operator_path,
					 operator_path,
					 workspace_path,
					 indexlist_path,
					 combination_table,
					 evaluation_order,
					 vector_settings,
					 training_vectors,
					 num_training_vectors,
					 max_loaded_memory);		
	}
	for (size_t i = 0; i<num_training_vectors; i++)
		free_vector(training_vectors[i]);
	free_settings(settings);
	return EXIT_SUCCESS;
}
