#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

static
void add_operator_blocks(matrix_block_setting_t block,
			 char **operator_paths,
			 double *coefficients,
			 size_t num_operators,
			 char *output_path,
			 double *workspace);

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments = parse_arguments(num_arguments,
						argument_list);
	if (should_show_usage(arguments))
	{
		show_usage(arguments);
		return EXIT_FAILURE;
	}
	combination_table_t combination_table = 
		new_combination_table(get_combination_table_argument(arguments),
				      get_num_protons_argument(arguments),
				      get_num_neutron_argument(arguments));
	size_t max_block_size = 0;
	while (has_next_matrix_block(combination_table))
	{
		matrix_block_setting_t current_block =
			next_matrix_block(combination_table);
		size_t needed_block_size =
		       	current_block.num_proton_combinations*
			current_block.num_neutron_combinations;	
		max_needed_block_size = max(max_needed_block_size,		
					    needed_block_size);
	}
	double *coefficients = get_coefficient_arguments(arguments);
	char **operator_paths = get_operator_path_arguments(arguments);
	size_t num_operators = num_operator_arguments(arguments);
	char *output_path = get_output_path_argument(arguments);
	size_t workspace_size_per_thread =
	       	max_needed_block_size*(num_operators+1);
	double *total_workspace = NULL;
	size_t total_workspace_size;
	reset_matrix_block_iterator(combination_table);
#pragma omp parallel
	{
#pragma omp single
		{
			total_workspace_size = 
				workspace_size_per_thread*
				omp_get_num_threads();
			total_workspace =
			       	(double*)malloc(total_workspace_size*
						sizeof(double));
		}
		size_t thread_id = omp_get_thread_num();
		double *workspace = 
			total_workspace_size + 
			thread_id*workspace_size_per_thread;
		while (has_next_matrix_block(combination_table))
		{
			matrix_block_setting_t current_block =
				next_matrix_block(combination_table);
			add_operator_blocks(current_block,
					    operator_paths,
					    coefficients,
					    num_operators,
					    output_path,
					    workspace);
		}
	}	
	free(total_workspace);
	free_arguments(arguments);
	return EXIT_SUCCESS;
}

static
void add_operator_blocks(matrix_block_setting_t block,
			 char **operator_paths,
			 double *coefficients,
			 size_t num_operators,
			 char *output_path,
			 double *workspace)
{
	size_t block_size = block.num_proton_combinations*
		block.num_neutron_combinations;
}
