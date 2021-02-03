#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include <matrix_block_setting/matrix_block_setting.h>
#include <combination_table/combination_table.h>
#include <arguments/arguments.h>
#include <error/error.h>
#include <errno.h>
#include <log/log.h>

#define max(a,b) ((a) > (b) ? (a) : (b))

static
void add_operator_blocks(matrix_block_setting_t block,
			 const char **operator_paths,
			 const double *coefficients,
			 size_t num_operators,
			 const char *output_path,
			 double *workspace);

static
void load_operator(double *buffer,
		   const char *operator_path,
		   size_t matrix_block_id,
		   size_t block_size);

static
void save_operator(const double *buffer,
		   const char *operator_path,
		   size_t matrix_block_id,
		   matrix_block_setting_t block);

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("VULCAN_LOGFILE",
			 "vulcan.log");	 
}

	__attribute__((destructor(900)))
void finalize()
{
	finalize_logging();
}

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
		max_block_size = max(max_block_size,		
				     needed_block_size);
	}
	const double *coefficients = get_coefficient_arguments(arguments);
	const char **operator_paths = get_operator_path_arguments(arguments);
	size_t num_operators = num_operator_arguments(arguments);
	const char *output_path = get_output_path_argument(arguments);
	size_t workspace_size_per_thread = max_block_size*2;
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
			if (total_workspace == NULL)
				error("Could not allocate needed work-space memory. %s\n",
				      strerror(errno));
		}
		size_t thread_id = omp_get_thread_num();
		double *workspace = 
			total_workspace+ 
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
			 const char **operator_paths,
			 const double *coefficients,
			 size_t num_operators,
			 const char *output_path,
			 double *workspace)
{
	size_t block_size = get_matrix_block_length(block);
	printf("Adding block %lu\n",block.matrix_block_id);
	double *output_block = workspace;
	double *input_block = workspace+block_size;
	memset(output_block,0,block_size*sizeof(double));
	for (size_t i = 0; i<num_operators; i++)
	{
		load_operator(input_block,
			      operator_paths[i],
			      block.matrix_block_id,
			      block_size);
		double coefficient = coefficients[i];
		printf("coefficient[%lu] = %lg\n",i,coefficient);
		for (size_t j = 0; j < block_size; j++)
			output_block[j]+=coefficient*input_block[j];
	}
	save_operator(output_block,
		      output_path,
		      block.matrix_block_id,
		      block);
}

static
void load_operator(double *buffer,
		   const char *operator_path,
		   size_t matrix_block_id,
		   size_t block_size)
{
	char filename_buffer[1024];
	sprintf(filename_buffer,
		"%s/%lu_matrix_elements",
		operator_path,
		matrix_block_id);
	FILE *matrix_file = fopen(filename_buffer,"r");
	if (matrix_file == 0)
	{
		memset(buffer,0,sizeof(double)*block_size);
	}
	else
	{
		fseek(matrix_file,2*sizeof(size_t),SEEK_SET);
		if (fread(buffer,sizeof(double),block_size,matrix_file) != block_size)
			error("Could not read block %lu from operator %s\n",
			      matrix_block_id,operator_path);
		fclose(matrix_file);
	}
}

static
void save_operator(const double *buffer,
		   const char *operator_path,
		   size_t matrix_block_id,
		   matrix_block_setting_t block)
{
	char filename_buffer[1024];
	sprintf(filename_buffer,
		"%s/%lu_matrix_elements",
		operator_path,
		matrix_block_id);
	FILE *matrix_file = fopen(filename_buffer,"w");
	if (fwrite(&block.num_neutron_combinations,
		   sizeof(size_t),
		   1,
		   matrix_file) != 1)
		error("Could not write neutron dimension to matrix file %lu "
		      "in operator %s\n",
		      matrix_block_id,
		      operator_path);
	if (fwrite(&block.num_proton_combinations,
		   sizeof(size_t),
		   1,
		   matrix_file) != 1)
		error("Could not write proton dimension to matrix file %lu "
		      "in operator %s\n",
		      matrix_block_id,
		      operator_path);
	size_t block_size = get_matrix_block_length(block);
	if (fwrite(buffer,
		   sizeof(double),
		   block_size,
		   matrix_file) != block_size)
		error("Could not write matrix elements to matrix block %lu for"
		      " operator %s\n",
		      matrix_block_id,
		      operator_path);
	fclose(matrix_file);
}
