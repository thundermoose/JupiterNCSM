#include <matrix_builder/matrix_builder.h>
#include <string_tools/string_tools.h>
#include <diagonalization/diagonalization.h>
#include <log/log.h>
#include <error/error.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <string.h>
#include <unit_testing/test.h>
#include <combination_table/combination_table.h>
#include <execution_order/execution_order.h>
#include <scheduler/scheduler.h>

struct _matrix_builder_
{
	matrix_builder_settings_t settings;
	combination_table_t combination_table;	
	execution_order_t execution_order;
	scheduler_t scheduler;
};

matrix_builder_t new_matrix_builder(matrix_builder_settings_t settings)
{
	matrix_builder_t builder =
		(matrix_builder_t)malloc(sizeof(struct _matrix_builder_));
	builder->settings = settings;
	builder->combination_table =
		new_combination_table(settings.combination_file_path,
				      settings.num_protons,
				      settings.num_neutrons);
	builder->execution_order =
		read_execution_order(settings.minerva_instruction_path,
				     builder->combination_table);
	builder->scheduler = 
		new_scheduler(builder->execution_order,
			      builder->combination_table,
			      settings.index_list_path,
			      settings.interaction_path);		
	return builder;	
}

matrix_t generate_matrix(matrix_builder_t builder)
{
	vector_settings_t settings_template =
		setup_vector_settings(builder->combination_table);
	vector_settings_t input_vector_settings = settings_template;
	vector_settings_t output_vector_settings = settings_template;
	input_vector_settings.directory_name =
	       	copy_string(builder->settings.input_vector_path);
	output_vector_settings.directory_name =
	       	copy_string(builder->settings.output_vector_path);
	vector_t input_vector = new_zero_vector(input_vector_settings);
	matrix_t matrix = new_zero_matrix(vector_dimension(input_vector),
					  vector_dimension(input_vector));
	for (size_t i = 0; i < vector_dimension(input_vector); i++)
	{
		vector_t output_vector =
		       	new_zero_vector(output_vector_settings);
		set_element(input_vector,i,1.0);
		save_vector(input_vector);
		run_matrix_vector_multiplication(output_vector_settings.directory_name,
						 input_vector_settings.directory_name,
						 builder->scheduler);
		set_element(input_vector,i,0.0);	
		for (size_t j = 0; j < vector_dimension(output_vector); j++)
		{
			log_entry("element[%lu] = %lg",
				  j,
				  get_element(output_vector,j));
			set_matrix_element(matrix,
					   j,
					   i,
					   get_element(output_vector,j));
		}
		free_vector(output_vector);
	}
	free_vector(input_vector);
	free(settings_template.block_sizes);
	free(input_vector_settings.directory_name);
	free(output_vector_settings.directory_name);
	return matrix;
}

void free_matrix_builder(matrix_builder_t builder)
{
	free_combination_table(builder->combination_table);
	free_execution_order(builder->execution_order);
	free_scheduler(builder->scheduler);
	free(builder);
}

#define BACCHUS_RUN "bacchus_run_data/he4/"
new_test(explicit_matrix_nmax0,
	 const char *matrix_out_path =
	 get_test_file_path("hamiltonian.npy");
	 const char *input_vector_path =
	 	get_test_file_path("input");
	 const char *output_vector_path =
	 	get_test_file_path("output");
	 matrix_builder_settings_t settings =
	 {
		.combination_file_path = 
		TEST_DATA BACCHUS_RUN "nmax0/comb.txt",
		.minerva_instruction_path =
	       		TEST_DATA BACCHUS_RUN "nmax0/greedy_2_16.txt",
		.index_list_path = 
		TEST_DATA BACCHUS_RUN "nmax0/index_lists",
		.interaction_path = 
		TEST_DATA BACCHUS_RUN "nmax0/interaction",
		.input_vector_path =
			copy_string(input_vector_path),
		.output_vector_path =
			copy_string(output_vector_path)
	 };
	 matrix_builder_t matrix_builder = new_matrix_builder(settings);
	 matrix_t matrix = generate_matrix(matrix_builder);
	 free_matrix_builder(matrix_builder);
	 FILE *matrix_file = fopen(matrix_out_path,"w");
	 save_numpy_matrix(matrix_file, matrix);
		 fclose(matrix_file);
	 eigen_system_t eigen_system = diagonalize_symmetric_matrix(matrix);
	 size_t num_eigen_values = get_num_eigen_values(eigen_system);
	 for (size_t i = 0; i<num_eigen_values; i++)
	 	printf("(%lu): %0.16lg\n",i,
		       get_eigen_value(eigen_system,i));
	 free_eigen_system(eigen_system);
	 free_matrix(matrix);
	 free(settings.input_vector_path);
	 free(settings.output_vector_path);
	 );

new_test(explicit_matrix_nmax2,
	 const char *matrix_out_path =
	 	get_test_file_path("hamiltonian.npy");
	 const char *input_vector_path =
	 	get_test_file_path("input");
	 const char *output_vector_path =
	 	get_test_file_path("output");
	 matrix_builder_settings_t settings =
	 {
		.combination_file_path = TEST_DATA BACCHUS_RUN "nmax2/comb.txt",
		.minerva_instruction_path =
	       		TEST_DATA BACCHUS_RUN "nmax2/greedy_3_16.txt",
		.index_list_path = 
		TEST_DATA BACCHUS_RUN "nmax2/index_lists",
		.interaction_path = 
		TEST_DATA BACCHUS_RUN "nmax2/interaction",
		.input_vector_path =
			copy_string(input_vector_path),
		.output_vector_path =
			copy_string(output_vector_path)
	 };
	 matrix_builder_t matrix_builder = new_matrix_builder(settings);
	 matrix_t matrix = generate_matrix(matrix_builder);
	 free_matrix_builder(matrix_builder);
	 FILE *matrix_file = fopen(matrix_out_path,"w");
	 save_numpy_matrix(matrix_file, matrix);
	 fclose(matrix_file);
	 eigen_system_t eigen_system = diagonalize_symmetric_matrix(matrix);
	 size_t num_eigen_values = get_num_eigen_values(eigen_system);
	 for (size_t i = 0; i<num_eigen_values; i++)
	 	printf("(%lu): %lg\n",i,
		       get_eigen_value(eigen_system,i));
	 free_eigen_system(eigen_system);
	 free_matrix(matrix);
	 free(settings.input_vector_path);
	 free(settings.output_vector_path);
	 );

new_test(explicit_matrix_nmax4,
	 const char *matrix_out_path =
	 	get_test_file_path("hamiltonian.npy");
	 const char *input_vector_path =
	 	get_test_file_path("input");
	 const char *output_vector_path =
	 	get_test_file_path("output");
	 matrix_builder_settings_t settings =
	 {
		.combination_file_path = TEST_DATA BACCHUS_RUN "nmax4/comb.txt",
		.minerva_instruction_path =
	       		TEST_DATA BACCHUS_RUN "nmax4/greedy_2_16.txt",
		.index_list_path = 
		TEST_DATA BACCHUS_RUN "nmax4/index_lists",
		.interaction_path = 
		TEST_DATA BACCHUS_RUN "nmax4/interaction",
		.input_vector_path =
			copy_string(input_vector_path),
		.output_vector_path =
			copy_string(output_vector_path)
	 };
	 matrix_builder_t matrix_builder = new_matrix_builder(settings);
	 matrix_t matrix = generate_matrix(matrix_builder);
	 free_matrix_builder(matrix_builder);
	 FILE *matrix_file = fopen(matrix_out_path,"w");
	 save_numpy_matrix(matrix_file, matrix);
	 fclose(matrix_file);
	 eigen_system_t eigen_system = diagonalize_symmetric_matrix(matrix);
	 size_t num_eigen_values = get_num_eigen_values(eigen_system);
	 for (size_t i = 0; i<num_eigen_values; i++)
	 	printf("(%lu): %lg\n",i,
		       get_eigen_value(eigen_system,i));
	 free_eigen_system(eigen_system);
	 free_matrix(matrix);
	 free(settings.input_vector_path);
	 free(settings.output_vector_path);
	 );

new_test(single_matrix_vector_multiplication,
	 const char *input_vector_path =
	 	get_test_file_path("input");
	 const char *output_vector_path =
	 	get_test_file_path("output");
	 const char *execution_order_path =
		TEST_DATA BACCHUS_RUN "nmax4/greedy_2_16.txt";
	 const char *index_list_path =
                TEST_DATA BACCHUS_RUN "nmax4/index_lists";
         const char *interaction_path =
                TEST_DATA BACCHUS_RUN "nmax4/interaction";	 
	 const char *combination_table_path = 
	 	TEST_DATA BACCHUS_RUN "nmax4/comb.txt";
	const size_t num_protons = 2;
	const size_t num_neutrons = 2;
	const size_t state_in_index = 73;
	const size_t state_out_index = 108;
	combination_table_t combination_table =
		new_combination_table(combination_table_path,
				      num_protons,
				      num_neutrons);
	execution_order_t execution_order =
		read_execution_order(execution_order_path,
				     combination_table);
	scheduler_t scheduler = new_scheduler(execution_order,
				      combination_table,
				      index_list_path,
				      interaction_path);
	vector_settings_t input_settings =
       		setup_vector_settings(combination_table);
	input_settings.directory_name = copy_string(input_vector_path);
	vector_t input_vector = new_zero_vector(input_settings);
	vector_settings_t output_settings =
		setup_vector_settings(combination_table);
	output_settings.directory_name = copy_string(output_vector_path);
	vector_t output_vector = new_zero_vector(output_settings);
	set_element(input_vector,state_in_index,1.0);
	save_vector(input_vector);
	save_vector(output_vector);
	run_matrix_vector_multiplication(output_vector_path,
					 input_vector_path,
					 scheduler);
	printf("Element is %lg\n",
	       get_element(output_vector,state_out_index));
	free_vector(output_vector);
	free_vector(input_vector);
	free_scheduler(scheduler);
	free_execution_order(execution_order);
	free_combination_table(combination_table);
	free(input_settings.directory_name);
	free(output_settings.directory_name);
	);

