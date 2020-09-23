#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <arguments/arguments.h>
#include <input/read_2nf_antoine_format.h>
#include <transform_block_settings/transform_block_settings.h>
#include <transform_2nf_block_manager/transform_2nf_block_manager.h>
#include <transform_3nf_block_manager/transform_3nf_block_manager.h>
#include <combination_table/combination_table.h>
#include <log/log.h>

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MERCURY_LOGFILE",
			 "mercury_J_scheme_to_internal.log");	 
}

static
void generate_1nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

static
void generate_2nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

static
void generate_3nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments =
	       	parse_argument_list(num_arguments, argument_list);	    
	if (to_few_arguments(arguments))
	{
		show_usage(arguments);
		free_arguments(arguments);
		return EXIT_SUCCESS;
	}
	combination_table_t combination_table =
		new_combination_table
		(get_combination_file_path_argument(arguments),
		 get_num_protons_argument(arguments),
		 get_num_neutrons_argument(arguments));
	generate_1nf_matrix_blocks(combination_table,
				   arguments);
	generate_2nf_matrix_blocks(combination_table,
				   arguments);
	generate_3nf_matrix_blocks(combination_table,
				   arguments);
	free_combination_table(combination_table);
	free_arguments(arguments);
}

static
void generate_1nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	printf("1NF blocks only:\n");
	const char *output_path_base = get_output_path_argument(arguments);
	const char *index_list_path = get_index_list_path_argument(arguments);
	while (has_next_1nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_1nf_block_iterator(combination_table);
		connection_list_t connection_list =
			read_connection_files(index_list_path,
					      current_matrix_block);
		mercury_matrix_block_t current_block =
			new_zero_mercury_matrix_block(connection_list);
		save_mercury_matrix_block(current_block,
					  output_path_base);
		free_mercury_matrix_block(current_block);
		free_connection_list(connection_list);
	}
}

static
void generate_2nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	printf("2NF blocks only:\n");
	log_entry("num_particle_argument = %lu",
		  get_num_particles_argument(arguments));
	antoine_2nf_file_t coupled_2nf_data =
		open_antoine_2nf_file
		(get_interaction_path_2nf_argument(arguments),
		 get_num_particles_argument(arguments),
		 get_single_particle_energy_argument(arguments),
		 get_two_particle_energy_argument(arguments));
	transform_2nf_block_manager_t manager =
		new_transform_2nf_block_manager
		(coupled_2nf_data,
		 get_index_list_path_argument(arguments),
		 get_index_list_path_argument(arguments),
		 get_single_particle_energy_argument(arguments));
	transform_block_settings_t transformed_block = {INT_MAX};	
	const char *output_path_base = get_output_path_argument(arguments);
	while (has_next_2nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_2nf_block_iterator(combination_table);
		transform_block_settings_t current_block =
			setup_transform_block(current_matrix_block);
		log_entry("<(pE %d nE %d)|(pE %d nE %d)> Tz = %d",
			  current_block.proton_energy_bra,
			  current_block.neutron_energy_bra,
			  current_block.proton_energy_ket,
			  current_block.neutron_energy_ket,
			  current_block.total_isospin);
		if (compare_transform_block_settings(&transformed_block,
						     &current_block) != 0)
		{
			log_entry("Needs to decouple a new block");
			decouple_transform_2nf_block(manager,
						     current_block);
			transformed_block = current_block;
		}
		mercury_matrix_block_t matrix_block = 
			get_transform_2nf_matrix_block(manager,
						       current_matrix_block);
		log_entry("Retrieved the current matrix_block");
		save_mercury_matrix_block(matrix_block,
					  output_path_base);
		log_entry("Saved the current matrix_block to file");
		free_mercury_matrix_block(matrix_block);
	}
	free_transform_2nf_block_manager(manager);
	free_antoine_2nf_file(coupled_2nf_data);
}

static
void generate_3nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	printf("3NF blocks only:\n");
	Data_File *coupled_3nf_data =
		open_data_file(get_interaction_path_3nf_argument(arguments));
	transform_3nf_block_manager_t manager =
	       	new_transform_3nf_block_manager
		(coupled_3nf_data,
		 get_index_list_path_argument(arguments),
		 get_single_particle_energy_argument(arguments));
	transform_block_settings_t transformed_block = {INT_MAX};
	const char *output_path_base = get_output_path_argument(arguments);					
	while (has_next_3nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_3nf_block_iterator(combination_table);
		transform_block_settings_t current_block =
			setup_transform_block(current_matrix_block);
		log_entry("%s p: %d %d %d n: %d %d %d,"
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
		if (compare_transform_block_settings(&current_block,
						     &transformed_block) != 0)
		{
			decouple_transform_3nf_block(manager,
						     current_block);
		}
		mercury_matrix_block_t matrix_block =
			get_transform_3nf_matrix_block(manager,
						       current_matrix_block);
		save_mercury_matrix_block(matrix_block,
					  output_path_base);
		free_mercury_matrix_block(matrix_block);
	}
	free_transform_3nf_block_manager(manager);
	free_data_file(coupled_3nf_data);
}
