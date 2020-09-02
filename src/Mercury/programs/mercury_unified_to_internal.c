#include <stdlib.h>
#include <stdio.h>
#include <arguments/arguments.h>
#include <interaction/interaction.h>
#include <combination_table/combination_table.h>
#include <connection_list/connection_list.h>
#include <single_particle_basis/single_particle_basis.h>
#include <mercury_matrix_block/mercury_matrix_block.h>
#include <log/log.h>

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MERCURY_LOGFILE",
			 "mercury.log");
}

static
void process_block(matrix_block_setting_t matrix_block_setting,
		   single_particle_basis_t sp_basis,
		   interaction_t interaction,
		   const char *index_list_path,
		   const char *output_path)
{
	printf("Processing matrix block %lu\r",
	       matrix_block_setting.matrix_block_id);
	fflush(stdout);
	log_entry("Current block: %s %d %d %d %d %d %d %lu %lu %lu",
		  block_type_to_string(matrix_block_setting.type),
		  matrix_block_setting.difference_energy_protons,
		  matrix_block_setting.difference_M_protons,
		  matrix_block_setting.depth_protons,
		  matrix_block_setting.difference_energy_neutrons,
		  matrix_block_setting.difference_M_neutrons,
		  matrix_block_setting.depth_neutrons,
		  matrix_block_setting.num_proton_combinations,
		  matrix_block_setting.num_neutron_combinations,
		  matrix_block_setting.matrix_block_id);
	connection_list_t connection_list =
		read_connection_files(index_list_path,
				      matrix_block_setting);
	mercury_matrix_block_t matrix_block =
		new_mercury_matrix_block(interaction,
					 connection_list,
					 sp_basis);
	save_mercury_matrix_block(matrix_block,
				  output_path);
	free_mercury_matrix_block(matrix_block);
	free_connection_list(connection_list);
	printf("Processing matrix block %lu: Done\n",
	       matrix_block_setting.matrix_block_id);
}

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments = parse_argument_list(num_arguments,
						    argument_list);
	if (to_few_arguments(arguments))
	{
		show_usage(arguments);
		free_arguments(arguments);
		return EXIT_FAILURE;
	}
	combination_table_t table =
		new_combination_table
		(get_combination_file_path_argument(arguments),
		 get_num_protons_argument(arguments),
		 get_num_neutrons_argument(arguments));
	log_entry("Parsed the combination table correctly\n"); 
	single_particle_basis_t sp_basis =
	       	new_single_particle_basis
		(get_single_particle_energy_argument(arguments));
	log_entry("Generated the single particle basis correctly\n");
	interaction_t interaction_2nf =
			new_interaction
			(get_interaction_path_2nf_argument(arguments));	
	interaction_t interaction_3nf = NULL;
	if (get_interaction_path_3nf_argument(arguments) != NULL)
		interaction_3nf = 
			new_interaction
			(get_interaction_path_3nf_argument(arguments));
	log_entry("Initiated the interaction correctly\n");
	if (single_block_mode(arguments))
	{
		matrix_block_setting_t matrix_block_setting =
			get_matrix_block_by_id(table,get_block_id(arguments));
		size_t num_particles =
			count_particles(matrix_block_setting.type);
		interaction_t interaction = 
			num_particles == 3 && 
			interaction_3nf != NULL?
			interaction_3nf : 
			interaction_2nf;
		process_block(matrix_block_setting,
			      sp_basis,
			      interaction,
			      get_index_list_path_argument(arguments),
			      get_output_path_argument(arguments));
	}
	else
	{
		while (has_next_matrix_block(table))
		{
			matrix_block_setting_t matrix_block_setting =
				next_matrix_block(table);
			size_t num_particles =
				count_particles(matrix_block_setting.type);
			interaction_t interaction = 
				num_particles == 3 &&
				interaction_3nf != NULL?
				interaction_3nf : 
				interaction_2nf;
			process_block(matrix_block_setting,
				      sp_basis,
				      interaction,
				      get_index_list_path_argument(arguments),
				      get_output_path_argument(arguments));
		}
	}
	free_interaction(interaction_2nf);
	if (interaction_3nf)
		free_interaction(interaction_3nf);
	free_combination_table(table);
	free_arguments(arguments);
	free_single_particle_basis(sp_basis);
	return EXIT_SUCCESS;
}
