#include <stdlib.h>
#include <stdio.h>
#include <arguments/arguments.h>
#include <transform_block_settings/transform_block_settings.h>
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
		return EXIT_SUCCESS;
	}
	combination_table_t combination_table =
		new_combination_table(get_combination_file_path(arguments),
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
	while (has_next_1nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_1nf_block_iterator(combination_table);
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
}

static
void generate_2nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	printf("2NF blocks only:\n");
	m_scheme_2p_basis_t current_ket_basis = NULL;
	m_scheme_2p_basis_t current_bra_basis = NULL;
	Dens_Matrix *current_transformed_matrix = NULL;
	transform_block_settings_t transformed_block = {0};
	size_t block_index = 0;
	while (has_next_2nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_2nf_block_iterator(combination_table);
		transform_block_settings_t current_block =
			setup_transform_block(current_matrix_block);
		printf("<(pE %d nE %d)|(pE %d nE %d)> Tz = %d\n",
		       current_block.proton_energy_bra,
		       current_block.neutron_energy_bra,
		       current_block.proton_energy_ket,
		       current_block.neutron_energy_ket,
		       current_block.total_isospin);
		if (compare_transform_block_settings
		    (&current_block,
		     &transformed_block) != 0 ||
		    block_index == 0)
		{
			if (current_transformed_matrix != NULL)
				free_dens_matrix(current_transformed_matrix);
			if (current_ket_basis != NULL)
				free_m_scheme_2p_basis(current_ket_basis);
			current_ket_basis = setup_ket_basis(current_block);
			current_bra_basis = setup_bra_basis(current_block);
			current_transformed_matrix = transform_matrix
		}
	}
}

static
void generate_3nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	printf("3NF blocks only:\n");
	while (has_next_3nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_3nf_block_iterator(combination_table);
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
}
