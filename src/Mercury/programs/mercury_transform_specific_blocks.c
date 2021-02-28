#include <stdlib.h>
#include <stdio.h>
#include <transform_block_settings/transform_block_settings.h>
#include <transform_3nf_block_manager/transform_3nf_block_manager.h>
#include <combination_table/combination_table.h>
#include <iterator/iterator.h>
#include <matrix_energy_block/matrix_energy_block.h>
#include <input/read_3nf_hdf5_file.h>

int main(int num_arguments,
	 char **argument_list)
{
	if (num_arguments < 13)
	{
		printf("Usage: %s <comb.txt> <Z> <N> <interaction> <index_lists> <output> <Nmax1> <Nmax2> <C1> <C3> <C4> <CD> <CE> <block>...\n",
		       *argument_list);
		return EXIT_FAILURE;
	}
	const char *combfile_path = argument_list[1];	
	size_t num_protons = atoll(argument_list[2]);
	size_t num_neutrons = atoll(argument_list[3]);
	const char *interaction_path = argument_list[4];
	const char *index_list_path = argument_list[5];
	const char *output_path = argument_list[6];
	size_t single_particle_energy = atoll(argument_list[7]);
	double lec_C1 = atof(argument_list[8]);
	double lec_C3 = atof(argument_list[9]);
	double lec_C4 = atof(argument_list[10]);
	double lec_CD = atof(argument_list[11]);
	double lec_CE = atof(argument_list[12]);
	size_t num_blocks = num_arguments-13;
	size_t *block_numbers = (size_t *)malloc(num_blocks*sizeof(size_t));
	printf("comb: %s Z: %lu N: %lu\n"
	       "interaction: %s\n"
	       "index_list_path: %s\n"
	       "output_path: %s\n"
	       "Nmax: %lu\n"
	       "lecs: %lg %lg %lg %lg %lg\n",
	       combfile_path,
	       num_protons,
	       num_neutrons,
	       interaction_path,
	       index_list_path,
	       output_path,
	       single_particle_energy,
	       lec_C1,
	       lec_C3,
	       lec_C4,
	       lec_CD,
	       lec_CE);

	for (size_t i = 0; i<num_blocks; i++)
		block_numbers[i] = atof(argument_list[i+13]);
	combination_table_t table = new_combination_table(combfile_path,
							  num_protons,
							  num_neutrons);
	Data_File *coupled_3nf_data =
		open_data_file(interaction_path);
	set_max_loaded_memory(coupled_3nf_data, (size_t)(70)<<30);
	set_weight(coupled_3nf_data,CE,lec_CE);
	set_weight(coupled_3nf_data,CD,lec_CD);
	set_weight(coupled_3nf_data,C1,lec_C1);
	set_weight(coupled_3nf_data,C3,lec_C3);
	set_weight(coupled_3nf_data,C4,lec_C4);
	transform_3nf_block_manager_t manager = 
		new_transform_3nf_block_manager(coupled_3nf_data,
						index_list_path,
						single_particle_energy);
	for (size_t i = 0; i<num_blocks; i++)
	{
		matrix_block_setting_t current_block =
			get_matrix_block_by_id(table,block_numbers[i]);
		transform_block_settings_t transform_settings =
			setup_transform_block(current_block);	
		decouple_transform_3nf_block(manager,transform_settings);
		mercury_matrix_block_t matrix_block =
			get_transform_3nf_matrix_block(manager,
						       current_block);
		save_mercury_matrix_block(matrix_block,output_path);
		free_mercury_matrix_block(matrix_block);
	}
	free_transform_3nf_block_manager(manager);
	free_data_file(coupled_3nf_data);
	free_combination_table(table);
	free(block_numbers);
	return EXIT_SUCCESS;
}
