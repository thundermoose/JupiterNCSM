#include <stdlib.h>
#include <stdio.h>
#include <index_list/index_list.h>
#include <combination_table/combination_table.h>
#include <log/log.h>

static
void transform_file(const char *index_list_path,
		    const char *output_path,
		    index_list_setting_t setting);

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("AURORA_LOGFILE",
			 "aurora.log");
}

int main(int num_arguments,
	 char **argument_list)
{
	if (num_arguments < 6)
	{
		printf("Usage %s <comb.txt> <index_list_path> <output_path> <Z> <N>\n",
		       *argument_list);
		return EXIT_FAILURE;
	}
	const char *comb_file_path = argument_list[1];
	const char *index_list_path = argument_list[2];
	const char *output_path = argument_list[3];
	size_t num_protons = atoi(argument_list[4]);
	size_t num_neutrons = atoi(argument_list[5]);
	combination_table_t table = new_combination_table(comb_file_path,
							  num_protons,
							  num_neutrons);
	while (has_next_index_list_setting(table))
	{
		index_list_setting_t setting = next_index_list_setting(table);
		transform_file(index_list_path,
			       output_path,
			       setting);
	}
	free_combination_table(table);
	return EXIT_SUCCESS;
}

static
void transform_file(const char *index_list_path,
		    const char *output_path,
		    index_list_setting_t setting)
{

	char index_list_file_name[4096];
	sprintf(index_list_file_name,
		"%s/%s_inds_index_lists/index_list_E_in%d_E_out%d_M_in%d_M_out%d_dE%d_dM%d_depth%d",
		index_list_path,
		block_type_to_string(setting.block_type),
		setting.energy_bra,
		setting.energy_ket,
		setting.M_bra,
		setting.M_ket,
		setting.energy_ket-setting.energy_bra,
		setting.M_ket-setting.M_bra,
		setting.depth);
	index_list_t index_list =
		parse_human_readable_index_list(index_list_file_name);
	sprintf(index_list_file_name,
		"%s/index_list_%lu",
		output_path,
		setting.index_list_id);
	save_index_list(index_list,
			index_list_file_name);
	free_index_list(index_list);
}
