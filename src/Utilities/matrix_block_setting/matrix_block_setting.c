#include <matrix_block_setting/matrix_block_setting.h>

size_t get_matrix_block_length(matrix_block_setting_t settings)
{
	switch (settings.type)
	{
		case N_block:
		case NN_block:
		case NNN_block:
			return settings.num_neutron_combinations;
		case P_block:
		case PP_block:
		case PPP_block:
			return settings.num_proton_combinations;
		case NP_block:
		case NNP_block:
		case NPP_block:
			return settings.num_neutron_combinations*
				settings.num_proton_combinations;
		default:
			return 0;
	}
}
