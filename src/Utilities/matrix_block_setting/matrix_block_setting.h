#ifndef __MATRIX_BLOCK_SETTING__
#define __MATRIX_BLOCK_SETTING__

#include <stdlib.h>
#include <block_type/block_type.h>

typedef struct
{
	block_type_t type;
	int difference_energy_protons;	
	int difference_M_protons;
	int depth_protons;
	int difference_energy_neutrons;	
	int difference_M_neutrons;
	int depth_neutrons;
	size_t num_proton_combinations;
	size_t num_neutron_combinations;
	size_t matrix_block_id;
} matrix_block_setting_t;

size_t get_matrix_block_length(matrix_block_setting_t settings);

#endif
