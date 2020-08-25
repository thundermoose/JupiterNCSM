#ifndef __TRANSFORM_BLOCK_SETTINGS__
#define __TRANSFORM_BLOCK_SETTINGS__

#include <matrix_block_setting/matrix_block_setting.h>

typedef struct
{
	int proton_energy_ket;
	int proton_energy_bra;
	int neutron_energy_ket;
	int neutron_energy_bra;
	int total_isospin;
} transform_block_settings_t;

transform_block_settings_t 
setup_transform_block(matrix_block_setting_t block_setting);

int compare_transform_block_settings(transform_block_settings_t *block_a,
				     transform_block_settings_t *block_b);
#endif
