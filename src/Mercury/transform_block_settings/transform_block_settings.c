#include <transform_block_settings/transform_block_settings.h>

transform_block_settings_t 
setup_transform_block(matrix_block_setting_t block_setting)
{
	transform_block_settings_t block =
	{

		.proton_energy_ket = block_setting.depth_protons,	
		.proton_energy_bra = block_setting.depth_protons +
		       	block_setting.difference_energy_protons,
		.neutron_energy_ket = block_setting.depth_neutrons,	
		.neutron_energy_bra = block_setting.depth_neutrons +
		       	block_setting.difference_energy_neutrons,
		.total_isospin = count_neutrons(block_setting.type) -
			count_protons(block_setting.type)
	};
	return block;
}

int compare_transform_block_settings(transform_block_settings_t *block_a,
				     transform_block_settings_t *block_b)
{
	int diff = block_a->proton_energy_ket - block_b->proton_energy_ket;
	if (diff)
		return diff;
	diff = block_a->proton_energy_bra - block_b->proton_energy_bra;
	if (diff)
		return diff;
	diff = block_a->neutron_energy_ket - block_b->proton_energy_ket;
	if (diff)
		return diff;
	diff = block_a->neutron_energy_bra - block_b->proton_energy_bra;
	if (diff)
		return diff;
	return block_a->total_isospin - block_b->total_isospin;
}
