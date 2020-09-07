#include <basis_block/basis_block.h>
#include <debug_mode/debug_mode.h>

basis_block_t new_basis_block(int proton_energy,
			      int proton_total_M,
			      int neutron_energy,
			      int neutron_total_M,
			      size_t num_proton_states,
			      size_t num_neutron_states,
			      size_t block_id)
{
	basis_block_t basis_block =
	{
		.Ep = proton_energy,
		.Mp = proton_total_M,
		.En = neutron_energy,
		.Mn = neutron_total_M,
		.num_proton_states = num_proton_states,
		.num_neutron_states = num_neutron_states,
		.block_id = block_id
	};
	return basis_block;
}

size_t get_basis_block_dimension(basis_block_t basis_block)
{
	return basis_block.num_proton_states*basis_block.num_neutron_states;
}

size_t get_num_protons(basis_block_t basis_block)
{
	return basis_block.num_protons;
}

size_t get_num_neutrons(basis_block_t basis_block)
{
	return basis_block.num_neutrons;
}
