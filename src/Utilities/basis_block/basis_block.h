#ifndef __BASIS_BLOCK__
#define __BASIS_BLOCK__

#include <stdlib.h>

typedef struct
{
	int Ep;
	int Mp;
	int En;
	int Mn;
	size_t num_proton_states;
	size_t num_neutron_states;
	size_t num_protons;
	size_t num_neutrons;
	size_t block_id;
} basis_block_t;

basis_block_t new_basis_block(int proton_energy,
			      int proton_total_M,
			      int neutron_energy,
			      int neutron_total_M,
			      size_t num_proton_states,
			      size_t num_neutron_states,
			      size_t block_id);

size_t get_basis_block_dimension(basis_block_t basis_block);

size_t get_num_protons(basis_block_t basis_block);

size_t get_num_neutrons(basis_block_t basis_block);

#endif
