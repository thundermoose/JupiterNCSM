#ifndef __M_SCHEME_2P_BASIS__
#define __M_SCHEME_2P_BASIS__

#include "sp_states.h"
#include <stdlib.h>

struct _m_scheme_2p_basis_;
typedef struct _m_scheme_2p_basis_ *m_scheme_2p_basis_t;

typedef struct
{
	sp_state_index a,b;
} m_scheme_2p_state_t;

m_scheme_2p_basis_t new_m_scheme_2p_basis(quantum_number e_max1,
					  quantum_number e_max2);

m_scheme_2p_basis_t new_m_scheme_2p_basis_fixed_isospin(quantum_number e_max1,
							quantum_number e_max2,
							quantum_number Tz);

m_scheme_2p_basis_t 
new_m_scheme_2p_basis_from_files(quantum_number e_max1,
				 const char *proton_basis_filename,
				 const char *neutron_basis_filename);

m_scheme_2p_basis_t generate_2p_block(m_scheme_2p_basis_t m_scheme_2p_basis,
				      quantum_number Tz,
				      quantum_number M,
				      quantum_number E);

// Warning: The following function assumes that the original basis is 
// sorted in blocks with constant M. If that is not the case
// use generate_2p_block instead
m_scheme_2p_basis_t cut_out_M_block(m_scheme_2p_basis_t m_scheme_2p_basis,
				    quantum_number M);

size_t get_m_scheme_2p_dimension(m_scheme_2p_basis_t m_scheme_2p_basis);

quantum_number get_m_scheme_2p_e_max1(m_scheme_2p_basis_t m_scheme_2p_basis);

quantum_number get_m_scheme_2p_e_max2(m_scheme_2p_basis_t m_scheme_2p_basis);

Shells *get_m_scheme_shells(m_scheme_2p_basis_t m_scheme_2p_basis);

SP_States *get_m_scheme_sp_states(m_scheme_2p_basis_t m_scheme_2p_basis);

m_scheme_2p_state_t 
get_m_scheme_2p_state(m_scheme_2p_basis_t m_scheme_2p_basis,
		      size_t index);

m_scheme_2p_state_t*
get_m_scheme_2p_states(m_scheme_2p_basis_t m_scheme_2p_basis);

size_t *m_scheme_2p_corresponding_indices(m_scheme_2p_basis_t basis);

void print_m_scheme_2p_basis(m_scheme_2p_basis_t basis);

#ifndef NLOGING
void log_m_scheme_2p_basis(m_scheme_2p_basis_t basis);
#else
#define log_m_scheme_2p_basis(basis)
#endif

size_t get_m_scheme_2p_state_index(m_scheme_2p_basis_t basis,
				   m_scheme_2p_state_t state);

void free_m_scheme_2p_basis(m_scheme_2p_basis_t m_scheme_2p_basis);
#endif
