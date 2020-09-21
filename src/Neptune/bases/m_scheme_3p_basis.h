#ifndef __M_SCHEME_3P_BASIS__
#define __M_SCHEME_3P_BASIS__
#include "sp_states.h"

typedef struct _m_scheme_3p_state_
{
	sp_state_index a,b,c;
} M_Scheme_3p_State;

M_Scheme_3p_State sort_on_shells(M_Scheme_3p_State s,
				 SP_States* sp_states);

typedef struct _m_scheme_3p_basis_
{
	M_Scheme_3p_State *states;
	size_t dimension;
	SP_States *sp_states;
	quantum_number e_max;
} M_Scheme_3p_Basis;

M_Scheme_3p_Basis* new_m_scheme_3p_basis_no_m_rest(quantum_number e_max,
						   SP_States *sp_states);

M_Scheme_3p_Basis* new_m_scheme_3p_basis(quantum_number e_max,
					 quantum_number m_tot,
					 SP_States *sp_states);

M_Scheme_3p_Basis* new_m_scheme_3p_basis2(quantum_number e_min,
					  quantum_number e_max,
					  quantum_number m_tot,
					  SP_States *sp_states);

M_Scheme_3p_Basis* new_m_scheme_3p_basis3(quantum_number e_max,
					  quantum_number m_tot,
					  quantum_number tz,
					  SP_States *sp_states);

M_Scheme_3p_Basis* 
new_m_scheme_3p_basis_from_basis_file(const char *proton_basis_file,
				      const char *neutron_basis_file,
				      size_t num_protons,
				      size_t num_neutrons,
				      SP_States *sp_states);


M_Scheme_3p_Basis* generate_block(M_Scheme_3p_Basis *mp_basis,
				  quantum_number tz,
				  quantum_number m_tot,
				  quantum_number E,
				  size_t *offset);

M_Scheme_3p_State* get_m_scheme_3p_states(M_Scheme_3p_Basis *mp_basis);

void list_m_scheme_3p_basis(M_Scheme_3p_Basis* mp_basis);

M_Scheme_3p_State new_m_scheme_3p_state(sp_state_index a,
					sp_state_index b,
					sp_state_index c);

void free_m_scheme_3p_basis(M_Scheme_3p_Basis* mp_basis);

#endif
