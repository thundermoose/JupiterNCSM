#ifndef __JJJ_COUPLED_3P__
#define __JJJ_COUPLED_3P__
#include <bases/shells.h>
#include <bases/m_scheme_3p_basis.h>


typedef struct _jjj_state_
{ 
  shell_index a,b,c;
  quantum_number j_ab,j_abc,tz,parity;
} JJJ_State;

typedef struct _jjj_basis_
{
  JJJ_State* states;
  size_t dimension;
  Shells* shells;
} JJJ_Basis;

JJJ_Basis* new_jjj_basis(quantum_number e_max,
			 Shells* shells);

JJJ_Basis* new_jjj_basis_m_scheme(M_Scheme_3p_Basis* ms_basis);

JJJ_Basis* new_jjj_basis_from_ascii(FILE* config_file,
				    Shells* shells);

/* Warning the following function do assumption on the order of the states*/
JJJ_Basis* get_block(JJJ_Basis* in_jjj_basis,
		     quantum_number j_abc,
		     quantum_number tz,
		     quantum_number parity);

void list_jjj_states(JJJ_Basis* jjj_basis);

ssize_t* matching_jjj_states(JJJ_Basis *jjj_basis_a,
			     JJJ_Basis *jjj_basis_b,
			     sshell_index* a_to_b_shells);

void free_jjj_basis(JJJ_Basis* jjj_basis);

#endif
