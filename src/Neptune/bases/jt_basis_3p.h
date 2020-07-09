#ifndef __JT_BASIS_3P__
#define __JT_BASIS_3P__

#include <bases/shells.h>
#include <bases/m_scheme_3p_basis.h>

typedef struct _jt_state_
{
	true_shell_index a,b,c;
	quantum_number jab,jabc,tab,tabc;
} JT_State;

typedef struct _jt_basis_
{
	JT_State *states;
	size_t dimension;
	Shells *shells;
} JT_Basis;

/* This method creates a jt_basis
 * with a minimum ho energy of e_min and
 * maximum e_max, with shells specified
 * in shells
 */
JT_Basis* new_jt_basis(quantum_number e_max,
		quantum_number e_min,
		Shells *shells);

/* This method creates a JT_Basis that
 * corresponds to a given m_scheme_basis
 */
JT_Basis* new_jt_basis_from_m_scheme(M_Scheme_3p_Basis* m_basis);


/* This method reads a JT_Basis from an ascii file
*/
JT_Basis* new_jt_basis_from_ascii_file(FILE* file,
		Shells *shells);

/* This method, extracts JT_States from jt_basis
 * that corresponds to j_abc and t_abc.
 */
JT_Basis* get_jt_block(JT_Basis* jt_basis,
		quantum_number j_abc,
		quantum_number t_abc);


/* This function lists the JT_States in a
 * given JT_Basis
 */
void print_jt_basis_3p(JT_Basis *jt_basis);

/* This function returns 1 if jt_basis contains
 * the same element more than one time otherwise 0
 */
int contains_duplicates(JT_Basis *jt_basis);


/* This method, generates an array with an
 * element for each state in jt_basis_a,
 * that is an index to the corresponding
 * state in jt_basis_b if such exists
 * otherwise it is -1
 */
ssize_t* matching_jt_states(JT_Basis* jt_basis_a,
		JT_Basis* jt_basis_b,
		ssize_t* a_to_b_shells);

void free_jt_basis_3p(JT_Basis* jt_basis);

#endif
