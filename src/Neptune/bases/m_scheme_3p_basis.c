#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "m_scheme_3p_basis.h"
#include <string_tools/string_tools.h>
#include <utils/permutation_tools.h>
#include <read_packed_states/read_packed_states.h>
#include <debug_mode/debug_mode.h>
#include <global_constants/global_constants.h>
#include <unit_testing/test.h>
#include <log/log.h>

// Evil global variables
static
SP_States *current_sp_states = NULL;

static
void read_single_particle_type_file(M_Scheme_3p_Basis *basis,
				    const char *basis_filename,
				    quantum_number iso_spin);

static
void read_two_particle_type_files(M_Scheme_3p_Basis *basis,
				  const char *proton_basis_filename,
				  const char *neutron_basis_filename,
				  const size_t num_protons,
				  const size_t num_neutrons);

static
void setup_state_to_index_map(M_Scheme_3p_Basis *basis);

static
int compare_m_scheme_3p_state(M_Scheme_3p_State *state_a,
			      M_Scheme_3p_State *state_b);

static
uint64_t compute_m_scheme_3p_hash(M_Scheme_3p_State *state);

static
void sort_state(M_Scheme_3p_State *state);

static
int compare_m_scheme_3p_states_m(M_Scheme_3p_State *state_a,
				 M_Scheme_3p_State *state_b);


M_Scheme_3p_State sort_on_shells(M_Scheme_3p_State s,
				 SP_States* sp_states)
{
	if (sp_states->sp_states[s.a].shell>sp_states->sp_states[s.b].shell)
	{
		SWAP(s.a,s.b);
	}
	if (sp_states->sp_states[s.b].shell>sp_states->sp_states[s.c].shell)
	{
		SWAP(s.b,s.c);
	}
	if (sp_states->sp_states[s.a].shell>sp_states->sp_states[s.b].shell)
	{
		SWAP(s.a,s.b);
	}
	return s;
}




int comp_m_scheme_3p_basis(const void* st1,
			   const void* st2,
			   void* sp_states_pt)
{
	SP_States* sp_states = (SP_States*)sp_states_pt;
	M_Scheme_3p_State* state_1 = (M_Scheme_3p_State*)st1;
	M_Scheme_3p_State* state_2 = (M_Scheme_3p_State*)st2;
	quantum_number Tz1 =
		(sp_states->shells->shells[sp_states->sp_states[state_1->a].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[state_1->b].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[state_1->c].shell].tz);
	quantum_number Tz2 =
		(sp_states->shells->shells[sp_states->sp_states[state_2->a].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[state_2->b].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[state_2->c].shell].tz);
	quantum_number diff = Tz1-Tz2;
	if (diff != 0)
		return diff;
	quantum_number M1 =
		(sp_states->sp_states[state_1->a].m+
		 sp_states->sp_states[state_1->b].m+
		 sp_states->sp_states[state_1->c].m);
	quantum_number M2 =
		(sp_states->sp_states[state_2->a].m+
		 sp_states->sp_states[state_2->b].m+
		 sp_states->sp_states[state_2->c].m);
	diff = M1-M2;
	if (diff != 0)
		return diff;
	quantum_number E1 =
		(sp_states->shells->shells[sp_states->sp_states[state_1->a].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[state_1->b].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[state_1->c].shell].e);
	quantum_number E2 =
		(sp_states->shells->shells[sp_states->sp_states[state_2->a].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[state_2->b].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[state_2->c].shell].e);
	diff = E1-E2;
	if (diff != 0)
		return diff;
	diff = state_1->a-state_2->a;
	if (diff != 0)
		return diff;
	diff = state_1->b-state_2->b;
	if (diff != 0)
		return diff;
	return state_1->c-state_2->c;
}


M_Scheme_3p_Basis* new_m_scheme_3p_basis_no_m_rest(quantum_number e_max,
						   SP_States* sp_states)
{
	M_Scheme_3p_Basis* mp_basis =
		(M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	mp_basis->sp_states = sp_states;
	mp_basis->e_max = e_max;
	mp_basis->states = NULL;
	mp_basis->dimension = 0;
	size_t max_dimension = 0;

	M_Scheme_3p_State current;
	for (current.a = 0;
	     current.a <sp_states->dimension-2;
	     current.a++)
	{
		SP_State a_state = sp_states->sp_states[current.a];
		Shell a_shell = sp_states->shells->shells[a_state.shell];
		if (a_shell.e>e_max)
			break;
		for (current.b = current.a+1;
		     current.b<sp_states->dimension-1;
		     current.b++)
		{
			SP_State b_state = sp_states->sp_states[current.b];
			Shell b_shell = sp_states->shells->shells[b_state.shell];
			if (b_shell.e+a_shell.e>e_max)
				break;
			for (current.c = current.b+1;
			     current.c <sp_states->dimension;
			     current.c++)
			{
				SP_State c_state = sp_states->sp_states[current.c];
				Shell c_shell = sp_states->shells->shells[c_state.shell];
				if (c_shell.e+b_shell.e+a_shell.e>e_max)
					break;



				if (max_dimension == mp_basis->dimension)
				{
					max_dimension=2*max_dimension+1;
					mp_basis->states =
						(M_Scheme_3p_State*)realloc(mp_basis->states,
									    sizeof(M_Scheme_3p_State)*max_dimension);
				}
				mp_basis->states[mp_basis->dimension++] = current;
			}
		}
	}
	mp_basis->states =
		(M_Scheme_3p_State*)realloc(mp_basis->states,
					    sizeof(M_Scheme_3p_State)*mp_basis->dimension);
	qsort_r(mp_basis->states,
		mp_basis->dimension,
		sizeof(M_Scheme_3p_State),
		comp_m_scheme_3p_basis,
		sp_states);
	setup_state_to_index_map(mp_basis);
	return mp_basis;
}


M_Scheme_3p_Basis* new_m_scheme_3p_basis(quantum_number e_max,
					 quantum_number m_tot,
					 SP_States* sp_states)
{
	M_Scheme_3p_Basis* mp_basis =
		(M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	mp_basis->sp_states = sp_states;

	mp_basis->states = NULL;
	mp_basis->dimension = 0;
	size_t max_dimension = 0;

	M_Scheme_3p_State current;
	for (current.a = 0;
	     current.a <sp_states->dimension-2;
	     current.a++)
	{
		SP_State a_state = sp_states->sp_states[current.a];
		Shell a_shell = sp_states->shells->shells[a_state.shell];
		if (a_shell.e>e_max)
			break;
		for (current.b = current.a+1;
		     current.b<sp_states->dimension-1;
		     current.b++)
		{
			SP_State b_state = sp_states->sp_states[current.b];
			Shell b_shell = sp_states->shells->shells[b_state.shell];
			if (b_shell.e+a_shell.e>e_max)
				break;
			for (current.c = current.b+1;
			     current.c <sp_states->dimension;
			     current.c++)
			{
				SP_State c_state = sp_states->sp_states[current.c];
				Shell c_shell = sp_states->shells->shells[c_state.shell];
				if (c_shell.e+b_shell.e+a_shell.e>e_max)
					break;
				if (a_state.m+b_state.m+c_state.m != m_tot)
					continue;

				if (max_dimension == mp_basis->dimension)
				{
					max_dimension=2*max_dimension+1;
					mp_basis->states =
						(M_Scheme_3p_State*)realloc(mp_basis->states,
									    sizeof(M_Scheme_3p_State)*max_dimension);
				}
				mp_basis->states[mp_basis->dimension++] = current;
			}
		}
	}
	mp_basis->states =
		(M_Scheme_3p_State*)realloc(mp_basis->states,
					    sizeof(M_Scheme_3p_State)*mp_basis->dimension);
	qsort_r(mp_basis->states,mp_basis->dimension,sizeof(M_Scheme_3p_State),
		comp_m_scheme_3p_basis,sp_states);
	setup_state_to_index_map(mp_basis);
	return mp_basis;
}

M_Scheme_3p_Basis* new_m_scheme_3p_basis2(quantum_number e_min,
					  quantum_number e_max,
					  quantum_number m_tot,
					  SP_States* sp_states)
{
	M_Scheme_3p_Basis* mp_basis =
		(M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	mp_basis->sp_states = sp_states;

	mp_basis->states = NULL;
	mp_basis->dimension = 0;
	size_t max_dimension = 0;

	M_Scheme_3p_State current;
	for (current.a = 0;
	     current.a <sp_states->dimension-2;
	     current.a++)
	{
		SP_State a_state = sp_states->sp_states[current.a];
		Shell a_shell = sp_states->shells->shells[a_state.shell];
		if (a_shell.e>e_max)
			break;
		for (current.b = current.a+1;
		     current.b<sp_states->dimension-1;
		     current.b++)
		{
			SP_State b_state = sp_states->sp_states[current.b];
			Shell b_shell = sp_states->shells->shells[b_state.shell];
			if (b_shell.e+a_shell.e>e_max)
				break;
			for (current.c = current.b+1;
			     current.c <sp_states->dimension;
			     current.c++)
			{
				SP_State c_state = sp_states->sp_states[current.c];
				Shell c_shell = sp_states->shells->shells[c_state.shell];
				if (c_shell.e+b_shell.e+a_shell.e>e_max)
					break;
				if (a_state.m+b_state.m+c_state.m != m_tot)
					continue;

				if (max_dimension == mp_basis->dimension)
				{
					max_dimension=2*max_dimension+1;
					mp_basis->states =
						(M_Scheme_3p_State*)realloc(mp_basis->states,
									    sizeof(M_Scheme_3p_State)*max_dimension);
				}
				mp_basis->states[mp_basis->dimension++] = current;
			}
		}
	}
	mp_basis->states =
		(M_Scheme_3p_State*)realloc(mp_basis->states,
					    sizeof(M_Scheme_3p_State)*mp_basis->dimension);
	qsort_r(mp_basis->states,mp_basis->dimension,sizeof(M_Scheme_3p_State),
		comp_m_scheme_3p_basis,sp_states);
	setup_state_to_index_map(mp_basis);
	return mp_basis;
}


void free_m_scheme_3p_basis(M_Scheme_3p_Basis* mp_basis)
{
	free(mp_basis->states);
	free_hash_map(mp_basis->state_to_index_map);
	free(mp_basis);
}

int get_3p_M(M_Scheme_3p_State state,
	     SP_States *sp_states)
{
	return sp_states->sp_states[state.a].m +
		sp_states->sp_states[state.b].m +
		sp_states->sp_states[state.c].m;
}


M_Scheme_3p_Basis* new_m_scheme_3p_basis3(quantum_number e_max,
					  quantum_number m_tot,
					  quantum_number tz,
					  SP_States* sp_states)
{
	M_Scheme_3p_Basis* mp_basis =
		(M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	mp_basis->sp_states = sp_states;

	mp_basis->states = NULL;
	mp_basis->dimension = 0;
	size_t max_dimension = 0;

	M_Scheme_3p_State current;
	for (current.a = 0;
	     current.a <sp_states->dimension-2;
	     current.a++)
	{
		SP_State a_state = sp_states->sp_states[current.a];
		Shell a_shell = sp_states->shells->shells[a_state.shell];
		if (a_shell.e>e_max)
			break;
		for (current.b = current.a+1;
		     current.b<sp_states->dimension-1;
		     current.b++)
		{
			SP_State b_state = sp_states->sp_states[current.b];
			Shell b_shell = sp_states->shells->shells[b_state.shell];
			if (b_shell.e+a_shell.e>e_max)
				break;
			for (current.c = current.b+1;
			     current.c <sp_states->dimension;
			     current.c++)
			{
				SP_State c_state = sp_states->sp_states[current.c];
				Shell c_shell = sp_states->shells->shells[c_state.shell];
				if (c_shell.e+b_shell.e+a_shell.e>e_max)
					break;
				if (a_state.m+b_state.m+c_state.m != m_tot)
					continue;
				if (a_shell.tz+b_shell.tz+c_shell.tz != tz)
					continue;

				if (max_dimension == mp_basis->dimension)
				{
					max_dimension=2*max_dimension+1;
					mp_basis->states =
						(M_Scheme_3p_State*)realloc(mp_basis->states,
									    sizeof(M_Scheme_3p_State)*max_dimension);
				}
				mp_basis->states[mp_basis->dimension++] = current;
			}
		}
	}
	mp_basis->states =
		(M_Scheme_3p_State*)realloc(mp_basis->states,
					    sizeof(M_Scheme_3p_State)*mp_basis->dimension);
	qsort_r(mp_basis->states,mp_basis->dimension,sizeof(M_Scheme_3p_State),
		comp_m_scheme_3p_basis,sp_states);
	setup_state_to_index_map(mp_basis);
	return mp_basis;
}


M_Scheme_3p_Basis* 
new_m_scheme_3p_basis_from_basis_file(const char *proton_basis_filename,
				      const char *neutron_basis_filename,
				      size_t num_protons,
				      size_t num_neutrons,
				      SP_States *sp_states)
{
	quantum_number iso_spin = (int)(num_neutrons) - (int)(num_protons);
	M_Scheme_3p_Basis *basis =
	       	(M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	basis->sp_states = sp_states;
	if (abs(iso_spin) == 3)
	{
		const char *basis_filename =
		       	iso_spin < 0 ? 
			proton_basis_filename :
		       	neutron_basis_filename;
		read_single_particle_type_file(basis,
					       basis_filename,
					       iso_spin);
	}
	else
	{
		read_two_particle_type_files(basis,
					     proton_basis_filename,
					     neutron_basis_filename,
					     num_protons,
					     num_neutrons);
	}
	setup_state_to_index_map(basis);
	return basis;
}


//inline
int diff_from_block(M_Scheme_3p_State current,
		    quantum_number tz,
		    quantum_number m_tot,
		    quantum_number E,
		    SP_States* sp_states)
{
	quantum_number Tz =
		(sp_states->shells->shells[sp_states->sp_states[current.a].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[current.b].shell].tz+
		 sp_states->shells->shells[sp_states->sp_states[current.c].shell].tz);
	quantum_number diff = Tz-tz;
	if (diff != 0)
		return diff;
	quantum_number M =
		(sp_states->sp_states[current.a].m+
		 sp_states->sp_states[current.b].m+
		 sp_states->sp_states[current.c].m);
	diff = M-m_tot;
	if (diff != 0)
		return diff;
	quantum_number Ec =
		(sp_states->shells->shells[sp_states->sp_states[current.a].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[current.b].shell].e+
		 sp_states->shells->shells[sp_states->sp_states[current.c].shell].e);
	return Ec-E;
}

int find_block_limits(M_Scheme_3p_Basis *mp_basis,
		      quantum_number tz,
		      quantum_number m_tot,
		      quantum_number E,
		      size_t *start,
		      size_t *stop)
{
	// Find an element in the block using binary search
	ssize_t known_element;
	ssize_t lower_limit = -1;
	ssize_t upper_limit = mp_basis->dimension;
	do
	{
		known_element = (upper_limit+lower_limit)>>1;
		M_Scheme_3p_State current_state
			= mp_basis->states[known_element];
		int diff = diff_from_block(current_state,
					   tz,m_tot,E,
					   mp_basis->sp_states);
		if (diff>0)
		{
			upper_limit = known_element;
		}
		else if(diff<0)
		{
			lower_limit = known_element;
		}
		else
		{
			break;
		}
	}
	while (upper_limit-lower_limit>1);
	M_Scheme_3p_State current_state
		= mp_basis->states[known_element];
	int diff = diff_from_block(current_state,
				   tz,m_tot,E,
				   mp_basis->sp_states);
	if (diff !=0)
	{
		return 1;
	}

	/* Using the known block element and the limits
	 * from the previous step to find the outer edges
	 * of the block using binary search
	 */ 

	// Find the lower edge
	ssize_t lower_edge;
	ssize_t lower_edge_max = known_element;
	do
	{
		lower_edge = (lower_edge_max+lower_limit)>>1;
		if (lower_edge<0)
		{
			lower_edge = 0;
		}
		M_Scheme_3p_State current_state
			= mp_basis->states[lower_edge];
		int diff = diff_from_block(current_state,
					   tz,m_tot,E,
					   mp_basis->sp_states);
		if (diff<0)
		{
			lower_limit = lower_edge;
		}
		else
		{
			lower_edge_max = lower_edge;
		}

	}
	while (lower_edge_max-lower_limit>1);
	*start= lower_edge_max;

	// Find the upper edge
	ssize_t upper_edge = known_element;
	while (upper_edge<mp_basis->dimension &&
	       diff_from_block(mp_basis->states[upper_edge],
			       tz,m_tot,E,
			       mp_basis->sp_states)==0)
	{
		upper_edge++;
	}

	/* Here resides the bug, it seems like that *stop
	 * sometimes is "inside" the block, but we assume
	 * it to be precisely "outside" it.
	 * I have tried several different ways to deal with
	 * this problem but for the moment nothing seems
	 * to work. Perhaps a partial rewrite is necessary
	 */


	*stop = upper_edge;

	return 0;
}

M_Scheme_3p_State* get_m_scheme_3p_states(M_Scheme_3p_Basis *mp_basis)
{
	M_Scheme_3p_State* states =
		(M_Scheme_3p_State*)
		malloc(mp_basis->dimension*sizeof(M_Scheme_3p_State));
	memcpy(states,
	       mp_basis->states,
	       sizeof(M_Scheme_3p_State)*mp_basis->dimension);
	return states;
}

size_t get_m_scheme_3p_dimension(M_Scheme_3p_Basis *basis)
{
	return basis->dimension;
}

void list_m_scheme_3p_basis(M_Scheme_3p_Basis* mp_basis)
{
	size_t i;
	for (i = 0; i < mp_basis->dimension; i++)
	{
		M_Scheme_3p_State s = mp_basis->states[i];
		printf("(%ld): %d#{%ld,%d} %d#{%ld,%d} %d#{%ld,%d}\n",
		       i,
		       s.a,
		       mp_basis->sp_states->sp_states[s.a].shell,
		       mp_basis->sp_states->sp_states[s.a].m,
		       s.b,
		       mp_basis->sp_states->sp_states[s.b].shell,
		       mp_basis->sp_states->sp_states[s.b].m,
		       s.c,
		       mp_basis->sp_states->sp_states[s.c].shell,
		       mp_basis->sp_states->sp_states[s.c].m);
	}

}


M_Scheme_3p_Basis* generate_block(M_Scheme_3p_Basis *mp_basis,
				  quantum_number tz,
				  quantum_number m_tot,
				  quantum_number E,
				  size_t* offset)
{
	// Find block limits
	size_t start;
	size_t stop;

	if (find_block_limits(mp_basis,
			      tz,m_tot,E,
			      &start,
			      &stop))
	{
		return NULL;
	}

	// Setup M_Scheme... structure 
	M_Scheme_3p_Basis* out 
		= (M_Scheme_3p_Basis*)malloc(sizeof(M_Scheme_3p_Basis));
	out->dimension = stop-start;
	out->states
		= (M_Scheme_3p_State*)malloc(sizeof(M_Scheme_3p_State)
					     *out->dimension);
	out->sp_states = mp_basis->sp_states;
	out->e_max = E;

	// Copy the states
	memcpy(out->states,
	       mp_basis->states+start,
	       sizeof(M_Scheme_3p_State)*
	       out->dimension);
	setup_state_to_index_map(out);
	*offset = start;
	return out;
}

size_t find_m_scheme_3p_state(M_Scheme_3p_Basis *basis,
			      M_Scheme_3p_State state)
{
	size_t index = no_index;
	if (hash_map_get(basis->state_to_index_map, &state, &index))
		return index;
	else
		return no_index;
}

M_Scheme_3p_State new_m_scheme_3p_state(sp_state_index a,
					sp_state_index b,
					sp_state_index c)
{
	M_Scheme_3p_State state =
	{
		.a = a,
		.b = b,
		.c = c
	};
	return state;
}

static
void read_single_particle_type_file(M_Scheme_3p_Basis *basis,
				    const char *basis_filename,
				    quantum_number iso_spin)
{
	short *packed_states = NULL;
	size_t num_states = 0;
	read_packed_states((void**)&packed_states,
			   &num_states,
			   3*sizeof(short),
			   basis_filename);
	basis->dimension = num_states;
	basis->states =
	       	(M_Scheme_3p_State*)
		malloc(num_states*sizeof(M_Scheme_3p_State));
	const sp_state_index particle_stride = iso_spin > 0 ? 1 : 0;
	for (size_t i = 0; i < num_states; i++)
	{
		basis->states[i].a = 2*packed_states[3*i] + particle_stride; 
		basis->states[i].b = 2*packed_states[3*i+1] + particle_stride; 
		basis->states[i].c = 2*packed_states[3*i+2] + particle_stride; 
	}		
	free(packed_states);
}

static
void read_two_particle_type_files(M_Scheme_3p_Basis *basis,
				  const char *proton_basis_filename,
				  const char *neutron_basis_filename,
				  const size_t num_protons,
				  const size_t num_neutrons)
{
	short *packed_proton_states = NULL;
	size_t num_proton_states = 0;
	read_packed_states((void**)&packed_proton_states,
			   &num_proton_states,
			   num_protons*sizeof(short),
			   proton_basis_filename);
	short *packed_neutron_states = NULL;
	size_t num_neutron_states = 0;
	read_packed_states((void**)&packed_neutron_states,
			   &num_neutron_states,
			   num_neutrons*sizeof(short),
			   neutron_basis_filename);	
	log_entry("Read %lu protons states and %lu neutron states",
		  num_proton_states,num_neutron_states);
	const size_t num_states = num_proton_states*num_neutron_states;
	basis->states = 
		(M_Scheme_3p_State*)
		malloc(num_states*sizeof(M_Scheme_3p_State));
	basis->dimension = num_states;
	size_t state_index = 0;
	uint64_t proton_mask = 
		num_protons > 1 ? 0x00000000FFFFFFFF : 0x000000000000FFFF;
//	uint64_t neutron_mask =
//		num_neutrons > 1 ? 0x0000FFFFFFFF0000 : 0x0000FFFF00000000;
	uint64_t particle_strids = 
		num_neutrons > 1 ?
		 0x0000000100010000 : 0x0000000100000000;
	for (size_t proton_state_index = 0;
	     proton_state_index < num_proton_states;
	     proton_state_index++)
	{
		short *current_proton_state = 
			(short*)(packed_proton_states +
				 num_protons*proton_state_index);
		uint64_t current_packed_state = 0;
		for (size_t i = 0; i < num_protons; i++)
			current_packed_state |= 
				(((uint64_t)
				  (current_proton_state[i])) <<
				(16*i)) << 1;
		for (size_t neutron_state_index = 0;
		     neutron_state_index < num_neutron_states;
		     neutron_state_index++)
		{
			short *current_neutron_state =
				(short*)(packed_neutron_states +
				       num_neutrons*neutron_state_index);
			current_packed_state &= proton_mask;	
			for (size_t i = 0; i < num_neutrons; i++)
				current_packed_state |=
					(((uint64_t)
					 (current_neutron_state[i])) <<
					(16*(i+num_protons)))<<1;
			current_packed_state |= particle_strids;
			log_entry("current_packed_state %lu: %lx",
				  state_index,
				  current_packed_state);
			basis->states[state_index].a = 
				(sp_state_index)(current_packed_state &
						 0x000000000000FFFF);
			basis->states[state_index].b = 
				(sp_state_index)((current_packed_state &
						  0x00000000FFFF0000) >> 16);
			basis->states[state_index].c = 
				(sp_state_index)((current_packed_state &
						  0x0000FFFF00000000) >> 32);
			sort_state(&basis->states[state_index]);
			state_index++;
		}
	}
	free(packed_proton_states);
	free(packed_neutron_states);
#pragma omp critical(sort_anicr_m_scheme_3p_basis)
	{
		current_sp_states = basis->sp_states;
		qsort(basis->states,
		      basis->dimension,
		      sizeof(M_Scheme_3p_State),
		      (__compar_fn_t)compare_m_scheme_3p_states_m);
		current_sp_states = NULL;
	}
}

static
void setup_state_to_index_map(M_Scheme_3p_Basis *basis)
{
	basis->state_to_index_map = 
		new_hash_map(basis->dimension,
			     2,
			     sizeof(size_t),
			     sizeof(M_Scheme_3p_State),
			     (__compar_fn_t)compare_m_scheme_3p_state,
			     (__hash_fn_t)compute_m_scheme_3p_hash);
	for (size_t i = 0; i < basis->dimension; i++)
	{
		hash_map_insert(basis->state_to_index_map,
				(void*)&basis->states[i],
				(void*)&i);
	}
}

static
int compare_m_scheme_3p_state(M_Scheme_3p_State *state_a,
			      M_Scheme_3p_State *state_b)
{
	int diff = state_a->a - state_b->a;
	if (diff)
		return diff;
	diff = state_a->b - state_b->b;
	if (diff)
		return diff;
	return state_a->c - state_b->c;
}

static
uint64_t compute_m_scheme_3p_hash(M_Scheme_3p_State *state)
{
	uint64_t h1 = (uint64_t)state->a ^ (((uint64_t)state->b)<<13);
	uint64_t h2 = (((uint64_t)state->b)<<3) ^ (((uint64_t)state->c)<<17);
	static const uint64_t magic_number = 0xFEDCBA9876543210;
	return h1 ^ __builtin_bswap64(h2) ^ magic_number;
}

static
void sort_state(M_Scheme_3p_State *state)
{
#define swap(a,b)\
	{\
		sp_state_index tmp = a;\
		a = b;\
		b = tmp;\
	}
	if (state->a > state->b)
		swap(state->a,state->b);
	if (state->b > state->c)
		swap(state->b,state->c);
	if (state->a > state->b)
		swap(state->a,state->b);
}

static
int compare_m_scheme_3p_states_m(M_Scheme_3p_State *state_a,
				 M_Scheme_3p_State *state_b)
{
	int M_a = 
		current_sp_states->sp_states[state_a->a].m +
		current_sp_states->sp_states[state_a->b].m +
		current_sp_states->sp_states[state_a->c].m; 
	int M_b = 
		current_sp_states->sp_states[state_b->a].m +
		current_sp_states->sp_states[state_b->b].m +
		current_sp_states->sp_states[state_b->c].m; 
	return M_a-M_b;
}

new_test(m_scheme_3p_nnp_basis_nmax_2,
	 const char *proton_basis_filename =
	 TEST_DATA
	 "anicr_basis_files/anicr_bases_protons_nmax2/p/basis_energy_1";
	 const char *neutron_basis_filename =
	 TEST_DATA
	 "anicr_basis_files/anicr_bases_neutrons_nmax2/nn/basis_energy_1";
	 Shells *shells = new_shells(2);
	 SP_States *sp_states = new_sp_states(shells);
	 M_Scheme_3p_Basis *basis =
	 new_m_scheme_3p_basis_from_basis_file(proton_basis_filename,
					       neutron_basis_filename,
					       1,2,
					       sp_states);
	list_m_scheme_3p_basis(basis);
       	free_m_scheme_3p_basis(basis);
	free_sp_states(sp_states);
	free_shells(shells);
	);	

new_test(m_scheme_3p_npp_basis_nmax_2,
	 const char *proton_basis_filename =
	 TEST_DATA
	 "anicr_basis_files/anicr_bases_protons_nmax2/pp/basis_energy_1";
	 const char *neutron_basis_filename =
	 TEST_DATA
	 "anicr_basis_files/anicr_bases_neutrons_nmax2/n/basis_energy_1";
	 Shells *shells = new_shells(2);
	 SP_States *sp_states = new_sp_states(shells);
	 M_Scheme_3p_Basis *basis =
	 new_m_scheme_3p_basis_from_basis_file(proton_basis_filename,
					       neutron_basis_filename,
					       2,1,
					       sp_states);
	list_m_scheme_3p_basis(basis);
       	free_m_scheme_3p_basis(basis);
	free_sp_states(sp_states);
	free_shells(shells);
	);	
