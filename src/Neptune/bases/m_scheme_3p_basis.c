#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <bases/m_scheme_3p_basis.h>
#include <string_tools/string_tools.h>
#include <utils/permutation_tools.h>


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
	return mp_basis;
}


void free_m_scheme_3p_basis(M_Scheme_3p_Basis* mp_basis)
{
	free(mp_basis->states);
	free(mp_basis);
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
	return mp_basis;
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
	{
		M_Scheme_3p_State current_state
			= mp_basis->states[known_element];
		int diff = diff_from_block(current_state,
					   tz,m_tot,E,
					   mp_basis->sp_states);
		if (diff !=0)
		{
			return 1;
		}
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
	*offset = start;
	return out;
}
