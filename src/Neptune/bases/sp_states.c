#include "sp_states.h"
#include <debug_mode/debug_mode.h>
#include <unit_testing/test.h>

Shells* _shells;

int comp_sp_states(SP_State* a,
		   SP_State* b)
{
	int diff =
		_shells->shells[a->shell].e-
		_shells->shells[b->shell].e;
	if (diff)
		return diff;
	diff =
		_shells->shells[a->shell].l-
		_shells->shells[b->shell].l;
	if (diff)
		return diff;
	diff =
		_shells->shells[a->shell].j-
		_shells->shells[b->shell].j;
	if (diff)
		return diff;
	diff =
		b->m-
		a->m;
	if (diff)
		return diff;
	diff =
		_shells->shells[a->shell].tz-
		_shells->shells[b->shell].tz;
	return diff;

}

/* Given a list of all shells
 * this generates a list of all
 * possible single particle states
 * hence it adds the m quantum number
 * (a.k.a jz)
 */
SP_States* new_sp_states(Shells* shells){
	SP_States* sp_states = (SP_States*)malloc(sizeof(SP_States));
	sp_states->shells = shells;
	// Determin dimension of the Hilbert space
	shell_index i;
	sp_state_index dim=0;
	for (i = 0; i<shells->num_of_shells; i++){
		dim+=(shells->shells[i].j+1);
	}
	// Allocate memory for the array
	sp_states->dimension = dim;
	sp_states->sp_states = (SP_State*)malloc(sizeof(SP_States)*dim);
	// Fill the array
	sp_state_index j = 0;
	for (i = 0; i<shells->num_of_shells; i++){
		SP_State s;
		s.shell = i;
		for (s.m =-shells->shells[i].j;
		     s.m<=shells->shells[i].j;
		     s.m+=2){

			sp_states->sp_states[j++]=s;

		}
	}
	_shells = shells;
	qsort(sp_states->sp_states,
	      sp_states->dimension,
	      sizeof(SP_State),
	      (__compar_fn_t)comp_sp_states);

	return sp_states;
}

void list_sp_states(SP_States* s)
{
	list_shells(s->shells);
	printf("sp_states:\n");
	size_t i;
	for (i = 0; i<s->dimension; i++)
	{
		Shell shell = s->shells->shells[s->sp_states[i].shell];
		printf("(%ld): %ld %d: (%d %d %d %d %d)\n",
		       i,
		       s->sp_states[i].shell,
		       s->sp_states[i].m,
		       shell.n,
		       shell.l,
		       shell.j,
		       s->sp_states[i].m,
		       shell.tz
		      );
	}
}


/* This frees all memory allocated by the function
 * above
 */
void free_sp_states(SP_States* sp_states){
	free(sp_states->sp_states);
	free(sp_states);
}

new_test(print_nmax2_antoine_basis,
	 Shells *shells = new_antoine_shells(2);
	 SP_States *states = new_sp_states(shells);
	 list_sp_states(states);
	 free_sp_states(states);
	 free_shells(shells);
	);
