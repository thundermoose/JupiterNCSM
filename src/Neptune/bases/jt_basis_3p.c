#include <stdlib.h>
#include <string.h>
#include "jt_basis_3p.h"

JT_Basis* new_jt_basis(quantum_number e_max,
		quantum_number e_min,
		Shells *shells)
{
	// First determin minimu, shell index corresponding to e_min
	true_shell_index min_shell = find_minimum_true_shell(shells,e_min);

	// then maximum shell index corresponding to e_max
	true_shell_index max_shell = find_maximum_true_shell(shells,e_max);

	// We need to set up the structure to store the final basis
	// as standard we initiate all fields to 0
	JT_Basis* jt_basis = (JT_Basis*)calloc(1,sizeof(JT_Basis));
	jt_basis->shells = shells;

	// At this point in the program it is unknown how many
	// states that will be created, therefore
	// the code has to be open to the possibility that
	// the guessed number of basis states are wrong
	size_t allocated_num_states = 0; // Our initial guess

	// compute the limits of the tab loop
	quantum_number tab_min = 0;
	quantum_number tab_max = 2;

	// now we loop over the shells, with one loop for each particle
	JT_State current_state;
	for (current_state.a  = min_shell;
			current_state.a <= max_shell;
			current_state.a++)
	{
		// retrive often used quantum numbers
		quantum_number e_a = shells->true_shells[current_state.a].e;
		quantum_number j_a = shells->true_shells[current_state.a].j;
		for (current_state.b  = current_state.a;
				current_state.b <= max_shell; // we follow the convention set up by Gustav Jansen et.al
				current_state.b++)
		{
			// retrive often used quantum numbers
			quantum_number e_b = shells->true_shells[current_state.b].e;
			quantum_number j_b = shells->true_shells[current_state.b].j;
			// The program need to make sure that the total energy is not to large
			if (e_a+e_b>e_max)
			{
				break; // Assumes that the shells are ordered after energy
			}
			// Computing the limits for the jab loop
			quantum_number jab_min = abs(j_a-j_b);
			quantum_number jab_max = j_a+j_b;

			for (current_state.c  = min_shell;
					current_state.c <= max_shell;
					current_state.c++)
			{
				// retrive often used quantum numbers
				quantum_number e_c = shells->true_shells[current_state.c].e;
				quantum_number j_c = shells->true_shells[current_state.c].j;
				// The program need to make sure that the total energy is not to large
				if (e_a+e_b+e_c>e_max)
				{
					break; // Assumes that the shells are ordered after energy
				}
				for (current_state.jab  = jab_min;
						current_state.jab <= jab_max;
						current_state.jab += 2)
				{
					// compute the limit for the Jabc loop
					quantum_number jabc_min = abs(current_state.jab-j_c);
					quantum_number jabc_max = current_state.jab-j_c;
					for (current_state.jabc  = jabc_min;
							current_state.jabc <= jabc_max;
							current_state.jabc += 2)
					{
						for (current_state.tab  = tab_min;
								current_state.tab <= tab_max;
								current_state.tab += 2)
						{
							// compute the limits for tabc
							quantum_number tabc_min = abs(current_state.tab-1);
							quantum_number tabc_max = current_state.tab+1;

							for (current_state.tabc  = tabc_min;
									current_state.tabc <= tabc_max;
									current_state.tabc += 2)
							{

								// The program might have underestimated
								// the number of basis states
								if (allocated_num_states==jt_basis->dimension)
								{
									allocated_num_states+=allocated_num_states+1;
									jt_basis->states = (JT_State*)realloc(jt_basis->states,
											sizeof(JT_State)*
											allocated_num_states);
								}
								// append the current state to the end of the list
								jt_basis->states[jt_basis->dimension++] = current_state;
							}
						}
					}
				}
			}
		}
	}
	// the program might have overestimated
	// the number of basis states
	if (allocated_num_states>jt_basis->dimension)
	{
		allocated_num_states = jt_basis->dimension;
		jt_basis->states = (JT_State*)realloc(jt_basis->states,
				sizeof(JT_State)*
				allocated_num_states);
	}
	return jt_basis;
}

int jt_state_exists(JT_State state,
		JT_Basis *basis)
{
	size_t i;
	for (i = 0; i<basis->dimension; i++)
	{
		if (basis->states[i].a == state.a &&
				basis->states[i].b == state.b &&
				basis->states[i].c == state.c &&
				basis->states[i].jab == state.jab &&
				basis->states[i].jabc == state.jabc &&
				basis->states[i].tab == state.tab &&
				basis->states[i].tabc == state.tabc)
		{
			return 1;
		}
	}
	return 0;
}

JT_Basis* new_jt_basis_from_m_scheme(M_Scheme_3p_Basis* m_basis)
{
	// First the program needs a JT_Basis to store the result in
	// As we start of we want all feilds to be 0
	JT_Basis* jt_basis = (JT_Basis*)calloc(1,sizeof(JT_Basis));

	// we can set the shells field immediately
	jt_basis->shells = m_basis->sp_states->shells;

	// we keep a pointer to the sp_states, for later use
	SP_States* sp_states = m_basis->sp_states;

	// The program do not know how many JT_States that
	// that will be created, therefore we have to
	// keep track on how much memory we have stored
	// so far we have allocated none
	size_t allocated_number_states = 0;

	// The limits for the tab loop will be constant so we can
	// set them on forehand
	const quantum_number tab_min = 0;
	const quantum_number tab_max = 2;

	// Now for each m-scheme state in m_basis
	// we compute which JT_States that can exist
	size_t i;
	for (i = 0; i<m_basis->dimension; i++)
	{
		// We fill the JT_State with occupied
		// shells in the m-scheme state
		JT_State jt_state;
		jt_state.a =jt_basis->shells->shells[sp_states->sp_states[m_basis->states[i].a].shell].tse;
		jt_state.b =jt_basis->shells->shells[sp_states->sp_states[m_basis->states[i].b].shell].tse;
		jt_state.c =jt_basis->shells->shells[sp_states->sp_states[m_basis->states[i].c].shell].tse;

		// To simplify the code I put handles for the specific shells here
		True_Shell a = jt_basis->shells->true_shells[jt_state.a];
		True_Shell b = jt_basis->shells->true_shells[jt_state.b];
		True_Shell c = jt_basis->shells->true_shells[jt_state.c];

		// Computing the limits for the jab loop
		const quantum_number jab_min = abs(a.j-b.j);
		const quantum_number jab_max = a.j+b.j;
		for (jt_state.jab  = jab_min;
				jt_state.jab <= jab_max;
				jt_state.jab += 2)
		{
			// compute the limits for the jabc loop
			const quantum_number jabc_min = abs(jt_state.jab-c.j);
			const quantum_number jabc_max = jt_state.jab+c.j;
			for (jt_state.jabc  = jabc_min;
					jt_state.jabc <= jabc_max;
					jt_state.jabc += 2)
			{
				for (jt_state.tab  = tab_min;
						jt_state.tab <= tab_max;
						jt_state.tab += 2)
				{
					/*
					   if (jt_state.a == jt_state.b &&
					   ((jt_state.jab + jt_state.tab)/2)%2==0)
					   continue;
					   */
					// compute the limits for the tabc loop
					const quantum_number tabc_min = abs(jt_state.tab-1);
					const quantum_number tabc_max = abs(jt_state.tab+1);
					for (jt_state.tabc  = tabc_min;
							jt_state.tabc <= tabc_max;
							jt_state.tabc += 2)
					{
						// It is possible that the state has already been created
						// if that is so we do not need to create it again
						if (jt_state_exists(jt_state,
									jt_basis))
						{
							continue;
						}
						// If the program has underestimated
						// the number of final states, we need
						// to allocate more memory
						if (allocated_number_states == jt_basis->dimension)
						{
							allocated_number_states+=allocated_number_states*2+1;
							jt_basis->states = (JT_State*)realloc(jt_basis->states,
									sizeof(JT_State)*
									allocated_number_states);
						}
						// Appending the state to the end of the array
						jt_basis->states[jt_basis->dimension++] = jt_state;

					}
				}
			}
		}
	}
	// At this point the program might have overestimated
	// the number of final states, so we trim down the array
	if (allocated_number_states > jt_basis->dimension)
	{
		allocated_number_states=jt_basis->dimension;
		jt_basis->states = (JT_State*)realloc(jt_basis->states,
				sizeof(JT_State)*
				allocated_number_states);
	}
	// now we are done
	return jt_basis;
}

JT_Basis* new_jt_basis_from_ascii_file(FILE* file,
		Shells *shells)
{
	// we are going to read the file row by row
	// using getline, therefore we need to
	// initiate getline's arguments
	char* row = NULL; // getline will allocate the necessary memory
	size_t row_len = 0;

	// skip ahead to the three particle states
	do
	{
		if (getline(&row,&row_len,file) < 0)
		{
			fprintf(stderr,"No three particle states found\n");
			exit(1);
		}
	}
	while (row != NULL &&
			strcmp(row,"---three-particle-states---\n") != 0);
	// Making sure that the format is correct
	if (getline(&row,&row_len,file) < 0)
	{
		fprintf(stderr,"File ended directly after three particle states heading\n");
		exit(1);
	}

	if (strcmp(row,"---i---a---b---c---Jab-Jabc-Tab-Tabc---\n") != 0)
	{
		fprintf(stderr,"Expected JT format\n");
		exit(1);
	}
	// Now we can start constructing the JT_Basis
	JT_Basis* jt_basis = (JT_Basis*)calloc(1,sizeof(JT_Basis));
	jt_basis->shells = shells;

	// We do not know how many states that are in the file
	// therefore we do the very qualified guess of 0
	size_t allocated_num_states = 0;

	do
	{
		if (getline(&row,&row_len,file) < 0)
		{
			// file ended which is ok
			break;
		}
		// trying to read the row as a JT_State,
		// the first index is ignored (this could be a problem)
		JT_State current;
		if (sscanf(row,"  %*d  %ld  %ld  %ld  %d  %d  %d  %d",
					&current.a,
					&current.b,
					&current.c,
					&current.jab,
					&current.jabc,
					&current.tab,
					&current.tabc) == 7)
		{


			// the program might have underestimated the number
			// states
			if (allocated_num_states == jt_basis->dimension)
			{
				allocated_num_states += allocated_num_states+1;
				jt_basis->states =
					(JT_State*)realloc(jt_basis->states,
							sizeof(JT_State)*
							allocated_num_states);
			}
			// append the current state to the end of the array
			jt_basis->states[jt_basis->dimension++] = current;
		}
	}
	while (strcmp(row,"\n") != 0);
	// The program might have overestimated the necessary size
	// of the array
	if (allocated_num_states > jt_basis->dimension)
	{
		allocated_num_states = jt_basis->dimension;
		jt_basis->states =
			(JT_State*)realloc(jt_basis->states,
					sizeof(JT_State)*
					allocated_num_states);
	}
	free(row);
	return jt_basis;
}


JT_Basis* get_jt_block(JT_Basis* jt_basis,
		quantum_number jabc,
		quantum_number tabc)
{
	// We need somewhere to store the results
	// we do also initiate all fields to 0
	JT_Basis* jt_sub_basis = (JT_Basis*)calloc(1,sizeof(JT_Basis));
	jt_sub_basis->shells = jt_basis->shells;

	// We do not know how many states with j_abc and t_abc
	// but it cannot be more than jt_basis->dimension
	size_t allocated_num_states = jt_basis->dimension;
	jt_sub_basis->states = (JT_State*)malloc(sizeof(JT_State)*allocated_num_states);

	size_t i;
	for (i = 0; i<jt_basis->dimension; i++)
	{
		// checking if the current state corresponds to
		// our criteria
		if (jt_basis->states[i].jabc == jabc &&
				jt_basis->states[i].tabc == tabc)
		{
			jt_sub_basis->states[jt_sub_basis->dimension++] = jt_basis->states[i];
		}
	}
	// if the program overestimated the number of states
	// we need to trim the array
	if (jt_sub_basis->dimension<allocated_num_states)
	{
		allocated_num_states = jt_sub_basis->dimension;
		jt_sub_basis->states = (JT_State*)realloc(jt_sub_basis->states,
				sizeof(JT_State)*
				allocated_num_states);
	}
	// now we are done
	return jt_sub_basis;
}

ssize_t find_trans_jt_state(JT_Basis* look_in,
		JT_State look_for,
		ssize_t* for_to_in_trans)
{
	if (for_to_in_trans[look_for.a]<0)
		return -1;
	look_for.a = for_to_in_trans[look_for.a];
	if (for_to_in_trans[look_for.b]<0)
		return -1;
	look_for.a = for_to_in_trans[look_for.b];
	if (for_to_in_trans[look_for.c]<0)
		return -1;
	look_for.a = for_to_in_trans[look_for.c];
	size_t i;
	for (i = 0; i<look_in->dimension; i++)
	{
		// This is what determines the order
		if (look_in->states[i].a == look_for.a &&
				look_in->states[i].b == look_for.b &&
				look_in->states[i].c == look_for.c &&
				look_in->states[i].jab == look_for.jab &&
				look_in->states[i].jabc == look_for.jabc &&
				look_in->states[i].tab == look_for.tab &&
				look_in->states[i].tabc == look_for.tabc)
		{
			return i;
		}
	}
	return -1;
}


void print_jt_basis_3p(JT_Basis *jt_basis)
{
	size_t i;
	for (i = 0; i<jt_basis->dimension; i++)
	{
		printf("(%ld): ((%ld %ld) %d %d, %ld) %d %d)\n",
				i,
				jt_basis->states[i].a,
				jt_basis->states[i].b,
				jt_basis->states[i].jab,
				jt_basis->states[i].tab,
				jt_basis->states[i].c,
				jt_basis->states[i].jabc,
				jt_basis->states[i].tabc);
	}
}


int contains_duplicates(JT_Basis *jt_basis)
{
	size_t i,j;
	for (i = 0; i<jt_basis->dimension-1; i++)
	{
		for (j = i+1; j<jt_basis->dimension; j++)
		{
#define QN(k,qn) jt_basis->states[k].qn
			if (QN(i,a) == QN(j,a) &&
					QN(i,b) == QN(j,b) &&
					QN(i,c) == QN(j,c) &&
					QN(i,jab) == QN(j,jab) &&
					QN(i,jabc) == QN(j,jabc) &&
					QN(i,tab) == QN(j,tab) &&
					QN(i,tabc) == QN(j,tabc))
			{
				return 1;
			}
		}
	}
	return 0;
}


ssize_t* matching_jt_states(JT_Basis* jt_basis_a,
		JT_Basis* jt_basis_b,
		ssize_t* a_to_b_shells)
{
	ssize_t* corresponding =
		(ssize_t*)malloc(sizeof(ssize_t)*jt_basis_a->dimension);

	size_t i;
	for (i = 0;
			i < jt_basis_a->dimension;
			i++)
	{
		corresponding[i] =
			find_trans_jt_state(jt_basis_b,
					jt_basis_a->states[i],
					a_to_b_shells);
	}
	return corresponding;
}


void free_jt_basis_3p(JT_Basis* jt_basis)
{
	free(jt_basis->states);
	free(jt_basis);
}
