#include <single_particle_basis/single_particle_basis.h>
#include <array_builder/array_builder.h>
#include <debug_mode/debug_mode.h>
#include <unit_testing/test.h>
#include <string.h>

struct _single_particle_basis_
{
	single_particle_state_t *states;	
	size_t num_states;
};

int compare_neptune_states(const single_particle_state_t *state_a,
			   const single_particle_state_t *state_b);

single_particle_basis_t new_single_particle_basis(int energy_max)
{
	single_particle_basis_t basis =
		(single_particle_basis_t)
		calloc(1,sizeof(struct _single_particle_basis_));
	array_builder_t basis_builder =
		new_array_builder((void**)&basis->states,
				  &basis->num_states,
				  sizeof(single_particle_state_t));
	single_particle_state_t current_state;
	current_state.neptune_index = 0;
	for (int energy = 0; energy <= energy_max; energy++)
	{
		for (current_state.l = energy % 2;
		     current_state.l <= energy;
		     current_state.l += 2)
		{
			current_state.n = (energy-current_state.l)/2;
			for (current_state.j = abs(current_state.l*2-1);
			     current_state.j <= current_state.l*2+1;
			     current_state.j += 2)
			{

				for (current_state.m = current_state.j;
				     current_state.m >= -current_state.j;
				     current_state.m -= 2)
				{
					current_state.tz = -1;
					append_array_element(basis_builder,
							     &current_state);
					current_state.neptune_index++;
					current_state.tz = 1;
					append_array_element(basis_builder,
							     &current_state);
					current_state.neptune_index++;
				}
			}
		}
	}
	free_array_builder(basis_builder);
	return basis;
}

single_particle_state_t get_state(single_particle_basis_t basis,
				  size_t index)
{
	return basis->states[index];
}

void free_single_particle_basis(single_particle_basis_t basis)
{
	free(basis->states);
	free(basis);
}

int get_m(single_particle_state_t state)
{
	return state.m;
}

int compare_neptune_states(const single_particle_state_t *state_a,
			   const single_particle_state_t *state_b)
{
	int diff;
	int Na = state_a->n*2+state_a->l;
	int Nb = state_b->n*2+state_b->l;
	diff = Na-Nb;
	if (diff)
		return diff;
#define compare(sign,quantum_number)\
	diff = sign*(state_a->quantum_number-state_b->quantum_number);\
	if (diff)\
		return diff;
	compare(1,l);
	compare(-1,j);
	compare(-1,m);
	compare(1,tz);
	return 0;
}

new_test(nmax2_single_particle_basis,
	 {
	 single_particle_basis_t basis = new_single_particle_basis(2);
	 for (size_t i = 0; i<basis->num_states; i++)
	 {
		printf("(%lu) %d %d %d %d %d (%lu)\n",
		       i,
		       basis->states[i].n,
		       basis->states[i].l,
		       basis->states[i].j,
		       basis->states[i].m,
		       basis->states[i].tz,
		       basis->states[i].neptune_index);
	 }
	 free_single_particle_basis(basis);
	 });
