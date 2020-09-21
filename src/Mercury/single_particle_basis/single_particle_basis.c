#include <single_particle_basis/single_particle_basis.h>
#include <array_builder/array_builder.h>
#include <debug_mode/debug_mode.h>
#include <unit_testing/test.h>
#include <string.h>

struct _single_particle_basis_
{
	SP_States *sp_states;
	Shells *shells;
};

single_particle_basis_t new_single_particle_basis(int energy_max)
{
	single_particle_basis_t basis =
		(single_particle_basis_t)
		malloc(sizeof(struct _single_particle_basis_));	
	basis->shells = new_shells(energy_max);
	basis->sp_states = new_sp_states(basis->shells);
	return basis;
}

single_particle_state_t get_state(single_particle_basis_t basis,
				  size_t index)
{
	SP_State sp_state = basis->sp_states->sp_states[index];
	Shell shell = basis->shells->shells[sp_state.shell];
	single_particle_state_t state =
	{
		.n = shell.n,
		.l = shell.l,
		.j = shell.j,
		.m = sp_state.m,
		.tz = shell.tz,
		.neptune_index = index
	};
	return state;
}

SP_States *get_sp_states(single_particle_basis_t basis)
{
	return basis->sp_states;
}

void free_single_particle_basis(single_particle_basis_t basis)
{
	free_sp_states(basis->sp_states);
	free_shells(basis->shells);
	free(basis);
}

int get_m(single_particle_state_t state)
{
	return state.m;
}

new_test(nmax2_single_particle_basis,
	 {
	 single_particle_basis_t basis = new_single_particle_basis(2);
	 for (size_t i = 0; i<basis->num_states; i++)
	 {
	 	single_particle_state_t state = get_state(basis,i);
		printf("(%lu) %d %d %d %d %d (%lu)\n",
		       i,
		       state.n,
		       state.l,
		       state.j,
		       state.m,
		       state.tz,
		       state.neptune_index);
	 }
	 free_single_particle_basis(basis);
	 });
