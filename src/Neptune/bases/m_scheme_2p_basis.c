#include "m_scheme_2p_basis.h"
#include <array_builder/array_builder.h>
#include <utils/helpful_macros.h>
#include <error/error.h>
#include <thundertester/test.h>
#include <assert.h>
#include <utils/debug_messages.h>
#include <errno.h>
#include <string.h>

struct _m_scheme_2p_basis_
{
	Shells *shells;
	SP_States *sp_states;
	size_t *corresponding_indices;
	m_scheme_2p_state_t *states;
	size_t dimension;
	quantum_number e_max1;
	quantum_number e_max2;
};

typedef struct
{
	quantum_number Tz;
	Shells *shells;
} fixed_tz_data_t;

typedef int (*state_condition_t)(SP_State state_a,
				 SP_State state_b,
				 void *data);
static 
int any_thing_goes(SP_State state_a,
		   SP_State state_b,
		   void *data);
static
int total_m_is_zero(SP_State state_a,
		    SP_State state_b,
		    void *data);
static
int total_tz_is_fixed(SP_State state_a,
		      SP_State state_b,
		      fixed_tz_data_t *data);

static
quantum_number total_tz(m_scheme_2p_state_t current_state,
			SP_States *sp_states);

static
quantum_number total_m(m_scheme_2p_state_t current_state,
		       SP_States *sp_states);
static
quantum_number total_energy(m_scheme_2p_state_t current_state,
			    SP_States *sp_states);

static
void read_two_particle_type_files(m_scheme_2p_basis_t basis,
				  const char *proton_basis_filename,
				  const char *neutron_basis_filename);

static
void read_single_particle_type_file(m_scheme_2p_basis_t basis,
				    const char *basis_filename,
				    const int iso_spin);


static
void read_states(void **states,
		 size_t *num_states,
		 size_t state_size,
		 const char *filename);

m_scheme_2p_basis_t new_m_scheme_2p_basis_condition(quantum_number e_max1,
						    quantum_number e_max2,
						    Shells *shells,
						    state_condition_t condition,
						    void *data)
{
	m_scheme_2p_basis_t m_scheme_2p_basis =
		(m_scheme_2p_basis_t)calloc(1,sizeof(struct _m_scheme_2p_basis_));
	m_scheme_2p_basis->e_max1 = e_max1;
	m_scheme_2p_basis->e_max2 = e_max2;
	m_scheme_2p_basis->shells = shells;
	SP_States *sp_states =
		m_scheme_2p_basis->sp_states =
		new_sp_states(m_scheme_2p_basis->shells);
	m_scheme_2p_basis->corresponding_indices = NULL;
	array_builder_t states_builder = 
		new_array_builder(
				  (void*)&m_scheme_2p_basis->states,
				  &m_scheme_2p_basis->dimension,
				  sizeof(m_scheme_2p_state_t));
	m_scheme_2p_state_t current_state = {0};
	for (current_state.a = 0;
	     current_state.a<sp_states->dimension-1;
	     current_state.a++)
	{
		SP_State state_a =
			sp_states->sp_states[current_state.a];
		Shell shell_a = 
			shells->shells[state_a.shell];
		if (shell_a.e * 2 >e_max2)
			break;
		for (current_state.b = current_state.a+1;
		     current_state.b < sp_states->dimension;
		     current_state.b++)
		{
			SP_State state_b =
				sp_states->sp_states[current_state.b];
			Shell shell_b = 
				shells->shells[state_b.shell];
			if (shell_a.e+shell_b.e > e_max2)
				break;
			if (!condition(state_a,state_b,data))
				continue;
			append_array_element(states_builder,
				       &current_state);
		}
	}
	free_array_builder(states_builder);
	return m_scheme_2p_basis;
}

m_scheme_2p_basis_t new_m_scheme_2p_basis(quantum_number e_max1,
					  quantum_number e_max2)
{
	Shells *shells = new_shells(e_max1);
	return new_m_scheme_2p_basis_condition(e_max1,e_max2,
					       shells,
					       any_thing_goes,NULL);
}

m_scheme_2p_basis_t new_m_scheme_2p_basis_fixed_isospin(quantum_number e_max1,
							quantum_number e_max2,
							quantum_number Tz)
{
	fixed_tz_data_t data =
	{
		.Tz = Tz,
		.shells = new_shells(e_max1)
	};
	return new_m_scheme_2p_basis_condition(e_max1,e_max2,
					       data.shells,
					       (state_condition_t)
					       total_tz_is_fixed,&data);
}

m_scheme_2p_basis_t 
new_m_scheme_2p_basis_from_files(quantum_number e_max1,
				 const char *proton_basis_filename,
				 const char *neutron_basis_filename)
{
	m_scheme_2p_basis_t basis =
	       	(m_scheme_2p_basis_t)
		malloc(sizeof(struct _m_scheme_2p_basis_));	
	basis->shells = new_antoine_shells(e_max1);
	basis->sp_states = new_sp_states(basis->shells);
	basis->corresponding_indices = NULL;
	quantum_number iso_spin = 
		(neutron_basis_filename != NULL) - 
		(proton_basis_filename != NULL);
	if (iso_spin == 0)
	{
		read_two_particle_type_files(basis,
					     proton_basis_filename,
					     neutron_basis_filename);
	}	
	else
	{
		const char *basis_filename =
			iso_spin < 0 ? 
			proton_basis_filename :
			neutron_basis_filename;
		read_single_particle_type_file(basis,
					       basis_filename,
					       iso_spin);
	}
	return basis;
}

m_scheme_2p_basis_t generate_2p_block(m_scheme_2p_basis_t m_scheme_2p_basis,
				      quantum_number Tz,
				      quantum_number M,
				      quantum_number energy)
{
	m_scheme_2p_basis_t basis =
		(m_scheme_2p_basis_t)
		calloc(1,sizeof(struct _m_scheme_2p_basis_));
	basis->shells = m_scheme_2p_basis->shells;
	basis->sp_states = m_scheme_2p_basis->sp_states;
	basis->e_max1 = m_scheme_2p_basis->e_max1;
	basis->e_max2 = energy;
	array_builder_t basis_builder =
		new_array_builder((void**)&basis->states,
				  &basis->dimension,
				  sizeof(m_scheme_2p_state_t));
	size_t num_correspondin_indices = 0;
	array_builder_t index_builder =
		new_array_builder((void**)&basis->corresponding_indices,
				  &num_correspondin_indices,
				  sizeof(size_t)); 
	for (size_t i = 0; i < m_scheme_2p_basis->dimension; i++)
	{
		m_scheme_2p_state_t current_state =
			m_scheme_2p_basis->states[i];	
		if (total_tz(current_state,
			     m_scheme_2p_basis->sp_states) != Tz ||
		    total_m(current_state,
			    m_scheme_2p_basis->sp_states) != M ||
		    total_energy(current_state,
				 m_scheme_2p_basis->sp_states) != energy)
			continue;
		DEBUG_MESS("add state: %lu, (%d %d)\n",i,
			current_state.a,
			current_state.b);
		append_array_element(basis_builder,&current_state);
		append_array_element(index_builder,&i);
	}
	assert(num_correspondin_indices == basis->dimension);
	free_array_builder(basis_builder);
	free_array_builder(index_builder);
	if (basis->dimension == 0)
	{
		free(basis);
		return NULL;
	}
	return basis;
}

size_t get_m_scheme_2p_dimension(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	return m_scheme_2p_basis->dimension;
}

quantum_number get_m_scheme_2p_e_max1(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	return m_scheme_2p_basis->e_max1;
}

quantum_number get_m_scheme_2p_e_max2(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	return m_scheme_2p_basis->e_max2;
}

void print_m_scheme_2p_basis(m_scheme_2p_basis_t basis)
{
	SP_States *sp_states = basis->sp_states;
	list_sp_states(sp_states);
	printf("m-scheme basis:\n");
	Shells *shells = sp_states->shells;
	if (basis->corresponding_indices == NULL)
		for (size_t i = 0; i<basis->dimension; i++)
		{
			sp_state_index state_a = basis->states[i].a;
			sp_state_index state_b = basis->states[i].b;
			shell_index shell_a = 
				sp_states->sp_states[state_a].shell;
			true_shell_index tshell_a =
			       	shells->shells[shell_a].tse;
			shell_index shell_b = 
				sp_states->sp_states[state_b].shell;
			true_shell_index tshell_b =
			       	shells->shells[shell_b].tse;
			printf("(%lu): (%d:%lu:%lu) (%d:%lu:%lu)\n",
					i,
					basis->states[i].a,
					shell_a,
					tshell_a,
					basis->states[i].b,
					shell_b,
					tshell_b);
		}
	else
		for (size_t i = 0; i<basis->dimension; i++)
			printf("(%lu:%lu): %d %d\n",
					i,basis->corresponding_indices[i],
					basis->states[i].a,
					basis->states[i].b);
}

Shells *get_m_scheme_shells(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	return m_scheme_2p_basis->sp_states->shells;
}

SP_States *get_m_scheme_sp_states(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	return m_scheme_2p_basis->sp_states;
}

m_scheme_2p_state_t get_m_scheme_2p_state(m_scheme_2p_basis_t m_scheme_2p_basis,
					  size_t index)
{
	return m_scheme_2p_basis->states[index];
}

m_scheme_2p_state_t *get_m_scheme_2p_states(m_scheme_2p_basis_t basis)
{
	m_scheme_2p_state_t *states =
		(m_scheme_2p_state_t*)malloc(basis->dimension*
					     sizeof(m_scheme_2p_state_t));
	memcpy(states,
	       basis->states,
	       basis->dimension*sizeof(m_scheme_2p_state_t));
	return states;
}

size_t *m_scheme_2p_corresponding_indices(m_scheme_2p_basis_t basis)
{
	return basis->corresponding_indices;
}


void free_m_scheme_2p_basis(m_scheme_2p_basis_t m_scheme_2p_basis)
{
	if (m_scheme_2p_basis->corresponding_indices == NULL)
	{
		free_sp_states(m_scheme_2p_basis->sp_states);
		free_shells(m_scheme_2p_basis->shells);
	}
	else
	{
		free(m_scheme_2p_basis->corresponding_indices);
	}
	free(m_scheme_2p_basis->states);
	free(m_scheme_2p_basis);
}

static 
int any_thing_goes(SP_State state_a,
		   SP_State state_b,
		   void *data)
{
	return 1;
}

	static
int total_m_is_zero(SP_State state_a, SP_State state_b,void *data)
{
	intentionaly_unused(data);
	return state_a.m == -state_b.m;
}

	static
int total_tz_is_fixed(SP_State state_a,
		      SP_State state_b,
		      fixed_tz_data_t *data)
{
	return total_m_is_zero(state_a,state_b,NULL) &&
		data->shells->shells[state_a.shell].tz +
		data->shells->shells[state_b.shell].tz ==
		2*data->Tz;
}

static
quantum_number total_tz(m_scheme_2p_state_t current_state,
			SP_States *sp_states)
{
	Shells *shells = sp_states->shells;
	return shells->shells[sp_states->sp_states[current_state.a].shell].tz + 
		shells->shells[sp_states->sp_states[current_state.b].shell].tz; 
}

static
quantum_number total_m(m_scheme_2p_state_t current_state,
		       SP_States *sp_states)
{
	return sp_states->sp_states[current_state.a].m + 
	       sp_states->sp_states[current_state.b].m; 
}
static
quantum_number total_energy(m_scheme_2p_state_t current_state,
			    SP_States *sp_states)
{
	Shells *shells = sp_states->shells;
	return shells->shells[sp_states->sp_states[current_state.a].shell].e + 
		shells->shells[sp_states->sp_states[current_state.b].shell].e; 
}

static
void read_single_particle_type_file(m_scheme_2p_basis_t basis,
				    const char *basis_filename,
				    const int iso_spin)
{
	unsigned int *packed_states = NULL;
	size_t num_states = 0;
	read_states((void**)&packed_states,
		    &num_states,
		    sizeof(unsigned int),
		    basis_filename);
	basis->dimension = num_states;
	basis->states =
		(m_scheme_2p_state_t*)
		malloc(num_states*sizeof(m_scheme_2p_state_t));
	const unsigned int particle_stride = 
		iso_spin > 0 ? 0x00010001 : 0;
	for (size_t i = 0; i < num_states; i++)
	{
		unsigned int packed_state =
			(packed_states[i]<<1) | particle_stride;
		basis->states[i].a = (packed_state & 0xFFFF0000) >> 16;
		basis->states[i].b = (packed_state & 0x0000FFFF);
	}
	free(packed_states);
}

static
void read_two_particle_type_files(m_scheme_2p_basis_t basis,
				  const char *proton_basis_filename,
				  const char *neutron_basis_filename)
{
	short *proton_packed_states = NULL;
	size_t num_proton_states = 0;
	read_states((void**)&proton_packed_states,
		    &num_proton_states,
		    sizeof(short),
		    proton_basis_filename);
	short *neutron_packed_states = NULL;
	size_t num_neutron_states = 0;
	read_states((void**)&neutron_packed_states,
		    &num_neutron_states,
		    sizeof(short),
		    neutron_basis_filename);

	basis->dimension = num_proton_states*num_neutron_states;
	basis->states =
	       	(m_scheme_2p_state_t*)
		malloc(basis->dimension*sizeof(m_scheme_2p_state_t));
	size_t state_index = 0;
	for (size_t proton_index = 0;
	     proton_index < num_proton_states;
	     proton_index++)
	{
		sp_state_index proton_state = 
			(sp_state_index)
			(proton_packed_states[proton_index])*2;
		for (size_t neutron_index = 0;
		     neutron_index < num_neutron_states;
		     neutron_index++)
		{
			sp_state_index neutron_state = 
				(sp_state_index)
				(neutron_packed_states[neutron_index])*2+1;

			m_scheme_2p_state_t state =
			{
				.a = min(proton_state,neutron_state),
				.b = max(proton_state,neutron_state)
			};
			basis->states[state_index++] = state;
		}
	}
	free(proton_packed_states);
	free(neutron_packed_states);
}

static
void read_states(void **states,
		 size_t *num_states,
		 size_t state_size,
		 const char *filename)
{
	FILE *basis_file = fopen(filename,"r");
	if (basis_file == NULL)
		error("Could not open file %s. %s\n",
		      filename,
		      strerror(errno));
	fseek(basis_file,0,SEEK_END);
	const size_t num_bytes = ftell(basis_file);
	fseek(basis_file,0,SEEK_SET);
	if (num_bytes % state_size != 0)
		error("Basis file %s do not include %lu bytes states.\n",
		      filename,
		      state_size);
	*num_states = num_bytes/state_size;
	*states = (unsigned int*)malloc(num_bytes);
	if (fread(*states,
		  state_size,
		  *num_states,
		  basis_file) != (*num_states))
		error("Could not read %lu bytes from %s.\n",
		      num_bytes,filename);
	fclose(basis_file);	
}

new_test(print_m_scheme_2p_basis_nmax2,
	 m_scheme_2p_basis_t basis = new_m_scheme_2p_basis(2,2);
	 assert_that(basis != NULL);
	 print_m_scheme_2p_basis(basis);
	 free_m_scheme_2p_basis(basis);
	);
new_test(print_m_scheme_2p_basis_nmax2_fixed_Tz,
	 m_scheme_2p_basis_t basis = 
	 new_m_scheme_2p_basis_fixed_isospin(2,2,0);
	 assert_that(basis != NULL);
	 print_m_scheme_2p_basis(basis);
	 free_m_scheme_2p_basis(basis);
	);

new_test(generate_2p_block_trivial,
	m_scheme_2p_basis_t basis = new_m_scheme_2p_basis(2,2);
	m_scheme_2p_basis_t block_basis = generate_2p_block(basis,
							2,0,2);
	print_m_scheme_2p_basis(block_basis);
	free_m_scheme_2p_basis(block_basis);
	free_m_scheme_2p_basis(basis);
);
