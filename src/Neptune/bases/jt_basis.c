#include <bases/jt_basis.h>
#include <stdint.h>
#include <array_builder/array_builder.h>
#include <utils/debug_messages.h>
#include <debug_mode/debug_mode.h>
#include <thundertester/test.h>

struct _jt_basis_
{
	jt_state_t *states;
	size_t dimension;
	Shells *shells;
	size_t *index_bins;
	size_t num_bins;
	int is_a_copy;
	size_t start_index;
};

typedef struct
{
	uint64_t parity:1;
	uint64_t T:1;
	uint64_t J:32;
} antoine_block_t;

typedef union
{
	size_t block_index;
	antoine_block_t block;
} antoine_block_iterator_t;

static
void append_block_basis(antoine_block_t block,
			Shells *shells,
			array_builder_t builder,
			quantum_number e_max2);

static
void fill_up_index_bins(jt_basis_t basis);

static
size_t find_state(jt_basis_t basis, jt_state_t state);

static
size_t bin_empty_index(jt_basis_t basis, jt_state_t state);

static inline
uint64_t compute_hash(jt_state_t state);

static inline
int compare_states(jt_state_t left_state, jt_state_t right_state);

static inline
quantum_number compute_state_parity(jt_state_t state,
				    Shells *shells);

static inline
int pauli_principle(jt_state_t state);

static inline
jt_state_t *find_block(jt_basis_t basis,
		       quantum_number J,
		       quantum_number T);

static inline
size_t find_block_length(jt_basis_t basis,
			 quantum_number J,
			 quantum_number T);

jt_basis_t new_antoine_basis(quantum_number e_max1,
			     quantum_number e_max2)
{
	jt_basis_t basis = (jt_basis_t)calloc(1,sizeof(struct _jt_basis_));
	basis->shells = new_antoine_shells(e_max1);
	array_builder_t builder = 
		new_array_builder((void**)&basis->states,
				  &basis->dimension,
				  sizeof(jt_state_t));
	const size_t num_blocks= 4*(e_max2+2);	
	for (antoine_block_iterator_t iterator = {0};
	     iterator.block_index < num_blocks;
	     iterator.block_index++)
	{
		DEBUG_MESS("block: %d %d %d\n",
			   iterator.block.J,
			   iterator.block.T,
			   iterator.block.parity);
		append_block_basis(iterator.block,
				   basis->shells,
				   builder,
				   e_max2);
	}
	free_array_builder(builder);
	fill_up_index_bins(basis);
	return basis;
}

jt_basis_t get_jt_block_basis(jt_basis_t origin,
			      quantum_number J,
			      quantum_number T)
{
	jt_basis_t basis = (jt_basis_t)calloc(1,sizeof(struct _jt_basis_));
	basis->is_a_copy = 1;
	basis->shells = origin->shells;
	basis->states = find_block(origin,J,T);
	if (basis->states == NULL)
	{
		free(basis);
		return NULL;
	}
	size_t start = basis->states - origin->states;
	basis->dimension = origin->dimension-start;
	basis->dimension = find_block_length(basis,J,T);
	if (basis->dimension == 0)
	{
		free(basis);
		return NULL;
	}
	fill_up_index_bins(basis);
	basis->start_index = start;
	for (size_t i = 0; i<basis->num_bins; i++)
	{
		if (basis->index_bins[i] != 0)
			basis->index_bins[i] += start;	
	}
	return basis;
}


index_hash_t compute_used_indices(jt_basis_t basis)
{
	const size_t max_num_bins = basis->dimension*basis->dimension;
	index_hash_t used_indices = new_index_hash(max_num_bins);
	size_t combined_index = 0;
	for (size_t left_index = 0; 
	     left_index < basis->dimension;
	     left_index++)
	{
		jt_state_t left_state = basis->states[left_index];
		quantum_number parity = 
			compute_state_parity(left_state,
					     basis->shells);
		for (size_t right_index = left_index;
		     right_index < basis->dimension &&
		     basis->states[right_index].J == left_state.J &&
		     basis->states[right_index].T == left_state.T;
		     right_index++)
		{
			jt_state_t right_state = basis->states[right_index];
			DEBUG_MESS("config(%lu) %lu %lu %lu %lu %d %d %d\n",
				   combined_index,
				   left_state.a,left_state.b,
				   right_state.a,right_state.b,
				   left_state.J,left_state.T,parity);
			if (parity != compute_state_parity(right_state,
							   basis->shells))
				break;
			DEBUG_MESS("config(%lu) %lu %lu %lu %lu %d %d %d\n",
				   combined_index,
				   left_state.a,left_state.b,
				   right_state.a,right_state.b,
				   left_state.J,left_state.T,parity);
			set_index(used_indices,
				  left_index,
				  right_index,
				  combined_index++);
		}
	}
	return used_indices;
}

size_t *get_sub_basis_indices(const jt_basis_t super_basis,
			      const jt_basis_t sub_basis)
{
	size_t *indices = (size_t*)malloc(sub_basis->dimension*sizeof(size_t));
	for (size_t i = 0; i<sub_basis->dimension; i++)
		indices[i] = find_state(super_basis,sub_basis->states[i]);	
	return indices;
}

size_t get_dimension(const jt_basis_t basis)
{
	return basis->dimension;
}

jt_state_t get_jt_state(const jt_basis_t basis,
			const size_t index)
{
	return basis->states[index];
}

size_t find_jt_state_index(jt_basis_t basis,
			   jt_state_t state)
{
	size_t index = find_state(basis,state);
	return index != no_index ? index-basis->start_index : no_index;	
}

void print_jt_basis(const jt_basis_t basis)
{
	printf("Printing jt_basis: %p\n",basis);
	for (size_t i = 0; i<basis->dimension; i++)
		printf("(%lu): (%lu %lu) %d %d\n",
		       i,
		       basis->states[i].a,
		       basis->states[i].b,
		       basis->states[i].J,
		       basis->states[i].T);
}

void print_used_indices(jt_basis_t basis)
{
	index_hash_t used_configurations = 
		compute_used_indices(basis);
	size_t k = 0;
	for (size_t i = 0; i<get_dimension(basis); i++)
		for (size_t j = 0; j<get_dimension(basis); j++)
			if ((k = get_index(used_configurations,
					   i,j)) != 
			    no_index)
				printf("(%lu) %2lu %2lu %2lu %2lu %2d %2d\n",
				       k,
				       basis->states[i].a+1,
				       basis->states[i].b+1,
				       basis->states[j].a+1,
				       basis->states[j].b+1,
				       basis->states[i].J,
				       basis->states[i].T);
	free_index_hash(used_configurations);
}

Shells *get_jt_basis_shells(jt_basis_t jt_basis)
{
	return jt_basis->shells;
}

void swap_shells(jt_state_t *state)
{
	true_shell_index temp_index = state->a;
	state->a = state->b;
	state->b = temp_index;
}

void free_jt_basis(jt_basis_t jt_basis)
{
	if (!jt_basis->is_a_copy)
	{
		free_shells(jt_basis->shells);
		free(jt_basis->states);
	}
	free(jt_basis->index_bins);
	free(jt_basis);
}

	static
void append_block_basis(antoine_block_t block,
			Shells *shells,
			array_builder_t builder,
			quantum_number e_max2)
{
	jt_state_t current_state =
	{
		.J = block.J,
		.T = block.T
	};
	for (current_state.a = 0;
	     current_state.a < shells->num_of_true_shells;
	     current_state.a++)
	{
		True_Shell shell_a = shells->true_shells[current_state.a];
		//DEBUG_MESS("shell_a.e = %d\n",shell_a.e);	
		if (2*shell_a.e > e_max2)
			break;
		for (current_state.b = current_state.a;
		     current_state.b < shells->num_of_true_shells;
		     current_state.b++)
		{
			True_Shell shell_b = 
				shells->true_shells[current_state.b];
			if (shell_a.e + shell_b.e > e_max2)
				break;
			if (!triangle_in_equality(shell_a,shell_b,block.J) ||
			    ((shell_a.l+shell_b.l)&1) != block.parity ||
			    !pauli_principle(current_state))
				continue;
			DEBUG_MESS("current_state = (%lu %lu) %d %d\n",
				   current_state.a,
				   current_state.b,
				   current_state.J,
				   current_state.T);
			append_array_element(builder,&current_state);
		}		
	}
}

	static
void fill_up_index_bins(jt_basis_t basis)
{
	basis->num_bins = basis->dimension*2;
	basis->index_bins = (size_t*)calloc(basis->num_bins,sizeof(size_t));
	for (size_t i = 0; i<basis->dimension; i++)
		basis->index_bins[bin_empty_index(basis,
						  basis->states[i])] = i+1;
}

	static
size_t find_state(jt_basis_t basis, jt_state_t state)
{
	uint64_t hash = compute_hash(state);
	size_t current_index = 0;
	while ((current_index = basis->index_bins[hash % basis->num_bins]) 
	       != 0 && 
	       (current_index == 0 || 
		compare_states(state,
			       basis->states[current_index-1-
			       basis->start_index]) != 0))
		hash++;
	return current_index-1;
}

	static
size_t bin_empty_index(jt_basis_t basis, jt_state_t state)
{
	uint64_t hash = compute_hash(state);
	while (basis->index_bins[hash % basis->num_bins] != 0)
		hash++;
	return hash % basis->num_bins;
}

	static inline
uint64_t compute_hash(jt_state_t state)
{
	uint64_t h1 = state.a ^ __builtin_bswap64(state.b);
	uint64_t h2 = state.J ^ __builtin_bswap64(state.J); 	
#define magic_number1 0x248C0BAF9023BB01
	uint64_t h3 = state.T*magic_number1;
#define magic_number2 0xF2890FAB54AA3CFF
	return h1^h2^h3^magic_number2;	
}

	static inline
int compare_states(jt_state_t left_state, jt_state_t right_state)
{
	int diff = left_state.J - right_state.J;
	if (diff)
		return diff;
	diff = left_state.T - right_state.T;
	if (diff)
		return diff;
	diff = left_state.a - right_state.a;
	if (diff)
		return diff;
	return left_state.b - right_state.b;
}

	static inline
quantum_number compute_state_parity(jt_state_t state,
				    Shells *shells)
{
	return (shells->true_shells[state.a].l +
		shells->true_shells[state.b].l)&1;
}

	static inline
int pauli_principle(jt_state_t state)
{
	return state.a != state.b || ((state.J+state.T)&1) == 1;
}

static inline
jt_state_t *find_block(jt_basis_t basis,
		       quantum_number J,
		       quantum_number T)
{
	for (size_t i = 0; i<basis->dimension; i++)
	{
		if (basis->states[i].J == J && basis->states[i].T == T)
			return basis->states+i;
	}
	return NULL;
}

static inline
size_t find_block_length(jt_basis_t basis,
			 quantum_number J,
			 quantum_number T)
{
	for (size_t i = 0; i<basis->dimension; i++)
		if (basis->states[i].J != J || basis->states[i].T != T)
			return i;
	return basis->dimension;
}

new_test(print_Nmax0_antoine_basis,
	 {
	 const quantum_number e_max1 = 0;
	 const quantum_number e_max2 = 0;
	 jt_basis_t basis = new_antoine_basis(e_max1,e_max2);
	 assert_that(basis != NULL);
	 print_jt_basis(basis);
	 free_jt_basis(basis);
	 });
new_test(print_Nmax2_antoine_basis,
	 {
	 const quantum_number e_max1 = 2;
	 const quantum_number e_max2 = 2;
	 jt_basis_t basis = new_antoine_basis(e_max1,e_max2);
	 assert_that(basis != NULL);
	 print_jt_basis(basis);
	 free_jt_basis(basis);
	 });

new_test(print_all_antoine_configuraitons,
	 {
	 const quantum_number e_max1 = 2;
	 const quantum_number e_max2 = 2;
	 Shells *shells = new_antoine_shells(e_max1);
	 list_shells(shells);	
	 free_shells(shells);
	 jt_basis_t basis = new_antoine_basis(e_max1,e_max2);
	 print_jt_basis(basis);
	 print_used_indices(basis);
	 free_jt_basis(basis);
	 });

new_test(self_match_test,
	 const quantum_number e_max1 = 2;
	 const quantum_number e_max2 = 2;
	 jt_basis_t basis = new_antoine_basis(e_max1,e_max2);
	 for (size_t i = 0; i < get_dimension(basis); i++)
	 {
	 printf("i_expect = %lu, i_found = %lu\n",
		i,find_jt_state_index(basis,basis->states[i]));
	 assert_that(i == find_jt_state_index(basis,basis->states[i]));
	 }
	);
