#ifndef __JT_BASIS__
#define __JT_BASIS__

#include <bases/shells.h>
#include <utils/index_hash.h>
struct _jt_basis_;
typedef struct _jt_basis_ *jt_basis_t;

typedef struct
{
	true_shell_index a,b;
	quantum_number J,T;
} jt_state_t;

jt_basis_t new_antoine_basis(quantum_number e_max1,
			     quantum_number e_max2);

jt_basis_t get_jt_block_basis(jt_basis_t origin,
			      quantum_number J,
			      quantum_number T);

index_hash_t compute_used_indices(jt_basis_t basis);

size_t *get_sub_basis_indices(const jt_basis_t super_basis,
			      const jt_basis_t sub_basis);
size_t get_dimension(const jt_basis_t basis);

jt_state_t get_jt_state(const jt_basis_t basis,
			size_t index);

size_t find_jt_state_index(jt_basis_t basis,
			   jt_state_t state);

void print_jt_basis(const jt_basis_t basis);

void print_used_indices(jt_basis_t basis);

Shells *get_jt_basis_shells(jt_basis_t basis);

void swap_shells(jt_state_t *state);

void free_jt_basis(jt_basis_t jt_basis);
#endif
