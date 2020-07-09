#ifndef __VECTOR_BLOCK__
#define __VECTOR_BLOCK__

#include <basis_block/basis_block.h>
#include <stdlib.h>

struct _vector_block_;
typedef struct _vector_block_ *vector_block_t;

vector_block_t new_vector_block(const char *base_directory,
				const basis_block_t basis_block);

void load_vector_block_elements(vector_block_t vector_block);

void save_vector_block_elements(vector_block_t vector_block);

size_t get_neutron_dimension(const vector_block_t vector_block);

size_t get_proton_dimension(const vector_block_t vector_block);

double *get_vector_block_elements(const vector_block_t vector_block);

void free_vector_block(vector_block_t vector_block);

#endif
