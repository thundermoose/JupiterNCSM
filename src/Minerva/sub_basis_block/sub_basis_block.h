#ifndef __SUB_BASIS_BLOCK__
#define __SUB_BASIS_BLOCK__

#include <stdlib.h>

typedef struct
{
	int E;
	int M;
	int depth;
	size_t dimension;
	size_t num_particles;
} sub_basis_block_t;

size_t get_sub_basis_block_dimension(sub_basis_block_t sub_basis_block);

size_t get_num_particles(sub_basis_block_t sub_basis_block);

#endif
