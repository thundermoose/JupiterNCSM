#ifndef __MATRIX_BLOCK__
#define __MATRIX_BLOCK__

#include <basis_block/basis_block.h>
#include <sub_basis_block/sub_basis_block.h>

struct _matrix_block_;
typedef struct _matrix_block_ *matrix_block_t;

matrix_block_t new_matrix_block(size_t block_id,
				const char *base_directory);

double *get_matrix_block_elements(const matrix_block_t matrix_block);

size_t get_neutron_matrix_dimension(const matrix_block_t matrix_block);

size_t get_proton_matrix_dimension(const matrix_block_t matrix_block);

void free_matrix_block(matrix_block_t matrix_block);

#endif
