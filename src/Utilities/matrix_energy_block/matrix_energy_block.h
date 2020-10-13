#ifndef __MATRIX_ENERGY_BLOCK__
#define __MATRIX_ENERGY_BLOCK__

struct _matrix_energy_block_;
typedef struct _matrix_energy_block_ *matrix_energy_block_t;

matrix_energy_block_t new_matrix_energy_block(matrix_block_setting_t *head,
					      size_t length_block);

int has_next_energy_matrix_block(matrix_energy_block_t block);

matrix_block_setting_t next_energy_matrix_block(matrix_energy_block_t block);

void free_matrix_energy_block(matrix_energy_block_t block);
#endif
