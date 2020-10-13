#ifndef __TRANSFORMED_BLOCK__
#define __TRANSFORMED_BLOCK__

struct _transformed_block_;
typedef struct _transformed_block_ *transformed_block_t;

transformed_block_t new_empty_transformed_block(matrix_energy_block_t block);

void set_transformed_block_ket_basis(matrix_energy_block_t block,
				     M_Scheme_3p_Basis *ket_basis);

void set_transformed_block_bra_basis(matrix_energy_block_t block,
				     M_Scheme_3p_Basis *bra_basis);

void set_transforemd_block_m_block(transformed_block_t block,
				   transformed_3nf_M_block_t M_block,
				   size_t index);

mercury_matrix_block_t get_3nf_mercury_matrix(transformed_block_t block);


void free_transformed_block(transformed_block_t block);

#endif
