#ifndef __TRANSFORMED_BLOCK__
#define __TRANSFORMED_BLOCK__

#include <single_particle_basis/single_particle_basis.h>
#include <bases/m_scheme_3p_basis.h>
#include <transformed_3nf_M_block/transformed_3nf_M_block.h>
#include <matrix_energy_block/matrix_energy_block.h>
#include <mercury_matrix_block/mercury_matrix_block.h>

struct _transformed_block_;
typedef struct _transformed_block_ *transformed_block_t;

transformed_block_t 
new_empty_transformed_block(matrix_energy_block_t energy_block,
			    char *index_list_path,
			    single_particle_basis_t single_particle_basis);

void set_transformed_block_ket_basis(transformed_block_t block,
				     M_Scheme_3p_Basis *ket_basis);

void set_transformed_block_bra_basis(transformed_block_t block,
				     M_Scheme_3p_Basis *bra_basis);

void set_transformed_block_m_block(transformed_block_t block,
				   transformed_3nf_M_block_t M_block,
				   size_t index);

mercury_matrix_block_t get_3nf_mercury_matrix(transformed_block_t block,
					      matrix_block_setting_t settings);


void free_transformed_block(transformed_block_t block);

#endif
