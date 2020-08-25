#include <transformed_block_manager/transformed_block_manager.h>
#include <bases/m_scheme_2p_basis.h>

struct _transformed_block_manager_
{
	m_scheme_2p_basis_t ket_basis;
	m_scheme_2p_basis_t bra_basis;	

};

transformed_block_manager_t 
new_transformed_block_manager(antoine_2nf_file_t coupled_2nf_data)
{
	return NULL;
}

void decouple_transform_block(transformed_block_manager_t manager,
			      transform_block_settings_t block_settings)
{
}

mercury_matrix_block_t 
get_transformed_matrix_block(transformed_block_manager_t manager,
			     matrix_block_setting_t settings)
{
	return NULL;
}

void free_transformed_block_manager(transformed_block_manager_t manager)
{
}
