#include <transformed_block/transformed_block.h>
#include <connection_list/connection_list.h>
#include <assert.h>

struct _transformed_block_
{
	matrix_energy_block_t energy_block;
	char *index_list_path;
	single_particle_basis_t single_particle_basis;
	M_Scheme_3p_Basis *ket_basis;
	M_Scheme_3p_Basis *bra_basis;
	int M_min;
	transformed_3nf_M_block_t *M_blocks;
	size_t num_M_blocks;
};

transformed_block_t 
new_empty_transformed_block(matrix_energy_block_t energy_block,
			    char *index_list_path,
			    single_particle_basis_t single_particle_basis)
{
	transformed_block_t block =
	       	(transformed_block_t)
		calloc(1,sizeof(struct _transformed_block_));
	
	block->energy_block = energy_block;
	block->index_list_path = index_list_path;
	block->single_particle_basis = single_particle_basis;
	return block;
}

void set_transformed_block_ket_basis(transformed_block_t block,
				     M_Scheme_3p_Basis *ket_basis)
{
	block->ket_basis = ket_basis;
}

void set_transformed_block_bra_basis(transformed_block_t block,
				     M_Scheme_3p_Basis *bra_basis)
{
	block->bra_basis = bra_basis;
}

void set_transformed_block_m_block(transformed_block_t block,
				   transformed_3nf_M_block_t M_block,
				   size_t index)
{
	assert(block->ket_basis);
	assert(block->bra_basis);
	if (block->M_blocks == NULL)
	{
		block->M_min = max(get_min_M(block->ket_basis),
				   get_min_M(block->bra_basis));
		int M_max = min(get_max_M(block->ket_basis),
				get_max_M(block->bra_basis));
		block->num_M_blocks = (M_max - block->M_min)/2+1;
		block->M_blocks = 
			(transformed_3nf_M_block_t*)
			calloc(block->num_M_blocks,
			       sizeof(transformed_3nf_M_block_t));

	}
	block->M_blocks[index] = M_block;
}

int get_transformed_block_min_M(transformed_block_t block)
{
	return block->min_M;
}

size_t get_num_M_blocks(transformed_block_t block)
{
	return block->num_M_blocks;
}

mercury_matrix_block_t get_3nf_mercury_matrix(transformed_block_t block,
					      matrix_block_setting_t settings)
{
	connection_list_t connection_list = 
		read_connection_files(block->index_list_path,
				      settings);
	size_t num_elements = num_connections(connection_list);
	double *elements = (double*)calloc(num_elements,sizeof(double));
	size_t element_index = 0;
	while (has_next_connection(connection_list))
	{
		connection_t current_connection =
			next_connection(connection_list);
		int M = compute_M(current_connection,
				  block->single_particle_basis);
		size_t i = (M-block->min_M)/2;
		M_Scheme_3p_Basis *ket_m_basis = block->blocks[i].ket_basis;
		M_Scheme_3p_Basis *bra_m_basis = block->blocks[i].bra_basis;
		Dens_Matrix *current_matrix = block->blocks[i].matrix;
		int phase = 1;
		size_t ket_index = get_ket_index(current_connection,
						 ket_m_basis,
						 &phase);
		size_t bra_index = get_bra_index(current_connection,
						 bra_m_basis,
						 &phase);
		elements[element_index++] =
			phase*get_dens_matrix_element(current_matrix,
						      bra_index,
						      ket_index);
	}
	free_connection_list(connection_list);
	return new_mercury_matrix_block_from_data(elements,
						  num_elements,
						  settings);
}


void free_transformed_block(transformed_block_t block)
{
	for (size_t i = 0; i < block->num_M_blocks; i++)
	{
		free_m_scheme_3p_basis(block->blocks[i].ket_basis);
		free_m_scheme_3p_basis(block->blocks[i].bra_basis);
		free_dens_matrix(block->blocks[i].matrix);
	}
	if (block->ket_basis)
		free_m_scheme_3p_basis(block->ket_basis);
	if (block->bra_basis)
		free_m_scheme_3p_basis(block->bra_basis);
	free(block);
}
