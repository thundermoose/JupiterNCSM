#include <matrix_energy_block/matrix_energy_block.h>

struct _matrix_energy_block_
{
	matrix_block_setting_t *matrix_blocks;
	size_t length_block;
	size_t iterator_index;
};

matrix_energy_block_t new_matrix_energy_block(matrix_block_setting_t *head,
					      size_t length_block)
{
	matrix_energy_block_t block =
	       	(matrix_energy_block_t)
		malloc(sizeof(struct _matrix_energy_block_));
	block->matrix_blocks = head;
	block->length_block = length_block;
	block->iterator_index = 0;
	return block;
}

int has_next_energy_matrix_block(matrix_energy_block_t block)
{
	return block->iterator_index < block->length_block;
}

matrix_block_setting_t next_energy_matrix_block(matrix_energy_block_t block)
{
	return block->matrix_blocks[block->iterator_index++];
}

int get_total_isospin(matrix_energy_block_t block)
{
	return count_neutrons(block->matrix_blocks->type) -
		count_protons(block->matrix_blocks->type);
}

int get_proton_energy_bra(matrix_energy_block_t block)
{
	return block->matrix_blocks->depth_protons +
	       	block->matrix_blocks->difference_energy_protons;
}

int get_proton_energy_ket(matrix_energy_block_t block)
{
	return block->matrix_blocks->depth_protons;
}

int get_neutron_energy_bra(matrix_energy_block_t block)
{
	return block->matrix_blocks->depth_neutrons +
	       	block->matrix_blocks->difference_energy_neutrons;
}

int get_neutron_energy_ket(matrix_energy_block_t block)
{
	return block->matrix_blocks->depth_neutrons;
}

int get_bra_energy(matrix_energy_block_t block)
{
	return get_proton_energy_bra(block) + get_neutron_energy_bra(block);
}

int get_ket_energy(matrix_energy_block_t block)
{
	return get_proton_energy_ket(block) + get_neutron_energy_ket(block);
}

void free_matrix_energy_block(matrix_energy_block_t block)
{
	free(block);
}
