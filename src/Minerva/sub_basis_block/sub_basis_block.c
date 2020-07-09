#include <sub_basis_block/sub_basis_block.h>

size_t get_sub_basis_block_dimension(sub_basis_block_t sub_basis_block)
{
	return sub_basis_block.dimension;
}

size_t get_num_particles(sub_basis_block_t sub_basis_block)
{
	return sub_basis_block.num_particles;
}
