#ifndef __MATRIX_VECTOR_MULTIPLICATION__
#define __MATRIX_VECTOR_MULTIPLICATION__

#include <vector_block/vector_block.h>
#include <matrix_block/matrix_block.h>
#include <index_list/index_list.h>

void multiplication_neutrons(vector_block_t out_block,
				 const vector_block_t in_block,
				 const matrix_block_t block,
				 const index_list_t neutron_list,
				 const int sign);

void multiplication_protons(vector_block_t out_block,
				const vector_block_t in_block,
				const matrix_block_t block,
				const index_list_t proton_list,
				const int sign);

void multiplication_neutrons_protons(vector_block_t out_block,
					 const vector_block_t in_block,
					 const matrix_block_t block,
					 const index_list_t neutron_list,
					 const index_list_t proton_list,
					 const int sign);
#endif
