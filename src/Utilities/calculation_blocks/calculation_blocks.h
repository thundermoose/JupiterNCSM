#ifndef __CALCULATION_BLOCKS__
#define __CALCULATION_BLOCKS__

#include <stdlib.h>
#include <stdio.h>

struct _calculation_blocks_;
typedef struct _calculation_blocks_ *calculation_blocks_t;

typedef struct _calculation_block_
{
	size_t num_particles;
	size_t vector_block_in;
	size_t vector_block_out;
	size_t primary_index_list;
	size_t secondary_index_list;
	size_t matrix_element_block;
	size_t num_multiplications;
} calculation_block_t;

calculation_blocks_t new_calculation_blocks(calculation_block_t *blocks,
					    size_t num_blocks);

calculation_blocks_t parse_calculation_blocks(FILE* combination_file);

size_t get_num_calculation_blocks(calculation_blocks_t blocks);

calculation_block_t *get_calculation_blocks_array(calculation_blocks_t block);

void reset_calculation_block_iterator(calculation_blocks_t blocks);

int has_next_calculation_block(calculation_blocks_t blocks);

calculation_block_t next_calculation_block(calculation_blocks_t blocks);

void free_calculation_blocks(calculation_blocks_t calculation_blocks);

size_t get_calculation_block_num_particles(calculation_block_t *block);

#endif
