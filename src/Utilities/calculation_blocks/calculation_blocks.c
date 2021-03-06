#include <calculation_blocks/calculation_blocks.h>
#include <array_builder/array_builder.h>
#include <global_constants/global_constants.h>
#include <string_tools/string_tools.h>
#include <string.h>
#include <error/error.h>

struct _calculation_blocks_
{
	size_t num_blocks;
	calculation_block_t *blocks;
	size_t next_block_index;
};

calculation_blocks_t new_calculation_blocks(calculation_block_t *blocks,
					    size_t num_blocks)
{
	calculation_blocks_t calculation_blocks = 
		(calculation_blocks_t)
		calloc(1,sizeof(struct _calculation_blocks_));
	calculation_blocks->blocks =
	       	(calculation_block_t*)
		malloc(num_blocks*sizeof(calculation_block_t));
	calculation_blocks->num_blocks = num_blocks;
	memcpy(calculation_blocks->blocks,
	       blocks,
	       calculation_blocks->num_blocks*sizeof(calculation_block_t));
	return calculation_blocks;
}

calculation_blocks_t parse_calculation_blocks(FILE* combination_file)
{
	calculation_blocks_t calculation_blocks = 
		(calculation_blocks_t)
		calloc(1,sizeof(struct _calculation_blocks_));
	array_builder_t blocks_builder = 
		new_array_builder((void**)&calculation_blocks->blocks,
				  &calculation_blocks->num_blocks,
				  sizeof(calculation_block_t));
	rewind(combination_file);
	size_t row_buffer_length = 0;
	char *row_buffer = NULL;
	while (!feof(combination_file))
	{
		if (getline(&row_buffer,&row_buffer_length,combination_file)<0)
			break;
		char *calc_block = strstr(row_buffer,"CALCBLOCK:");
		if (calc_block == NULL)
			continue;
		char *size_string = jump_back_word(row_buffer,
						   jump_back_word(row_buffer,
								  calc_block));
		calculation_block_t current_block;
		if (sscanf(calc_block,
			   "CALCBLOCK:%lu: %lu, %lu, %lu, %lu, %lu",
			   &current_block.num_particles,
			   &current_block.vector_block_in,
			   &current_block.vector_block_out,
			   &current_block.primary_index_list,
			   &current_block.secondary_index_list,
			   &current_block.matrix_element_block) != 6)
		{
			current_block.secondary_index_list = no_index;
			if (sscanf(calc_block,
				   "CALCBLOCK:%lu: %lu, %lu, %lu, %lu",
				   &current_block.num_particles,
				   &current_block.vector_block_in,
				   &current_block.vector_block_out,
				   &current_block.primary_index_list,
				   &current_block.matrix_element_block) != 5)
				error("CALCBLOCK \"%s\" is not formated correctly\n",
				      calc_block);
		}
		if (sscanf(size_string,"%lu",&current_block.num_multiplications) != 1)
				error("CALCBLOCK \"%s\" is not formated correctly\n",
				      size_string);
			
		
		append_array_element(blocks_builder,&current_block);
	}
	free(row_buffer);
	free_array_builder(blocks_builder);
	return calculation_blocks;
}

size_t get_num_calculation_blocks(calculation_blocks_t blocks)
{
	return blocks->num_blocks;
}

calculation_block_t *get_calculation_blocks_array(calculation_blocks_t blocks)
{
	calculation_block_t *calculation_blocks =
	       	(calculation_block_t*)
		malloc(blocks->num_blocks*sizeof(calculation_block_t));
	memcpy(calculation_blocks,
	       blocks->blocks,
	       blocks->num_blocks*sizeof(calculation_block_t));
	return calculation_blocks;
}

void reset_calculation_block_iterator(calculation_blocks_t blocks)
{
	blocks->next_block_index = 0;
}

int has_next_calculation_block(calculation_blocks_t blocks)
{
	return blocks->next_block_index < blocks->num_blocks;
}

calculation_block_t next_calculation_block(calculation_blocks_t blocks)
{
	return blocks->blocks[blocks->next_block_index++];
}

void free_calculation_blocks(calculation_blocks_t calculation_blocks)
{
	free(calculation_blocks->blocks);
	free(calculation_blocks);
}

size_t get_calculation_block_num_particles(calculation_block_t *block)
{
	return block->num_particles;
}
