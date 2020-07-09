#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

#include <stdlib.h>
#include <vector_block/vector_block.h>
#include <index_list/index_list.h>
#include <matrix_block/matrix_block.h>
#include <combination_table/combination_table.h>

struct _memory_manager_;
typedef struct _memory_manager_ *memory_manager_t;

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table);

vector_block_t load_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id);

vector_block_t load_output_vector_block(memory_manager_t manager,
					size_t vector_block_id);

index_list_t load_index_list(memory_manager_t manager,
			     size_t index_list_id,
			     int sign);

matrix_block_t load_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id);

void unload_input_vector_block(memory_manager_t manager,
			       size_t vector_block_id);

void unload_output_vector_block(memory_manager_t manager,
				size_t vector_block_id);

void unload_index_list(memory_manager_t manager,
		       size_t index_list_id,
		       int sign);

void unload_matrix_block(memory_manager_t manager,
			 size_t matrix_block_id);

void free_memory_manager(memory_manager_t manager);

#endif
