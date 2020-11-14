#ifndef __MEMORY_MANAGER__
#define __MEMORY_MANAGER__

#include <stdlib.h>
#include <vector_block/vector_block.h>
#include <index_list/index_list.h>
#include <matrix_block/matrix_block.h>
#include <combination_table/combination_table.h>
#include <execution_order/execution_order.h>

struct _memory_manager_;
typedef struct _memory_manager_ *memory_manager_t;

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table,
				    execution_order_t execution_order);

void initialize_multi_thread_environment(memory_manager_t manager);

void begin_instruction(memory_manager_t manager,
		       execution_instruction_t instruction);


vector_block_t request_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id);

vector_block_t request_output_vector_block(memory_manager_t manager,
					size_t vector_block_id);

index_list_t request_index_list(memory_manager_t manager,
			     size_t index_list_id);

matrix_block_t request_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id);

void release_input_vector(memory_manager_t manager, size_t array_id);

void release_output_vector(memory_manager_t manager, size_t array_id);

void release_index_list(memory_manager_t manager, size_t array_id);

void release_matrix_block(memory_manager_t manager, size_t array_id);

void free_memory_manager(memory_manager_t manager);

#endif
