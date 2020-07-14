#include <memory_manager/memory_manager.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>

typedef enum
{
	UNKNOWN = 0,
	VECTOR_BLOCK = 1,
	INDEX_LIST = 2,
	MATRIX_BLOCK = 3
} array_type_t;
typedef struct
{
	array_type_t type;
	void *primary_array;
	void *secondary_array;
} array_t;
struct _memory_manager_
{	
	array_t *all_arrays;
	size_t num_arrays;
	char *input_vector_base_directory;
	char *output_vector_base_directory;
	char *index_list_base_directory;
	char *matrix_base_directory;
	combination_table_t combination_table;
};

static
void unload_array(memory_manager_t manager,
		  size_t id);

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table)
{
	memory_manager_t manager =
		(memory_manager_t)malloc(sizeof(struct _memory_manager_));
	const size_t num_arrays = get_num_arrays(combination_table);
	manager->all_arrays = (array_t*)calloc(num_arrays,
					       sizeof(array_t));
	manager->num_arrays = num_arrays;
	manager->input_vector_base_directory =
		copy_string(input_vector_base_directory);
	manager->output_vector_base_directory =
		copy_string(output_vector_base_directory);
	manager->index_list_base_directory =
		copy_string(index_list_base_directory);
	manager->matrix_base_directory =
		copy_string(matrix_base_directory);
	manager->combination_table = combination_table;
	return manager;
}

vector_block_t load_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id)
{
	assert(vector_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[vector_block_id-1];
	if (current_array->primary_array == NULL)
	{
		basis_block_t basis_block =
			get_basis_block(manager->combination_table,
					vector_block_id);
		current_array->primary_array = 
			(void*)
			new_vector_block(manager->input_vector_base_directory,
					 basis_block);
		current_array->type = VECTOR_BLOCK;
		return (vector_block_t) current_array->primary_array;
	}
	else if (current_array->type != VECTOR_BLOCK)
	{
		error("Array %lu is not a vector block\n",
		      vector_block_id);	      
	}
	else
	{
		return (vector_block_t) current_array->primary_array;	
	}
}

vector_block_t load_output_vector_block(memory_manager_t manager,
					size_t vector_block_id)
{
	assert(vector_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[vector_block_id-1];
	if (current_array->secondary_array == NULL)
	{
		basis_block_t basis_block =
			get_basis_block(manager->combination_table,
					vector_block_id);
		current_array->secondary_array = 
			(void*)
			new_vector_block(manager->output_vector_base_directory,
					 basis_block);
		current_array->type = VECTOR_BLOCK;
		return (vector_block_t) current_array->secondary_array;
	}
	else if (current_array->type != VECTOR_BLOCK)
	{
		error("Array %lu is not a vector block\n",
		      vector_block_id);	      
	}
	else
	{
		return (vector_block_t) current_array->secondary_array;	
	}
}

index_list_t load_index_list(memory_manager_t manager,
			     const size_t index_list_id)
{
	log_entry("load_index_list(%lu)",
		  index_list_id);
	assert(index_list_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[index_list_id-1];
	index_list_t index_list = current_array->primary_array; 
	log_entry("current_array->primary_array = %p",
		  current_array->primary_array);
	log_entry("current_array->secondary_array = %p",
		  current_array->secondary_array);
	const char *base_directory = manager->index_list_base_directory;
	if (index_list == NULL)
		index_list = 
			new_index_list_from_id(base_directory,
					       index_list_id);
	void **array = &current_array->primary_array; 
	if (*array == NULL)
	{
		*array = (void*)index_list;
		current_array->type = INDEX_LIST;
	}
	else if (current_array->type != INDEX_LIST)
	{
		error("Array %lu is not an index list\n",
		      index_list_id);	      
	}
	log_entry("current_array->primary_array = %p",
		  current_array->primary_array);
	return (index_list_t)index_list;	
}

matrix_block_t load_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id)
{
	assert(matrix_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[matrix_block_id-1];	
	if (current_array->primary_array == NULL)
	{
		matrix_block_t matrix_block =
			new_matrix_block(matrix_block_id,
					 manager->matrix_base_directory);
		current_array->primary_array = (void*)matrix_block;
		current_array->type = MATRIX_BLOCK;
		return matrix_block;
	}
	else if (current_array->type != MATRIX_BLOCK)
	{
		error("Array %lu is not a matrix block\n",
		      matrix_block_id);
	}
	else
	{
		return (matrix_block_t)
			current_array->primary_array;
	}
}

void unload_input_vector_block(memory_manager_t manager,
			       size_t vector_block_id)
{
	assert(vector_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[vector_block_id-1];
	if (current_array->type != VECTOR_BLOCK)
		error("The array %lu is not a vector block\n",
		      vector_block_id);
	if (current_array->primary_array == NULL)
		error("The input vector block %lu is not loaded\n",
		      vector_block_id);	      
	free_vector_block((vector_block_t)current_array->primary_array);
	current_array->primary_array = NULL;
}

void unload_output_vector_block(memory_manager_t manager,
				size_t vector_block_id)
{
	assert(vector_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[vector_block_id-1];
	if (current_array->type != VECTOR_BLOCK)
		error("The array %lu is not a vector block\n",
		      vector_block_id);
	if (current_array->secondary_array == NULL)
		error("The input vector block %lu is not loaded\n",
		      vector_block_id);	      
	vector_block_t vector_block =
	       	(vector_block_t)current_array->secondary_array;
	save_vector_block_elements(vector_block);
	free_vector_block(vector_block);
	current_array->secondary_array = NULL;
}

void unload_index_list(memory_manager_t manager,
		       size_t index_list_id)
{
	log_entry("unload_index_list(%lu)",
		  index_list_id);
	assert(index_list_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[index_list_id-1];
	if (current_array->type != INDEX_LIST)
		error("Array %lu is not an index list\n",
		      index_list_id);
	if (current_array->primary_array == NULL)
		error("Index list %lu is not loaded\n",
		      index_list_id);
	free_index_list((index_list_t)current_array->primary_array);
	current_array->primary_array = NULL;
}

void unload_matrix_block(memory_manager_t manager,
			 size_t matrix_block_id)
{
	assert(matrix_block_id <= manager->num_arrays);
	array_t *current_array = &manager->all_arrays[matrix_block_id-1];
	if (current_array->type != MATRIX_BLOCK)
		error("Array %lu is not a matrix block\n",
		      matrix_block_id);
	if (current_array->primary_array == NULL)
		error("Matrix block %lu is not loaded\n",
		      matrix_block_id);
	free_matrix_block((matrix_block_t)current_array->primary_array);
	current_array->primary_array = NULL;
}

void free_memory_manager(memory_manager_t manager)
{
	log_entry("Freeing the memory manager");
	for (size_t i = 0; i<manager->num_arrays; i++)
	{
		unload_array(manager,i+1);
	}
	free(manager->all_arrays);
	free(manager->input_vector_base_directory);
	free(manager->output_vector_base_directory);
	free(manager->index_list_base_directory);
	free(manager->matrix_base_directory);
	free(manager);
}

static
void unload_array(memory_manager_t manager,
		  size_t id)
{
	array_t array = manager->all_arrays[id-1];
	switch (array.type)
	{
		case VECTOR_BLOCK:
			if (array.primary_array != NULL)
				unload_input_vector_block(manager,id);
			if (array.secondary_array != NULL)
				unload_output_vector_block(manager,id);
			break;
		case INDEX_LIST:
			if (array.primary_array != NULL)
				unload_index_list(manager,id);
			break;
		case MATRIX_BLOCK:
			if (array.primary_array != NULL)
				unload_matrix_block(manager,id);
			break;
		default:
			break;	
	}
	array.primary_array = NULL;
	array.secondary_array = NULL;
}
