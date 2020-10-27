#include <memory_manager/memory_manager.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>

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
	omp_lock_t array_lock;
	size_t num_current_users_primary_array;
	omp_lock_t num_current_users_primary_array_lock;
	size_t num_current_users_secondary_array;
	omp_lock_t num_current_users_secondary_array_lock;
} array_t;

struct _memory_manager_
{	
	array_t *all_arrays;
	size_t *array_sizes;
	size_t num_arrays;
	char *input_vector_base_directory;
	char *output_vector_base_directory;
	char *index_list_base_directory;
	char *matrix_base_directory;
	combination_table_t combination_table;
	execution_order_t execution_order;
};

static
void load_necessary_data(memory_manager_t manager,
			 execution_instruction_t instruction);

static
size_t compute_needed_memory_to_load(memory_manager_t manager,
				     execution_instruction_t instruction);

static
void unload_non_needed_arrays(memory_manager_t manager,
			      size_t size_to_unload,
			      execution_instruction_t instruction);

static
array_t *fetch_array(memory_manager_t manager,size_t array_id);

static
void unload_array(memory_manager_t manager,
		  size_t id);

static
void wait_for_primary_release(array_t *array);

static
void wait_for_secondary_release(array_t *array);

static
void use_primary_array(array_t *array);

static
void use_secondary_array(array_t *array);

static
void unuse_primary_array(array_t *array);

static
void unuse_secondary_array(array_t *array);

static
int is_primary_array_inuse(array_t *array);

static
int is_secondary_array_inuse(array_t *array);

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table,
				    execution_order_t execution_order)
{
	memory_manager_t manager =
		(memory_manager_t)malloc(sizeof(struct _memory_manager_));
	const size_t num_arrays = get_num_arrays(combination_table);
	manager->all_arrays = (array_t*)calloc(num_arrays,
					       sizeof(array_t));
	manager->num_arrays = num_arrays;
	for (size_t i = 0; i < manager->num_arrays; i++)
		omp_init_lock(&manager->all_arrays[i].array_lock);
	manager->input_vector_base_directory =
		copy_string(input_vector_base_directory);
	manager->output_vector_base_directory =
		copy_string(output_vector_base_directory);
	manager->index_list_base_directory =
		copy_string(index_list_base_directory);
	manager->matrix_base_directory =
		copy_string(matrix_base_directory);
	manager->combination_table = combination_table;
	manager->execution_order = execution_order;
	return manager;
}

void launch_memory_manager_thread(memory_manager_t manager)
{
	size_t thread_id = omp_get_thread_num();
	if (thread_id > 0)
	{
		// Give the memory manager time to load
		// data before the data is needed
		usleep(50);
		return;
	}
	execution_order_iterator_t instructions_to_load = 
		get_execution_order_iterator(manager->execution_order);
	while (has_next_instruction(instructions_to_load))
	{
		execution_instruction_t instruction =
			next_instruction(instruction);
		load_necessary_data(manager,instruction);
	}
}

vector_block_t request_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id)
{
	array_t *current_array = fetch_array(manager,vector_block_id);
	omp_set_lock(&current_array->array_lock);
	use_primary_array(current_array);
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
		omp_unset_lock(&current_array->array_lock);
		return (vector_block_t) current_array->primary_array;
	}
	else if (current_array->type != VECTOR_BLOCK)
	{
		error("Array %lu is not a vector block\n",
		      vector_block_id);	      
	}
	else
	{
		omp_unset_lock(&current_array->array_lock);
		return (vector_block_t) current_array->primary_array;	
	}
}

vector_block_t request_output_vector_block(memory_manager_t manager,
					size_t vector_block_id)
{
	array_t *current_array = fetch_array(manager,vector_block_id);
	omp_set_lock(&current_array->array_lock);
	use_secondary_array(current_array);
	if (current_array->secondary_array == NULL)
	{
		basis_block_t basis_block =
			get_basis_block(manager->combination_table,
					vector_block_id);
		current_array->secondary_array = 
			(void*)
			new_output_vector_block(manager->output_vector_base_directory,
					       basis_block);
		current_array->type = VECTOR_BLOCK;
		omp_unset_lock(&current_array->array_lock);
		return (vector_block_t) current_array->secondary_array;
	}
	else if (current_array->type != VECTOR_BLOCK)
	{
		error("Array %lu is not a vector block\n",
		      vector_block_id);	      
	}
	else
	{
		omp_unset_lock(&current_array->array_lock);
		return (vector_block_t) current_array->secondary_array;	
	}
}

index_list_t request_index_list(memory_manager_t manager,
			     const size_t index_list_id)
{
	log_entry("load_index_list(%lu)",
		  index_list_id);
	array_t *current_array = fetch_array(manager,index_list_id);
	omp_set_lock(&current_array->array_lock);
	use_primary_array(current_array);
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
	omp_unset_lock(&current_array->array_lock);
	return (index_list_t)index_list;	
}

matrix_block_t request_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id)
{
	array_t *current_array = fetch_array(manager,matrix_block_id);
	omp_set_lock(&current_array->array_lock);
	use_primary_array(current_array);
	if (current_array->primary_array == NULL)
	{
		matrix_block_t matrix_block =
			new_matrix_block(matrix_block_id,
					 manager->matrix_base_directory);
		current_array->primary_array = (void*)matrix_block;
		current_array->type = MATRIX_BLOCK;
		omp_unset_lock(&current_array->array_lock);
		return matrix_block;
	}
	else if (current_array->type != MATRIX_BLOCK)
	{
		error("Array %lu is not a matrix block\n",
		      matrix_block_id);
	}
	else
	{
		omp_unset_lock(&current_array->array_lock);
		return (matrix_block_t)
			current_array->primary_array;
	}
}

void release_input_vector(memory_manager_t manager, size_t array_id)
{
	unuse_primary_array(fetch_array(manager,array_id));
}

void release_output_vector(memory_manager_t manager, size_t array_id)
{
	unuse_secondary_array(fetch_array(manager,array_id));
}

void release_index_list(memory_manager_t manager, size_t array_id)
{
	unuse_primary_array(fetch_array(manager,array_id));
}

void release_matrix_block(memory_manager_t manager, size_t array_id)
{
	unuse_primary_array(fetch_array(manager,array_id));
}

void unload_input_vector_block(memory_manager_t manager, 
			       size_t vector_block_id)
{
	array_t *current_array = fetch_array(manager,vector_block_id);
	omp_set_lock(&current_array->array_lock);
	wait_for_primary_release(current_array);
	if (current_array->type != VECTOR_BLOCK)
		error("The array %lu is not a vector block\n",
		      vector_block_id);
	if (current_array->primary_array == NULL)
		error("The input vector block %lu is not loaded\n",
		      vector_block_id);	      
	free_vector_block((vector_block_t)current_array->primary_array);
	current_array->primary_array = NULL;
	omp_unset_lock(&current_array->array_lock);
}

void unload_output_vector_block(memory_manager_t manager,
				size_t vector_block_id)
{
	array_t *current_array = fetch_array(manager,vector_block_id);
	omp_set_lock(&current_array->array_lock);
	wait_for_secondary_release(current_array);
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
	omp_unset_lock(&current_array->array_lock);
}

void unload_index_list(memory_manager_t manager, 
		       size_t index_list_id)
{
	log_entry("unload_index_list(%lu)",
		  index_list_id);
	array_t *current_array = fetch_array(manager,index_list_id);
	omp_set_lock(&current_array->array_lock);
	wait_for_primary_release(current_array);
	if (current_array->type != INDEX_LIST)
		error("Array %lu is not an index list\n",
		      index_list_id);
	if (current_array->primary_array == NULL)
		error("Index list %lu is not loaded\n",
		      index_list_id);
	free_index_list((index_list_t)current_array->primary_array);
	current_array->primary_array = NULL;
	omp_unset_lock(&current_array->array_lock);
}

void unload_matrix_block(memory_manager_t manager, 
			 size_t matrix_block_id)
{
	array_t *current_array = fetch_array(manager,matrix_block_id);
	omp_set_lock(&current_array->array_lock);
	wait_for_primary_release(current_array);
	if (current_array->type != MATRIX_BLOCK)
		error("Array %lu is not a matrix block\n",
		      matrix_block_id);
	if (current_array->primary_array == NULL)
		error("Matrix block %lu is not loaded\n",
		      matrix_block_id);
	free_matrix_block((matrix_block_t)current_array->primary_array);
	current_array->primary_array = NULL;
	omp_unset_lock(&current_array->array_lock);
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
void load_necessary_data(memory_manager_t manager,
			 execution_instruction_t instruction)
{
	size_t total_needed_memory_to_load = 
		compute_needed_memory_to_load(manager,instruction);
	if (manager->size_current_loaded_memory + 
	    total_needed_memory_to_load >= manager->maximum_loaded_memory)
	{
		size_t memory_to_unload = 
			manager->size_current_loaded_memory + 
			total_needed_memory_to_load -
			manager->maximum_loaded_memory;
		unload_non_needed_arrays(manager,memory_to_unload,instruction);
	}
	if (!is_array_loaded(manager,instruction.vector_block_in))
		load_vector_array(manager,instruction.vector_block_in);
	if (!is_array_loaded(manager,instruction.vector_block_out))
		load_vector_array(manager,instruction.vector_block_out);
	if (!is_array_loaded(manager,instruction.matrix_element_file))
		load_matrix_array(manager,instruction.matrix_element_file);
	if (!is_array_loaded(manager,instruction.neutron_index))
		load_index_array(manager,instruction.neutron_index);
	if (!is_array_loaded(manager,instruction.proton_index))
		load_index_array(manager,instruction.proton_index);
}

static
size_t compute_needed_memory_to_load(memory_manager_t manager,
				     execution_instruction_t instruction)
{
	size_t needed_memory = 0;
	if (!is_array_loaded(manager,instruction.vector_block_in))
		needed_memory +=
			manager->array_sizes[instruction.vector_block_in-1];	    
	if (instruction.vector_block_out != instruction.vector_block_in &&
	    !is_array_loaded(instruction.vector_block_out))
		needed_memory +=
			manager->array_sizes[instruction.vector_block_out-1];
	if (!is_array_loaded(manager,instruction.matrix_element_file))
		needed_memory +=
			manager->array_sizes[instruction.matrix_element_file-1];
	if (!is_array_loaded(manager,instruction.neutron_index) &&
	    instruction.neutron_index != no_index)
		needed_memory +=
			manager->array_sizes[instruction.neutron_index-1];
	if (!is_array_loaded(manager,instruction.proton_index) &&
	    instruction.proton_index != no_index)
		needed_memory +=
			manager->array_sizes[instruction.proton_index-1];
	return needed_memory;
}

static
void unload_non_needed_arrays(memory_manager_t manager,
			      size_t size_to_unload,
			      execution_instruction_t instruction)
{

}

static
array_t *fetch_array(memory_manager_t manager,size_t array_id)
{
	assert(array_id <= manager->num_arrays);
	array_t *fetched_array = NULL;
#pragma omp critical (fetch_array)
	{
		fetched_array = &manager->all_arrays[array_id - 1];
	}
	return fetched_array;
}

static
void unload_array(memory_manager_t manager,
		  size_t id)
{
	array_t array = *fetch_array(manager,id);
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

static
void wait_for_primary_release(array_t *array)
{
	while (is_primary_array_inuse(array))
		usleep(5);
}

static
void wait_for_secondary_release(array_t *array)
{
	while (is_secondary_array_inuse(array))
		usleep(5);
}

static
void use_primary_array(array_t *array)
{
	omp_set_lock(&array->num_current_users_primary_array_lock);
	array->num_current_users_primary_array++;
	omp_unset_lock(&array->num_current_users_primary_array_lock);
}

static
void unuse_primary_array(array_t *array)
{
	omp_set_lock(&array->num_current_users_primary_array_lock);
	array->num_current_users_primary_array--;
	omp_unset_lock(&array->num_current_users_primary_array_lock);
}

static
void use_secondary_array(array_t *array)
{
	omp_set_lock(&array->num_current_users_secondary_array_lock);
	array->num_current_users_secondary_array++;
	omp_unset_lock(&array->num_current_users_secondary_array_lock);
}

static
void unuse_secondary_array(array_t *array)
{
	omp_set_lock(&array->num_current_users_secondary_array_lock);
	array->num_current_users_secondary_array--;
	omp_unset_lock(&array->num_current_users_secondary_array_lock);
}

static
int is_primary_array_inuse(array_t *array)
{
	int inuse = 0;
	omp_set_lock(&array->num_current_users_primary_array_lock);
	inuse = array->num_current_users_primary_array > 0;
	omp_unset_lock(&array->num_current_users_primary_array_lock);
	return inuse;
}

static
int is_secondary_array_inuse(array_t *array)
{
	int inuse = 0;
	omp_set_lock(&array->num_current_users_secondary_array_lock);
	inuse = array->num_current_users_secondary_array > 0;
	omp_unset_lock(&array->num_current_users_secondary_array_lock);
	return inuse;
}
