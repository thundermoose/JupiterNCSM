#include <memory_manager/memory_manager.h>
#include <string_tools/string_tools.h>
#include <radix_sort/radix_sort.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>

typedef enum
{
	UNKNOWN,
	VECTOR_BLOCK,
	MATRIX_BLOCK,
	INDEX_LIST
} array_type_t;

typedef struct
{
	array_type_t type;
	void *primary_array;
	void *secondary_array;
	size_t needed_by_instruction;
	size_t size_array;
	size_t array_id;
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
	execution_order_t execution_order;
	size_t size_current_loaded_memory;
	size_t maximum_loaded_memory;
	size_t *calculation_threads_positions;
	size_t *candidate_arrays_workspace;
	omp_lock_t calculation_threads_positoins_lock;
	omp_lock_t request_lock;
};

static
void initialize_arrays(memory_manager_t manager);

static
void initialize_multi_thread_environment(memory_manager_t manager);

static
size_t get_size_of_unloaded_arrays(memory_manager_t manager,
				   execution_instruction_t instruction);

static
void unload_not_needed_arrays(memory_manager_t manager,
			      size_t min_size_to_unload);

static
void load_needed_arrays(memory_manager_t manager,
			execution_instruction_t instruction);

static
void load_array(memory_manager_t manager,
		size_t array_id);

static
void unload_array(memory_manager_t manager,
		  size_t array_id);

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table,
				    execution_order_t execution_order)
{
	memory_manager_t manager =
	       	(memory_manager_t)malloc(sizeof(struct _memory_manager_));
	manager->num_arrays = get_num_arrays(combination_table);
	manager->all_arrays = (array_t*)calloc(manager->num_arrays,
					       sizeof(array_t));
	manager->candidate_arrays_workspace =
	       	(size_t*)malloc(manager->num_arrays*sizeof(size_t));
	manager->input_vector_base_directory = input_vector_base_directory;
	manager->output_vector_base_directory = output_vector_base_directory;
	manager->index_list_base_directory = index_list_base_directory;
	manager->matrix_base_directory = matrix_base_directory;
	manager->combination_table = combination_table;
	manager->execution_order = execution_order;
	omp_init_lock(&manager->calculation_threads_positoins_lock);
	omp_set_lock(&manager->calculation_threads_positoins_lock);
	initialize_arrays(manager);
	return manager;
}

void launch_memory_manager_thread(memory_manager_t manager)
{
	initialize_multi_thread_environment(manager);
	execution_order_iterator_t instructions =
		get_execution_order_iterator(manager->execution_order);
	while (has_next_instruction(instructions))
	{
		execution_instruction_t instruction =
			next_instruction(instructions);
		size_t size_to_load = 
			get_size_of_unloaded_arrays(manager,instruction);
		if (size_to_load + manager->size_current_loaded_memory >
		    manager->maximum_loaded_memory)
		{
			size_t min_size_to_unload = 
				size_to_load +
			       	manager->size_current_loaded_memory -
				manager->maximum_loaded_memory;
			unload_not_needed_arrays(manager,min_size_to_unload);
		}
		load_needed_arrays(manager,instruction);
	}

}

void begin_instruction(memory_manager_t manager,
		       execution_instruction_t instruction)
{
	omp_set_lock(&manager->calculation_threads_positoins_lock);
	size_t thread_id = omp_get_thread_num();
	manager->calculation_threads_positions[thread_id] =
	       	instruction.instruction_index;
	omp_unset_lock(&manager->calculation_threads_positoins_lock);
}

vector_block_t request_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id)
{
}

vector_block_t request_output_vector_block(memory_manager_t manager,
					size_t vector_block_id)
{
}

index_list_t request_index_list(memory_manager_t manager,
			     const size_t index_list_id)
{
}

matrix_block_t request_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id)
{
}

void release_input_vector(memory_manager_t manager, size_t array_id)
{
}

void release_output_vector(memory_manager_t manager, size_t array_id)
{
}

void release_index_list(memory_manager_t manager, size_t array_id)
{
}

void release_matrix_block(memory_manager_t manager, size_t array_id)
{
}

void free_memory_manager(memory_manager_t manager)
{
}

static
void initialize_arrays(memory_manager_t manager)
{
	size_t *array_sizes = get_array_sizes(manager->combination_table);
	for (size_t i = 0; i<manager->num_arrays; i++)
	{
		array_t *current_array = &manager->all_arrays[i];
		current_array->array_id = i+1;
		current_array->size_array = array_sizes[i];
	}	
	free(array_sizes);
	reset_basis_block_iterator(manager->combination_table);
	while (has_next_basis_block(manager->combination_table))
	{
		basis_block_t basis_block =
			next_basis_block(manager->combination_table);

	}
}

static
void initialize_multi_thread_environment(memory_manager_t manager)
{
	size_t num_threads = omp_get_num_threads();
	manager->calculation_threads_positions =
	       	(size_t*)calloc(num_threads*sizeof(size_t));	
	omp_unset_lock(&manager->calculation_threads_positoins_lock);
}

static
size_t get_size_of_unloaded_arrays(memory_manager_t manager,
				   execution_instruction_t instruction)
{
	size_t size_of_unloaded_arrays = 0;
	if (!is_array_loaded(manager,instruction.vector_block_in))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.vector_block_in);
	if (!is_array_loaded(manager,instruction.vector_block_out))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.vector_block_out);
	if (!is_array_loaded(manager,instruction.matrix_element_file))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.matrix_element_file);
	if (instruction.neutron_index != no_index &&
	    !is_array_loaded(manager,instruction.neutron_index))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.neutron_index);
	if (instruction.proton_index != no_index &&
	    !is_array_loaded(manager,instruction.proton_index))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.proton_index);
	return size_of_unloaded_arrays;
}

static
void unload_not_needed_arrays(memory_manager_t manager,
			      size_t min_size_to_unload)
{
	size_t size_unloaded_memory = 0;
	while (size_unloaded_memory < min_size_to_unload)
	{

		omp_set_lock(&manager->request_lock);
		size_t last_thread = get_last_thread(manager);
		size_t *candidate_arrays = manager->candidate_arrays_workspace;
		size_t num_candidate_arrays = 0;
		for (size_t i = 0; i < manager->num_arrays; i++)
		{
			array_t *candidate = manager->all_arrays+i;
			if ((candidate->primary_array != NULL || 
			     candidate->secondary_array != NULL) &&
			    candidate->needed_by_instruction < last_thread)
				candidate_arrays[num_candidate_arrays++] =
					candidate->array_id;
		}
		omp_unset_lock(&manager->request_lock);
		for (size_t i = 0; i < num_candidate_arrays; i++)
		{
			size_t size_array =
			       	get_size_array(manager, candidate_arrays[i]);
			unload_array(manager,candidate_arrays[i]);
			size_unloaded_memory += size_array;
			if (size_unloaded_memory > min_size_to_unload)
				break;
		}
		usleep(10);
	}
}

static
void load_needed_arrays(memory_manager_t manager,
			execution_instruction_t instruction)
{
	if (!is_array_loaded(manager,instruction.vector_block_in))
		load_array(manager,instruction.vector_block_in);
	if (!is_array_loaded(manager,instruction.vector_block_out))
		load_array(manager,instruction.vector_block_out);
	if (!is_array_loaded(manager,instruction.matrix_element_file))
		load_array(manager,instruction.matrix_element_file);
	if (instruction.neutron_index != no_index && 
	    !is_array_loaded(manager,instruction.neutron_index))
		load_array(manager,instruction.neutron_index);
	if (instruction.proton_index != no_index && 
	    !is_array_loaded(manager,instruction.proton_index))
		load_array(manager,instruction.proton_index);
}

static
void load_array(memory_manager_t manager,
		size_t array_id)
{

}
