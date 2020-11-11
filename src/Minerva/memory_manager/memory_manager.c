#include <memory_manager/memory_manager.h>
#include <string_tools/string_tools.h>
#include <radix_sort/radix_sort.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>

#define min(a,b) ((a) < (b) ? (a) : (b)) 

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
	omp_lock_t array_lock;
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
	omp_lock_t calculation_threads_positions_lock;
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
int is_array_loaded(memory_manager_t manager,
		    size_t array_id);

static
void load_array(memory_manager_t manager,
		size_t array_id);

static
void unload_array(memory_manager_t manager,
		  size_t array_id);

static
void set_needed_by(memory_manager_t manager,
		   size_t array_id,
		   size_t needed_by);

static
size_t get_last_thread(memory_manager_t manager);

static
size_t get_array_size(memory_manager_t manager,
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
	manager->input_vector_base_directory = copy_string(input_vector_base_directory);
	manager->output_vector_base_directory = copy_string(output_vector_base_directory);
	manager->index_list_base_directory = copy_string(index_list_base_directory);
	manager->matrix_base_directory = copy_string(matrix_base_directory);
	manager->combination_table = combination_table;
	manager->execution_order = execution_order;
	omp_init_lock(&manager->calculation_threads_positions_lock);
	omp_set_lock(&manager->calculation_threads_positions_lock);
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
	omp_set_lock(&manager->calculation_threads_positions_lock);
	size_t thread_id = omp_get_thread_num();
	manager->calculation_threads_positions[thread_id] =
	       	instruction.instruction_index;
	omp_unset_lock(&manager->calculation_threads_positions_lock);
}

vector_block_t request_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id)
{
	// No new requests are allowed while unloading
	omp_set_lock(&manager->request_lock);
	omp_unset_lock(&manager->request_lock);
	array_t *array = &manager->all_arrays[vector_block_id-1];
	omp_set_lock(&array->array_lock);
	if (array->type != VECTOR_BLOCK)
		error("%lu is not a vector block\n", vector_block_id);
	vector_block_t output_vector_block = 
		(vector_block_t)array->primary_array;
	omp_unset_lock(&array->array_lock);
	return output_vector_block;
}

vector_block_t request_output_vector_block(memory_manager_t manager,
					size_t vector_block_id)
{
	// No new requests are allowed while unloading
	omp_set_lock(&manager->request_lock);
	omp_unset_lock(&manager->request_lock);
	array_t *array = &manager->all_arrays[vector_block_id-1];
	omp_set_lock(&array->array_lock);
	if (array->type != VECTOR_BLOCK)
		error("%lu is not a vector block\n", vector_block_id);
	vector_block_t output_vector_block = 
		(vector_block_t)array->secondary_array;
	omp_unset_lock(&array->array_lock);
	return output_vector_block;
}

index_list_t request_index_list(memory_manager_t manager,
			     const size_t index_list_id)
{
	// No new requests are allowed while unloading
	omp_set_lock(&manager->request_lock);
	omp_unset_lock(&manager->request_lock);
	array_t *array = &manager->all_arrays[index_list_id-1];
	omp_set_lock(&array->array_lock);
	if (array->type != INDEX_LIST)
		error("%lu is not an index list\n", index_list_id);
	index_list_t output_index_list = 
		(index_list_t)array->primary_array;
	omp_unset_lock(&array->array_lock);
	return output_index_list;
}

matrix_block_t request_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id)
{
	// No new requests are allowed while unloading
	omp_set_lock(&manager->request_lock);
	omp_unset_lock(&manager->request_lock);
	array_t *array = &manager->all_arrays[matrix_block_id-1];
	omp_set_lock(&array->array_lock);
	if (array->type != MATRIX_BLOCK)
		error("%lu is not a matrix block\n", matrix_block_id);
	matrix_block_t output_index_list = 
		(matrix_block_t)array->primary_array;
	omp_unset_lock(&array->array_lock);
	return output_index_list;
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
	for (size_t i = 0; i < manager->num_arrays; i++)
		omp_destroy_lock(&manager->all_arrays[i].array_lock);
	free(manager->all_arrays);
	free(manager->input_vector_base_directory);
	free(manager->output_vector_base_directory);
	free(manager->index_list_base_directory);
	free(manager->matrix_base_directory);
	omp_destroy_lock(&manager->calculation_threads_positions_lock);
	omp_destroy_lock(&manager->request_lock);
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
		omp_init_lock(&current_array->array_lock);
	}	
	free(array_sizes);
	iterator_t basis_blocks =
	       	new_basis_block_iterator(manager->combination_table);
	basis_block_t current_basis_block;
	for (initialize(basis_blocks,&current_basis_block);
	     has_next_element(basis_blocks);
	     next_element(basis_blocks,&current_basis_block))
		manager->all_arrays[current_basis_block.block_id].type = 
			VECTOR_BLOCK;
	free_iterator(basis_blocks);
	iterator_t index_lists = 
		new_index_list_setting_iterator(manager->combination_table);
	index_list_setting_t index_list;
	for (initialize(index_lists,&index_list);
	     has_next_element(index_lists);
	     next_element(index_lists,&index_list))
		manager->all_arrays[index_list.index_list_id].type = 
			INDEX_LIST;
	free_iterator(index_lists);
	iterator_t matrix_blocks =
		new_matrix_block_settings_iterator(manager->combination_table);
	matrix_block_setting_t matrix_block;
	for (initialize(matrix_blocks,&matrix_block);
	     has_next_element(matrix_blocks);
	     next_element(matrix_blocks,&matrix_block))
		manager->all_arrays[matrix_block.matrix_block_id].type =
			MATRIX_BLOCK;
	free_iterator(matrix_blocks);
}

static
void initialize_multi_thread_environment(memory_manager_t manager)
{
	size_t num_threads = omp_get_num_threads();
	manager->calculation_threads_positions =
	       	(size_t*)calloc(num_threads,sizeof(size_t));	
	omp_unset_lock(&manager->calculation_threads_positions_lock);
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
			       	get_array_size(manager, candidate_arrays[i]);
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
	set_needed_by(manager,
		      instruction.vector_block_in,
		      instruction.instruction_index);
	if (!is_array_loaded(manager,instruction.vector_block_out))
		load_array(manager,instruction.vector_block_out);
	set_needed_by(manager,
		      instruction.vector_block_out,
		      instruction.instruction_index);
	if (!is_array_loaded(manager,instruction.matrix_element_file))
		load_array(manager,instruction.matrix_element_file);
	set_needed_by(manager,
		      instruction.matrix_element_file,
		      instruction.instruction_index);
	if (instruction.neutron_index != no_index && 
	    !is_array_loaded(manager,instruction.neutron_index))
		load_array(manager,instruction.neutron_index);
	set_needed_by(manager,
		      instruction.neutron_index,
		      instruction.instruction_index);
	if (instruction.proton_index != no_index && 
	    !is_array_loaded(manager,instruction.proton_index))
		load_array(manager,instruction.proton_index);
	set_needed_by(manager,
		      instruction.proton_index,
		      instruction.instruction_index);
}

static
int is_array_loaded(memory_manager_t manager,
		    size_t array_id)
{
	return manager->all_arrays[array_id-1].primary_array != NULL ||
		manager->all_arrays[array_id-1].secondary_array != NULL;
}

static
void load_array(memory_manager_t manager,
		size_t array_id)
{
	array_t *array = &manager->all_arrays[array_id-1];
	omp_set_lock(&array->array_lock);
	switch(array->type)
	{
	case VECTOR_BLOCK:
		;
		basis_block_t basis_block =
		       	get_basis_block(manager->combination_table,
					array_id);
		array->primary_array =
			(void*)
		       	new_vector_block(manager->input_vector_base_directory,
					 basis_block);
		array->secondary_array = 
			(void*)
			new_output_vector_block
			(manager->output_vector_base_directory,
			 basis_block);
		break;
	case INDEX_LIST:
		array->primary_array =
			(void*)
			new_index_list_from_id
			(manager->index_list_base_directory,
			 array_id);
		break;
	case MATRIX_BLOCK:
		array->primary_array =
			(void*)
			new_matrix_block(array_id,
					 manager->matrix_base_directory);		
		break;
	default:
		error("Can't load unkown array\n");
	}	
	manager->size_current_loaded_memory+=array->size_array;
	omp_unset_lock(&array->array_lock);
}

static
void unload_array(memory_manager_t manager,
		  size_t array_id)
{
	array_t *array = &manager->all_arrays[array_id-1];
	omp_set_lock(&array->array_lock);
	switch(array->type)
	{
	case VECTOR_BLOCK:
		free_vector_block((vector_block_t)array->primary_array);
		save_vector_block_elements((vector_block_t)
					   array->secondary_array);
		free_vector_block((vector_block_t)array->secondary_array);
		break;
	case INDEX_LIST:
		free_index_list((index_list_t)array->primary_array);
		break;
	case MATRIX_BLOCK:
		free_matrix_block((matrix_block_t)array->primary_array);
		break;
	default:
		error("Can't load unkown array\n");
	}	
	manager->size_current_loaded_memory-=array->size_array;
	omp_unset_lock(&array->array_lock);
}

static
void set_needed_by(memory_manager_t manager,
		   size_t array_id,
		   size_t needed_by)
{
	manager->all_arrays[array_id-1].needed_by_instruction = needed_by;
}

static
size_t get_last_thread(memory_manager_t manager)
{
	omp_set_lock(&manager->calculation_threads_positions_lock);
	size_t last_thread = get_num_instructions(manager->execution_order);
	size_t num_threads = omp_get_num_threads();	
	for (size_t i = 1; i < num_threads; i++)
	{
		last_thread = min(last_thread,
				  manager->calculation_threads_positions[i]);
	}
	omp_unset_lock(&manager->calculation_threads_positions_lock);
	return last_thread;
}

static
size_t get_array_size(memory_manager_t manager,
		      size_t array_id)
{
	return manager->all_arrays[array_id-1].size_array;	
}
