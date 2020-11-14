#include <memory_manager/memory_manager.h>
#include <string_tools/string_tools.h>
#include <radix_sort/radix_sort.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <omp.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <debug_mode/debug_mode.h>

#define min(a,b) ((a) < (b) ? (a) : (b)) 
#define max(a,b) ((a) > (b) ? (a) : (b)) 

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
	size_t in_use;
	size_t needed_by_instruction;
	size_t size_array;
	size_t array_id;
	// To make sure that no two threads loads the same array
	omp_nest_lock_t loading_array_lock;
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
	double total_wating_time;
	size_t num_waits;
	double min_waiting_time;
	double max_waiting_time;
	omp_lock_t calculation_threads_positions_lock;
};

static
void initialize_arrays(memory_manager_t manager);


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
void wait_til_array_is_loaded(memory_manager_t manager, size_t array_id);

static
int is_array_loaded(memory_manager_t manager, size_t array_id);

static
void load_array(memory_manager_t manager, size_t array_id);

static
void unload_array(memory_manager_t manager, size_t array_id);

static
void set_all_needed_by(memory_manager_t manager,
		       execution_instruction_t instruction);
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
	manager->maximum_loaded_memory = (size_t)(16)<<30;
	manager->size_current_loaded_memory = 0;
	manager->num_waits = 0;
	manager->total_wating_time = 0.0;
	manager->min_waiting_time = INFINITY;
	manager->max_waiting_time = 0.0;
	omp_init_lock(&manager->calculation_threads_positions_lock);
	initialize_arrays(manager);
	return manager;
}

void begin_instruction(memory_manager_t manager,
		       execution_instruction_t instruction)
{
	size_t thread_id = omp_get_thread_num();
	omp_set_lock(&manager->calculation_threads_positions_lock);
	manager->calculation_threads_positions[thread_id] =
		instruction.instruction_index;
	omp_unset_lock(&manager->calculation_threads_positions_lock);
	set_all_needed_by(manager,instruction);
	size_t size_to_load = 
		get_size_of_unloaded_arrays(manager,instruction);
	if (size_to_load == 0)
		return;
	if (size_to_load + manager->size_current_loaded_memory >
	    manager->maximum_loaded_memory)
	{
		size_t min_size_to_unload = 
			size_to_load +
			manager->size_current_loaded_memory -
			manager->maximum_loaded_memory;
		log_entry("Min size to unload: %lu B\n",
			  min_size_to_unload);
		unload_not_needed_arrays(manager,min_size_to_unload);
	}
	load_needed_arrays(manager,instruction);
}

vector_block_t request_input_vector_block(memory_manager_t manager,
				       size_t vector_block_id)
{
	// No new requests are allowed while unloading
	wait_til_array_is_loaded(manager,vector_block_id);
	array_t *array = &manager->all_arrays[vector_block_id-1];
	if (array->type != VECTOR_BLOCK)
		error("%lu is not a vector block\n", vector_block_id);
	vector_block_t output_vector_block = 
		(vector_block_t)array->primary_array;
	array->in_use++;
	return output_vector_block;
}

vector_block_t request_output_vector_block(memory_manager_t manager,
					size_t vector_block_id)
{
	// No new requests are allowed while unloading
	wait_til_array_is_loaded(manager,vector_block_id);
	array_t *array = &manager->all_arrays[vector_block_id-1];
	if (array->type != VECTOR_BLOCK)
		error("%lu is not a vector block\n", vector_block_id);
	vector_block_t output_vector_block = 
		(vector_block_t)array->secondary_array;
	array->in_use++;
	return output_vector_block;
}

index_list_t request_index_list(memory_manager_t manager,
			     const size_t index_list_id)
{
	// No new requests are allowed while unloading
	wait_til_array_is_loaded(manager,index_list_id);
	array_t *array = &manager->all_arrays[index_list_id-1];
	if (array->type != INDEX_LIST)
		error("%lu is not an index list\n", index_list_id);
	index_list_t output_index_list = 
		(index_list_t)array->primary_array;
	array->in_use++;
	return output_index_list;
}

matrix_block_t request_matrix_block(memory_manager_t manager,
				 size_t matrix_block_id)
{
	// No new requests are allowed while unloading
	wait_til_array_is_loaded(manager,matrix_block_id);
	array_t *array = &manager->all_arrays[matrix_block_id-1];
	if (array->type != MATRIX_BLOCK)
		error("%lu is not a matrix block\n", matrix_block_id);
	matrix_block_t output_index_list = 
		(matrix_block_t)array->primary_array;
	array->in_use++;
	return output_index_list;
}

void release_input_vector(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	current_array->in_use--;
}

void release_output_vector(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	current_array->in_use--;
}

void release_index_list(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	current_array->in_use--;
}

void release_matrix_block(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	current_array->in_use--;
}

void free_memory_manager(memory_manager_t manager)
{
	printf("Average wait time: %lg µs\n",
	       manager->total_wating_time/manager->num_waits);
	printf("Min wait time: %lg µs\n",
	       manager->min_waiting_time);
	printf("Max wait time: %lg µs\n",
	       manager->max_waiting_time);
	for (size_t i = 0; i < manager->num_arrays; i++)
	{
		if (is_array_loaded(manager,i+1))
			unload_array(manager,i+1);
		omp_destroy_nest_lock(&manager->all_arrays[i].loading_array_lock);
	}
	free(manager->all_arrays);
	free(manager->input_vector_base_directory);
	free(manager->output_vector_base_directory);
	free(manager->index_list_base_directory);
	free(manager->matrix_base_directory);
	omp_destroy_lock(&manager->calculation_threads_positions_lock);
	free(manager);
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
		current_array->in_use = 0;
		omp_init_nest_lock(&current_array->loading_array_lock);
	}	
	free(array_sizes);
	iterator_t basis_blocks =
	       	new_basis_block_iterator(manager->combination_table);
	basis_block_t current_basis_block;
	for (initialize(basis_blocks,&current_basis_block);
	     has_next_element(basis_blocks);
	     next_element(basis_blocks,&current_basis_block))
	{
		log_entry("Array %lu is a vector block\n",
		       current_basis_block.block_id);
		manager->all_arrays[current_basis_block.block_id-1].type = 
			VECTOR_BLOCK;
	}
	free_iterator(basis_blocks);
	iterator_t index_lists = 
		new_index_list_setting_iterator(manager->combination_table);
	index_list_setting_t index_list;
	for (initialize(index_lists,&index_list);
	     has_next_element(index_lists);
	     next_element(index_lists,&index_list))
		manager->all_arrays[index_list.index_list_id-1].type = 
			INDEX_LIST;
	free_iterator(index_lists);
	iterator_t matrix_blocks =
		new_matrix_block_settings_iterator(manager->combination_table);
	matrix_block_setting_t matrix_block;
	for (initialize(matrix_blocks,&matrix_block);
	     has_next_element(matrix_blocks);
	     next_element(matrix_blocks,&matrix_block))
		manager->all_arrays[matrix_block.matrix_block_id-1].type =
			MATRIX_BLOCK;
	free_iterator(matrix_blocks);
}

void initialize_multi_thread_environment(memory_manager_t manager)
{
#pragma omp single
	{
		size_t num_threads = omp_get_num_threads();
		manager->calculation_threads_positions =
			(size_t*)calloc(num_threads,sizeof(size_t));	
	}
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

		size_t last_thread = get_last_thread(manager);
		size_t *candidate_arrays = manager->candidate_arrays_workspace;
		size_t num_candidate_arrays = 0;
		size_t can_unload = 0;
		for (size_t i = 0; i < manager->num_arrays; i++)
		{
			array_t *candidate = manager->all_arrays+i;
			if ((candidate->primary_array != NULL || 
			     candidate->secondary_array != NULL) &&
			    candidate->needed_by_instruction < last_thread &&
			    candidate->in_use == 0)
			{
				log_entry("Can unload %lu of size %lu B\n",
				       i+1,candidate->size_array);
				candidate_arrays[num_candidate_arrays++] =
					candidate->array_id;
				can_unload += candidate->size_array;
			}
		}
		for (size_t i = 0; i < num_candidate_arrays; i++)
		{
			size_t size_array =
			       	get_array_size(manager, candidate_arrays[i]);
			log_entry("Unloading array %lu\n",
			       candidate_arrays[i]);
			unload_array(manager,candidate_arrays[i]);
			size_unloaded_memory += size_array;
			if (size_unloaded_memory > min_size_to_unload)
				break;
		}
		usleep(1);
	}
}

static
void load_needed_arrays(memory_manager_t manager,
			execution_instruction_t instruction)
{
	load_array(manager,instruction.vector_block_in);
	load_array(manager,instruction.vector_block_out);
	load_array(manager,instruction.matrix_element_file);
	if (instruction.neutron_index != no_index) 
		load_array(manager,instruction.neutron_index);
	if (instruction.proton_index != no_index)
		load_array(manager,instruction.proton_index);
}

static
void wait_til_array_is_loaded(memory_manager_t manager, size_t array_id)
{
	struct timespec t1,t2;
	clock_gettime(CLOCK_REALTIME,&t1);
	while (!is_array_loaded(manager,array_id))
		usleep(1);
	clock_gettime(CLOCK_REALTIME,&t2);
	double waiting_time = (t2.tv_sec - t1.tv_sec)*1e6 +
		(t2.tv_nsec-t1.tv_nsec)*1e-3;
	manager->num_waits++;
	manager->total_wating_time += waiting_time;
	manager->min_waiting_time = 
		min(manager->min_waiting_time,waiting_time);
	manager->max_waiting_time =
		max(manager->max_waiting_time,waiting_time);
}
static
int is_array_loaded(memory_manager_t manager,
		    size_t array_id)
{
	array_t *array = &manager->all_arrays[array_id - 1];
	omp_set_nest_lock(&array->loading_array_lock);
	int value = array->primary_array != NULL ||
		array->secondary_array != NULL;
	omp_unset_nest_lock(&array->loading_array_lock);
	return value;
}

static
void load_array(memory_manager_t manager,
		size_t array_id)
{
	array_t *array = &manager->all_arrays[array_id-1];
	if (!omp_test_nest_lock(&array->loading_array_lock))
		return;
	if (is_array_loaded(manager,array_id))
	{
		omp_unset_nest_lock(&array->loading_array_lock);
		return;
	}
	switch(array->type)
	{
	case VECTOR_BLOCK:
		log_entry("It is a vector block\n");
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
		log_entry("It is an index list\n");
		array->primary_array =
			(void*)
			new_index_list_from_id
			(manager->index_list_base_directory,
			 array_id);
		break;
	case MATRIX_BLOCK:
		log_entry("It is a matrix block\n");
		array->primary_array =
			(void*)
			new_matrix_block(array_id,
					 manager->matrix_base_directory);		
		break;
	default:
		error("Can't load unkown array\n");
	}	
	manager->size_current_loaded_memory+=array->size_array;
	omp_unset_nest_lock(&array->loading_array_lock);
}

static
void unload_array(memory_manager_t manager,
		  size_t array_id)
{
	array_t *array = &manager->all_arrays[array_id-1];
	switch(array->type)
	{
	case VECTOR_BLOCK:
		log_entry("Unloading vector %p\n",
			  array->primary_array);
		free_vector_block((vector_block_t)array->primary_array);
		save_vector_block_elements((vector_block_t)
					   array->secondary_array);
		log_entry("Unloading vector %p\n",
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
	array->primary_array = NULL;
	array->secondary_array = NULL;
	manager->size_current_loaded_memory-=array->size_array;
}

static
void set_all_needed_by(memory_manager_t manager,
		       execution_instruction_t instruction)
{
	set_needed_by(manager,
		      instruction.vector_block_in,
		      instruction.instruction_index);
	set_needed_by(manager,
		      instruction.vector_block_out,
		      instruction.instruction_index);
	set_needed_by(manager,
		      instruction.matrix_element_file,
		      instruction.instruction_index);
	set_needed_by(manager,
		      instruction.neutron_index,
		      instruction.instruction_index);
	set_needed_by(manager,
		      instruction.proton_index,
		      instruction.instruction_index);
}

static
void set_needed_by(memory_manager_t manager,
		   size_t array_id,
		   size_t needed_by)
{
	if (array_id == no_index)
		return;
#pragma omp critical(set_needed_by)
	manager->all_arrays[array_id-1].needed_by_instruction = 
		max(manager->all_arrays[array_id-1].needed_by_instruction,
		    needed_by);
}

static
size_t get_last_thread(memory_manager_t manager)
{
	size_t last_thread = get_num_instructions(manager->execution_order);
	size_t num_threads = omp_get_num_threads();	
	omp_set_lock(&manager->calculation_threads_positions_lock);
	for (size_t i = 0; i < num_threads; i++)
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
