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
	omp_lock_t array_loading_in_other_thread_lock;
	// To make use no to threads use in_use at the same time
	omp_lock_t in_use_lock;
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
	evaluation_order_t evaluation_order;
	size_t size_current_loaded_memory;
	size_t maximum_loaded_memory;
	size_t *candidate_arrays_workspace;
	double total_wating_time;
	size_t num_waits;
	double min_waiting_time;
	double max_waiting_time;
	omp_lock_t size_current_loaded_memory_lock;
};

static
void initialize_arrays(memory_manager_t manager);


static
size_t get_size_of_unloaded_arrays(memory_manager_t manager,
				   evaluation_instruction_t instruction);

static
void make_space_for(memory_manager_t manager,
		    evaluation_instruction_t instruction);

static
void load_needed_arrays(memory_manager_t manager,
			evaluation_instruction_t instruction);
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
		       evaluation_instruction_t instruction);

static
void set_all_in_use(memory_manager_t manager,
		    evaluation_instruction_t instruction);

static
void set_needed_by(memory_manager_t manager,
		   size_t array_id,
		   size_t needed_by);

static
void set_in_use(memory_manager_t manager,
		size_t array_id);

static
int is_array_in_use(memory_manager_t manager,
		    size_t array_id);

static
size_t get_array_size(memory_manager_t manager,
		      size_t array_id);

static
uint64_t needed_by_key_function(size_t *array_index, memory_manager_t manager);

memory_manager_t new_memory_manager(const char *input_vector_base_directory,
				    const char *output_vector_base_directory,
				    const char *index_list_base_directory,
				    const char *matrix_base_directory,
				    combination_table_t combination_table,
				    evaluation_order_t evaluation_order,
				    size_t maximum_loaded_memory)
{
	memory_manager_t manager =
	       	(memory_manager_t)calloc(1,sizeof(struct _memory_manager_));
	manager->num_arrays = get_num_arrays(combination_table);
	manager->all_arrays = 
		(array_t*)calloc(manager->num_arrays, sizeof(array_t));
	manager->candidate_arrays_workspace =
	       	(size_t*)calloc(manager->num_arrays,sizeof(size_t));
	manager->input_vector_base_directory = copy_string(input_vector_base_directory);
	manager->output_vector_base_directory = copy_string(output_vector_base_directory);
	manager->index_list_base_directory = copy_string(index_list_base_directory);
	manager->matrix_base_directory = copy_string(matrix_base_directory);
	manager->combination_table = combination_table;
	manager->evaluation_order = evaluation_order;
	manager->maximum_loaded_memory = maximum_loaded_memory;
	manager->size_current_loaded_memory = 0;
	manager->num_waits = 0;
	manager->total_wating_time = 0.0;
	manager->min_waiting_time = INFINITY;
	manager->max_waiting_time = 0.0;
	omp_init_lock(&manager->size_current_loaded_memory_lock);
	initialize_arrays(manager);
	return manager;
}

void begin_instruction(memory_manager_t manager,
		       evaluation_instruction_t instruction)
{
	set_all_needed_by(manager,instruction);
#pragma omp critical (unloading)
	{
		set_all_in_use(manager,instruction);
		make_space_for(manager,instruction);
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
	return output_index_list;
}

void release_input_vector(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	omp_set_lock(&current_array->in_use_lock);
	current_array->in_use--;
	omp_unset_lock(&current_array->in_use_lock);
}

void release_output_vector(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	omp_set_lock(&current_array->in_use_lock);
	current_array->in_use--;
	omp_unset_lock(&current_array->in_use_lock);
}

void release_index_list(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	omp_set_lock(&current_array->in_use_lock);
	current_array->in_use--;
	omp_unset_lock(&current_array->in_use_lock);
}

void release_matrix_block(memory_manager_t manager, size_t array_id)
{
	array_t *current_array = &manager->all_arrays[array_id-1];
	omp_set_lock(&current_array->in_use_lock);
	current_array->in_use--;
	omp_unset_lock(&current_array->in_use_lock);
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
		array_t *array = &manager->all_arrays[i];
		omp_destroy_nest_lock(&array->loading_array_lock);
		omp_destroy_lock(&array->array_loading_in_other_thread_lock);
		omp_destroy_lock(&array->in_use_lock);
	}
	free(manager->all_arrays);
	free(manager->input_vector_base_directory);
	free(manager->output_vector_base_directory);
	free(manager->index_list_base_directory);
	free(manager->matrix_base_directory);
	free(manager->candidate_arrays_workspace);
	omp_destroy_lock(&manager->size_current_loaded_memory_lock);
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
		omp_init_lock(&current_array->array_loading_in_other_thread_lock);
		omp_init_lock(&current_array->in_use_lock);
	}	
	free(array_sizes);
	size_t num_threads = 0;
#pragma omp parallel
	num_threads = omp_get_num_threads();
	iterator_t basis_blocks =
	       	new_basis_block_iterator(manager->combination_table);
	basis_block_t current_basis_block;
	for (initialize(basis_blocks,&current_basis_block);
	     has_next_element(basis_blocks);
	     next_element(basis_blocks,&current_basis_block))
	{
		log_entry("Array %lu is a vector block\n",
		       current_basis_block.block_id);
		array_t *current_array =
			&manager->all_arrays[current_basis_block.block_id-1];
		current_array->type = 
			VECTOR_BLOCK;
		// Since there is one output vector per thread
		// and one input vector
		current_array->size_array *= (num_threads+1);
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

static
size_t get_size_of_unloaded_arrays(memory_manager_t manager,
				   evaluation_instruction_t instruction)
{
	size_t size_of_unloaded_arrays = 0;
	if (!is_array_loaded(manager,instruction.vector_block_in))
		size_of_unloaded_arrays += 
			get_array_size(manager,instruction.vector_block_in);
	if (instruction.vector_block_in != 
	    instruction.vector_block_out &&
	    !is_array_loaded(manager,instruction.vector_block_out))
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
void make_space_for(memory_manager_t manager,
		    evaluation_instruction_t instruction)
{
	size_t needed_memory = 
		get_size_of_unloaded_arrays(manager,instruction);
	if (needed_memory > manager->maximum_loaded_memory)
		error("Block %lu needs %lu B to be loaded,"
		      "But maximaly allowed loaded memory is %lu.\n",
		      instruction.instruction_index,
		      needed_memory,
		      manager->maximum_loaded_memory);


	if (manager->size_current_loaded_memory + needed_memory <
	    manager->maximum_loaded_memory)
		return;
	size_t needed_memory_to_unload = 
		(manager->size_current_loaded_memory + needed_memory) -
		manager->maximum_loaded_memory;
	size_t *candidates = manager->candidate_arrays_workspace;
	size_t num_candidates = 0;
	size_t can_unload = 0;
	do{
		can_unload = 0;	
		num_candidates = 0;
		size_t num_arrays_in_use = 0;
		for (size_t i = 0; i < manager->num_arrays; i++)
		{
			if (is_array_loaded(manager,i+1) &&
			    !is_array_in_use(manager,i+1))
			{
				can_unload += get_array_size(manager,i+1);
				candidates[num_candidates++] = i+1;
			}		
			num_arrays_in_use += is_array_in_use(manager,i+1);
		}
	}
	while (can_unload < needed_memory_to_unload);
	// Sort in order of ascending needed_by key word
	// to minimize risk that other threads 
#pragma omp critical(set_needed_by)
	rsort_r(candidates,
		num_candidates,
		sizeof(size_t),
		(__key_function_r_t)needed_by_key_function,
		manager);
	size_t unloaded_memory = 0;
	for (size_t i = 0; i<num_candidates; i++)
	{
		unload_array(manager,candidates[i]);
		unloaded_memory += get_array_size(manager,candidates[i]);
		if (unloaded_memory > needed_memory_to_unload)
			break;
	}
}

static
void load_needed_arrays(memory_manager_t manager,
			evaluation_instruction_t instruction)
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
	if (array_id == no_index)
		return 0;
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
	if (!omp_test_lock(&array->array_loading_in_other_thread_lock))
		return;
	omp_set_nest_lock(&array->loading_array_lock);
	if (is_array_loaded(manager,array_id))
	{
		omp_unset_nest_lock(&array->loading_array_lock);
		omp_unset_lock(&array->array_loading_in_other_thread_lock);
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
		error("Can't load unknown array\n");
	}	
	omp_set_lock(&manager->size_current_loaded_memory_lock);
	manager->size_current_loaded_memory+=array->size_array;
	omp_unset_lock(&manager->size_current_loaded_memory_lock);
	omp_unset_nest_lock(&array->loading_array_lock);
	omp_unset_lock(&array->array_loading_in_other_thread_lock);
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
		error("Can't load unknown array\n");
	}	
	array->primary_array = NULL;
	array->secondary_array = NULL;
	omp_set_lock(&manager->size_current_loaded_memory_lock);
	manager->size_current_loaded_memory-=array->size_array;
	omp_unset_lock(&manager->size_current_loaded_memory_lock);
}

static
void set_all_needed_by(memory_manager_t manager,
		       evaluation_instruction_t instruction)
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
void set_all_in_use(memory_manager_t manager,
		    evaluation_instruction_t instruction)
{
	set_in_use(manager,instruction.vector_block_in);
	set_in_use(manager,instruction.vector_block_out);
	set_in_use(manager,instruction.matrix_element_file);
	set_in_use(manager,instruction.neutron_index);
	set_in_use(manager,instruction.proton_index);
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
void set_in_use(memory_manager_t manager,
		size_t array_id)
{
	if (array_id == no_index)
		return;
	array_t *array = &manager->all_arrays[array_id - 1];
	omp_set_lock(&array->in_use_lock);
	array->in_use++;
	omp_unset_lock(&array->in_use_lock);
}

static
int is_array_in_use(memory_manager_t manager,
		    size_t array_id)
{
	if (array_id == no_index)
		return 0;	
	array_t *array = &manager->all_arrays[array_id - 1];
	omp_set_lock(&array->in_use_lock);
	int state = array->in_use > 0;
	omp_unset_lock(&array->in_use_lock);
	return state;
}

static
size_t get_array_size(memory_manager_t manager,
		      size_t array_id)
{
	return manager->all_arrays[array_id-1].size_array;	
}

static
uint64_t needed_by_key_function(size_t *array_id, memory_manager_t manager)
{
	return manager->all_arrays[(*array_id)-1].needed_by_instruction;
}
