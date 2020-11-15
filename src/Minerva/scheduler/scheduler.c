#include <scheduler/scheduler.h>
#include <memory_manager/memory_manager.h>
#include <matrix_vector_multiplication/matrix_vector_multiplication.h>
#include <string_tools/string_tools.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <error/error.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

struct _scheduler_
{
	execution_order_t execution_order;
	char *index_lists_base_directory;
	char *matrix_file_base_directory;
	combination_table_t combination_table;
	size_t maximum_loaded_memory;
};

static
void execute_instruction(execution_instruction_t instruction,
			 memory_manager_t memory_manager,
			 scheduler_t scheduler);

static
void neutron_case(memory_manager_t memory_manager,
		  execution_instruction_t instruction);

static
void proton_case(memory_manager_t memory_manager,
		 execution_instruction_t instruction);
static
void neutron_proton_case(memory_manager_t memory_manager,
			 execution_instruction_t instruction);

static
void diagonal_neutron_case(memory_manager_t memory_manager,
			   execution_instruction_t instruction);

static
void diagonal_proton_case(memory_manager_t memory_manager,
			  execution_instruction_t instruction);
static
void diagonal_neutron_proton_case(memory_manager_t memory_manager,
				  execution_instruction_t instruction);

static
void off_diagonal_neutron_case(memory_manager_t memory_manager,
			       execution_instruction_t instruction);

static
void off_diagonal_proton_case(memory_manager_t memory_manager,
			      execution_instruction_t instruction);
static
void off_diagonal_neutron_proton_case(memory_manager_t memory_manager,
				      execution_instruction_t instruction);

scheduler_t new_scheduler(execution_order_t execution_order,
			  combination_table_t combination_table,
			  const char *index_lists_base_directory,
			  const char *matrix_file_base_directory,
			  size_t maximum_loaded_memory)
{
	scheduler_t scheduler =
		(scheduler_t)malloc(sizeof(struct _scheduler_));
	scheduler->execution_order = execution_order;
	scheduler->combination_table = combination_table;
	scheduler->index_lists_base_directory =
		copy_string(index_lists_base_directory);
	scheduler->matrix_file_base_directory =
		copy_string(matrix_file_base_directory);
	scheduler->maximum_loaded_memory = maximum_loaded_memory;
	return scheduler;
}

void run_matrix_vector_multiplication(const char *output_vector_base_directory,
				      const char *input_vector_base_directory,
				      scheduler_t scheduler)
{
	memory_manager_t memory_manager = 
		new_memory_manager(input_vector_base_directory,
				   output_vector_base_directory,
				   scheduler->index_lists_base_directory,
				   scheduler->matrix_file_base_directory,
				   scheduler->combination_table,
				   scheduler->execution_order,
				   scheduler->maximum_loaded_memory);
	execution_order_iterator_t instruction_iterator =
		get_execution_order_iterator(scheduler->execution_order);
	double fastest_block_time = INFINITY;
	double slowes_block_time = -INFINITY;
	double total_block_time = 0;
#pragma omp parallel shared(memory_manager,scheduler,instruction_iterator)
	{
#pragma omp single
		{
			printf("There are %d threads\n",
			       omp_get_num_threads());
		}
		size_t thread_id = omp_get_thread_num();
		printf("thread_id = %lu\n",thread_id);
		//if (thread_id < 4)
		//	launch_memory_manager_thread(memory_manager);
		//else
		while (has_next_instruction(instruction_iterator))
		{
			execution_instruction_t instruction =
				next_instruction(instruction_iterator);
			begin_instruction(memory_manager,instruction);
			struct timespec t_start,t_end;
			clock_gettime(CLOCK_REALTIME,&t_start);
			execute_instruction(instruction,
					    memory_manager,
					    scheduler);
			clock_gettime(CLOCK_REALTIME,&t_end);
			double block_time = 
				(t_end.tv_sec - t_start.tv_sec)*1e6+
				(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
			if (instruction.type != unload)
			{
#pragma omp critical
				fastest_block_time = 
					min(fastest_block_time,
					    block_time);
#pragma omp critical
				slowes_block_time= 
					max(slowes_block_time,
					    block_time);
#pragma omp critical
				total_block_time += block_time;
			}
		}
	}
	free_memory_manager(memory_manager);
	free_execution_order_iterator(instruction_iterator);
	printf("Fastest block: %lg µs\n",fastest_block_time);
	printf("Slowest block: %lg µs\n",slowes_block_time);
	printf("Average block: %lg µs\n",
	       total_block_time /
	       get_num_instructions(scheduler->execution_order));
}

void free_scheduler(scheduler_t scheduler)
{
	free(scheduler->index_lists_base_directory);
	free(scheduler->matrix_file_base_directory);
	free(scheduler);
}

	static
void execute_instruction(execution_instruction_t instruction,
			 memory_manager_t memory_manager,
			 scheduler_t scheduler)
{
	log_entry("instruction = \n"
		  "{\n"
		  "\t.type = %d,\n"
		  "\t.vector_block_in = %lu,\n"
		  "\t.vector_block_out = %lu,\n"
		  "\t.matrix_element_file = %lu,\n"
		  "\t.neutron_index = %lu,\n"
		  "\t.proton_index = %lu\n"
		  "}\n",
		  instruction.type,
		  instruction.vector_block_in,
		  instruction.vector_block_out,
		  instruction.matrix_element_file,
		  instruction.neutron_index,
		  instruction.proton_index);
	switch (instruction.type)
	{
		case neutron_block:
			neutron_case(memory_manager, instruction);
			break;
		case proton_block:
			proton_case(memory_manager, instruction);
			break;
		case neutron_proton_block:
			neutron_proton_case(memory_manager, instruction);
			break;
		case unload:
			log_entry("Ignoring unloading\n");
			break;
		default:
			error("instruction type is not known\n");
			break;
	}	
}

static
void neutron_case(memory_manager_t memory_manager,
		  execution_instruction_t instruction)
{
	if (instruction.vector_block_in == instruction.vector_block_out)
		diagonal_neutron_case(memory_manager,instruction);
	else
		off_diagonal_neutron_case(memory_manager,instruction);
}

static
void proton_case(memory_manager_t memory_manager,
		 execution_instruction_t instruction)
{
	if (instruction.vector_block_in == instruction.vector_block_out)
		diagonal_proton_case(memory_manager,instruction);
	else
		off_diagonal_proton_case(memory_manager,instruction);
}

static
void neutron_proton_case(memory_manager_t memory_manager,
			 execution_instruction_t instruction)
{
	if (instruction.vector_block_in == instruction.vector_block_out)
		diagonal_neutron_proton_case(memory_manager,instruction);
	else
		off_diagonal_neutron_proton_case(memory_manager,instruction);
}

	static
void diagonal_neutron_case(memory_manager_t memory_manager,
			   execution_instruction_t instruction)
{
	log_entry("Running the neutron only case");
	vector_block_t input_vector_block =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t output_vector_block =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t list = 
		request_index_list(memory_manager,
				instruction.neutron_index);
	multiplication_neutrons(output_vector_block,
				input_vector_block,
				matrix_block,
				list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.neutron_index);
}

	static
void diagonal_proton_case(memory_manager_t memory_manager,
			  execution_instruction_t instruction)
{
	log_entry("Running the proton only case");
	vector_block_t input_vector_block =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t output_vector_block =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t list = 
		request_index_list(memory_manager,
				instruction.proton_index);
	multiplication_protons(output_vector_block,
			       input_vector_block,
			       matrix_block,
			       list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.proton_index);
}

	static
void diagonal_neutron_proton_case(memory_manager_t memory_manager,
				  execution_instruction_t instruction)
{
	log_entry("Running the neutron and proton case");
	vector_block_t input_vector_block =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t output_vector_block =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t  neutron_list =
		request_index_list(memory_manager,
				instruction.neutron_index);
	index_list_t  proton_list =
		request_index_list(memory_manager,
				instruction.proton_index);
	multiplication_neutrons_protons(output_vector_block,
					input_vector_block,
					matrix_block,
					neutron_list,
					proton_list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.neutron_index);
	release_index_list(memory_manager,instruction.proton_index);
}

	static
void off_diagonal_neutron_case(memory_manager_t memory_manager,
			   execution_instruction_t instruction)
{
	log_entry("Running the neutron only case");
	vector_block_t input_vector_block_left =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t input_vector_block_right =
		request_input_vector_block(memory_manager,
					instruction.vector_block_out);
	vector_block_t output_vector_block_left =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	vector_block_t output_vector_block_right =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_in);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t list = 
		request_index_list(memory_manager,
				instruction.neutron_index);
	multiplication_neutrons_off_diag(output_vector_block_left,
					 output_vector_block_right,
					 input_vector_block_left,
					 input_vector_block_right,
					 matrix_block,
					 list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_input_vector(memory_manager,instruction.vector_block_out);
	release_output_vector(memory_manager,instruction.vector_block_in);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.neutron_index);
}

	static
void off_diagonal_proton_case(memory_manager_t memory_manager,
			  execution_instruction_t instruction)
{
	log_entry("Running the proton only case");
	vector_block_t input_vector_block_left =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t input_vector_block_right =
		request_input_vector_block(memory_manager,
					instruction.vector_block_out);	
	vector_block_t output_vector_block_left =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	vector_block_t output_vector_block_right =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_in);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t list = 
		request_index_list(memory_manager,
				instruction.proton_index);
	multiplication_protons_off_diag(output_vector_block_left,
					output_vector_block_right,
					input_vector_block_left,
					input_vector_block_right,
					matrix_block,
					list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_input_vector(memory_manager,instruction.vector_block_out);
	release_output_vector(memory_manager,instruction.vector_block_in);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.proton_index);
}

	static
void off_diagonal_neutron_proton_case(memory_manager_t memory_manager,
				  execution_instruction_t instruction)
{
	log_entry("Running the neutron and proton case");
	vector_block_t input_vector_block_left =
		request_input_vector_block(memory_manager,
					instruction.vector_block_in);	
	vector_block_t input_vector_block_right =
		request_input_vector_block(memory_manager,
					instruction.vector_block_out);	
	vector_block_t output_vector_block_left =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_out);	
	vector_block_t output_vector_block_right =
		request_output_vector_block(memory_manager,
					 instruction.vector_block_in);	
	matrix_block_t matrix_block =
		request_matrix_block(memory_manager,
				  instruction.matrix_element_file);
	index_list_t  neutron_list =
		request_index_list(memory_manager,
				instruction.neutron_index);
	index_list_t  proton_list =
		request_index_list(memory_manager,
				instruction.proton_index);
	multiplication_neutrons_protons_off_diag(output_vector_block_left,
						 output_vector_block_right,
					input_vector_block_left,
					input_vector_block_right,
					matrix_block,
					neutron_list,
					proton_list);
	release_input_vector(memory_manager,instruction.vector_block_in);
	release_output_vector(memory_manager,instruction.vector_block_out);
	release_input_vector(memory_manager,instruction.vector_block_out);
	release_output_vector(memory_manager,instruction.vector_block_in);
	release_matrix_block(memory_manager,instruction.matrix_element_file);
	release_index_list(memory_manager,instruction.neutron_index);
	release_index_list(memory_manager,instruction.proton_index);
}
