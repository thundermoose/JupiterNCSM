#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <arguments/arguments.h>
#include <input/read_2nf_antoine_format.h>
#include <transform_block_settings/transform_block_settings.h>
#include <transform_2nf_block_manager/transform_2nf_block_manager.h>
#include <transform_3nf_block_manager/transform_3nf_block_manager.h>
#include <combination_table/combination_table.h>
#include <iterator/iterator.h>
#include <matrix_energy_block/matrix_energy_block.h>
#include <radix_sort/radix_sort.h>
#include <error/error.h>
#include <log/log.h>
#include <time.h>
#include <omp.h>

	__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MERCURY_LOGFILE",
			 "mercury_J_scheme_to_internal.log");	 
}

static
void generate_1nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

static
void generate_2nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

static
void generate_2nf_zero_matrix_blocks(combination_table_t combination_table,
				     arguments_t arguments);


void generate_3nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments);

static
void generate_3nf_matrix_blocks_parallel(combination_table_t combination_table,
					 arguments_t arguments);

static
void process_matrix_energy_block(matrix_energy_block_t current_block,
				 transform_3nf_block_manager_t manager,
				 const char *output_path_base,
				 size_t block_index,
				 FILE *finished_block_file);

static
void mark_block_as_finished(FILE* finished_block_file,size_t block_index);

static
void load_finished_blocks(size_t **finished_blocks,
			  size_t *num_finished_blocks,
			  const char *finished_energy_blocks_path);

static
uint64_t size_t_key(const size_t *element);

int main(int num_arguments,
	 char **argument_list)
{
	struct timespec t_start,t_end;
	printf("Mercury starts:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	arguments_t arguments =
	       	parse_argument_list(num_arguments, argument_list);	    
	if (to_few_arguments(arguments))
	{
		show_usage(arguments);
		free_arguments(arguments);
		return EXIT_SUCCESS;
	}
	combination_table_t combination_table =
		new_combination_table
		(get_combination_file_path_argument(arguments),
		 get_num_protons_argument(arguments),
		 get_num_neutrons_argument(arguments));
	generate_1nf_matrix_blocks(combination_table,
				   arguments);
	if (no_2nf_argument(arguments))
		generate_2nf_zero_matrix_blocks(combination_table, arguments);
	else
		generate_2nf_matrix_blocks(combination_table, arguments);
	generate_3nf_matrix_blocks_parallel(combination_table, arguments);
	free_combination_table(combination_table);
	free_arguments(arguments);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double mercury_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Mercury ends after %lg µs\n",mercury_time);
	return EXIT_SUCCESS;
}

static
void generate_1nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	struct timespec t_start,t_end;
	printf("Generate 1nf matrix blocks:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	const char *output_path_base = get_output_path_argument(arguments);
	const char *index_list_path = get_index_list_path_argument(arguments);
	while (has_next_1nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_1nf_block_iterator(combination_table);
		connection_list_t connection_list =
			read_connection_files(index_list_path,
					      current_matrix_block);
		mercury_matrix_block_t current_block =
			new_zero_mercury_matrix_block(connection_list);
		save_mercury_matrix_block(current_block,
					  output_path_base);
		free_mercury_matrix_block(current_block);
		free_connection_list(connection_list);
	}
	clock_gettime(CLOCK_REALTIME,&t_end);
	double generate_1nf_matrix_blocks_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Geneerate 1nf matrix blocks ends after %lg µs\n",
	       generate_1nf_matrix_blocks_time);
}

static
void generate_2nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	struct timespec t_start,t_end;
	printf("Generate 2nf matrix blocks:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	log_entry("num_particle_argument = %lu",
		  get_num_particles_argument(arguments));
	antoine_2nf_file_t coupled_2nf_data =
		open_antoine_2nf_file
		(get_interaction_path_2nf_argument(arguments),
		 get_num_particles_argument(arguments),
		 get_single_particle_energy_argument(arguments),
		 get_two_particle_energy_argument(arguments));
	transform_2nf_block_manager_t manager =
		new_transform_2nf_block_manager
		(coupled_2nf_data,
		 get_index_list_path_argument(arguments),
		 get_index_list_path_argument(arguments),
		 get_single_particle_energy_argument(arguments));
	transform_block_settings_t transformed_block = {INT_MAX};	
	const char *output_path_base = get_output_path_argument(arguments);
	while (has_next_2nf_block(combination_table))
	{
		matrix_block_setting_t current_matrix_block = 
			next_2nf_block_iterator(combination_table);
		transform_block_settings_t current_block =
			setup_transform_block(current_matrix_block);
		log_entry("<(pE %d nE %d)|(pE %d nE %d)> Tz = %d",
			  current_block.proton_energy_bra,
			  current_block.neutron_energy_bra,
			  current_block.proton_energy_ket,
			  current_block.neutron_energy_ket,
			  current_block.total_isospin);
		if (compare_transform_block_settings(&transformed_block,
						     &current_block) != 0)
		{
			log_entry("Needs to decouple a new block");
			decouple_transform_2nf_block(manager,
						     current_block);
			transformed_block = current_block;
		}
		mercury_matrix_block_t matrix_block = 
			get_transform_2nf_matrix_block(manager,
						       current_matrix_block);
		log_entry("Retrieved the current matrix_block");
		save_mercury_matrix_block(matrix_block,
					  output_path_base);
		log_entry("Saved the current matrix_block to file");
		free_mercury_matrix_block(matrix_block);
	}
	free_transform_2nf_block_manager(manager);
	free_antoine_2nf_file(coupled_2nf_data);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double generate_2nf_matrix_blocks_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Geneerate 2nf matrix blocks ends after %lg µs\n",
	       generate_2nf_matrix_blocks_time);
}

static
void generate_2nf_zero_matrix_blocks(combination_table_t combination_table,
				     arguments_t arguments)
{
	struct timespec t_start,t_end;
	printf("Generate 2nf matrix blocks:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	const char *output_path_base = get_output_path_argument(arguments);
	const char *index_list_path = get_index_list_path_argument(arguments);
	iterator_t iterator_2nf_blocks = 
		new_2nf_matrix_block_settings_iterator(combination_table);
	matrix_block_setting_t current_matrix_block;
	for (initialize(iterator_2nf_blocks,&current_matrix_block);
	     has_next_element(iterator_2nf_blocks);
	     next_element(iterator_2nf_blocks,&current_matrix_block))
	{
		connection_list_t connection_list =
			read_connection_files(index_list_path,
					      current_matrix_block);
		mercury_matrix_block_t current_block =
			new_zero_mercury_matrix_block(connection_list);
		save_mercury_matrix_block(current_block,
					  output_path_base);
		free_mercury_matrix_block(current_block);
		free_connection_list(connection_list);
	}
	clock_gettime(CLOCK_REALTIME,&t_end);
	double generate_2nf_matrix_blocks_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Geneerate 2nf matrix blocks ends after %lg µs\n",
	       generate_2nf_matrix_blocks_time);
}

void generate_3nf_matrix_blocks(combination_table_t combination_table,
				arguments_t arguments)
{
	if (get_interaction_path_3nf_argument(arguments) == NULL)
		return;
	struct timespec t_start,t_end;
	printf("Generate 3nf matrix blocks:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	printf("3NF blocks only:\n");
	Data_File *coupled_3nf_data =
		open_data_file(get_interaction_path_3nf_argument(arguments));
	transform_3nf_block_manager_t manager =
	       	new_transform_3nf_block_manager
		(coupled_3nf_data,
		 get_index_list_path_argument(arguments),
		 get_single_particle_energy_argument(arguments));
	transform_block_settings_t transformed_block = {INT_MAX};
	const char *output_path_base = get_output_path_argument(arguments);					
	struct timespec time_start;
	struct timespec time_end;
	while (has_next_3nf_block(combination_table))
	{
		clock_gettime(CLOCK_REALTIME,&time_start);
		matrix_block_setting_t current_matrix_block = 
			next_3nf_block_iterator(combination_table);
		printf("Compute 3NF block %lu\n",
		       current_matrix_block.matrix_block_id);
		transform_block_settings_t current_block =
			setup_transform_block(current_matrix_block);
		log_entry("%s p: %d %d %d n: %d %d %d,"
			  "#PC: %lu #NC: %lu id: %lu\n",
			  block_type_to_string(current_matrix_block.type),
			  current_matrix_block.difference_energy_protons,
			  current_matrix_block.difference_M_protons,
			  current_matrix_block.depth_protons,
			  current_matrix_block.difference_energy_neutrons,
			  current_matrix_block.difference_M_neutrons,
			  current_matrix_block.depth_neutrons,
			  current_matrix_block.num_proton_combinations,
			  current_matrix_block.num_neutron_combinations,
			  current_matrix_block.matrix_block_id);
		if (compare_transform_block_settings(&current_block,
						     &transformed_block) != 0)
		{
			decouple_transform_3nf_block(manager,
						     current_block);
			transformed_block = current_block;
		}
		mercury_matrix_block_t matrix_block =
			get_transform_3nf_matrix_block(manager,
						       current_matrix_block);
		save_mercury_matrix_block(matrix_block,
					  output_path_base);
		free_mercury_matrix_block(matrix_block);
		clock_gettime(CLOCK_REALTIME,&time_end);
		double elapsed_time = (time_end.tv_sec-time_start.tv_sec)*1e6 +
			(time_end.tv_nsec-time_start.tv_nsec)*1e-3;
		printf("Time: %lg µs\n",elapsed_time);
	}
	free_transform_3nf_block_manager(manager);
	free_data_file(coupled_3nf_data);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double generate_3nf_matrix_blocks_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Geneerate 3nf matrix blocks ends after %lg µs\n",
	       generate_3nf_matrix_blocks_time);
}

static
void generate_3nf_matrix_blocks_parallel(combination_table_t combination_table,
					 arguments_t arguments)
{
	if (get_interaction_path_3nf_argument(arguments) == NULL)
		return;
	struct timespec t_start,t_end;
	printf("Generate 3nf matrix blocks:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	printf("3NF blocks only:\n");
	Data_File *coupled_3nf_data =
		open_data_file(get_interaction_path_3nf_argument(arguments));
	set_max_loaded_memory(coupled_3nf_data,
			      get_max_loaded_memory_argument(arguments));
	transform_3nf_block_manager_t manager =
		new_transform_3nf_block_manager
		(coupled_3nf_data,
		 get_index_list_path_argument(arguments),
		 get_single_particle_energy_argument(arguments));
	const char *output_path_base = get_output_path_argument(arguments);					
	size_t current_block_index = 0;
	size_t *finished_blocks = NULL;
	size_t num_finished_blocks = 0;
	load_finished_blocks(&finished_blocks,
			     &num_finished_blocks,
			     get_finished_energy_blocks_argument(arguments));
	size_t current_finished_block = 0;
	FILE *finished_block_file =
	       	fopen(get_finished_energy_blocks_argument(arguments),"a");
#pragma omp parallel
	{
#pragma omp single
		{
			while (has_next_3nf_matrix_energy_block(combination_table))
			{
				matrix_energy_block_t current_energy_block =
					next_3nf_matrix_energy_block(combination_table);
				if (current_finished_block < num_finished_blocks && 
				    current_block_index == finished_blocks[current_finished_block])
				{
					current_finished_block++;	
				}
				else
				{
#pragma omp task
					process_matrix_energy_block
						(current_energy_block,
						 manager,
						 output_path_base,
						 current_block_index,
						 finished_block_file);
				}
				current_block_index++;
			}
		}
	}
	fclose(finished_block_file);
	free_transform_3nf_block_manager(manager);
	free_data_file(coupled_3nf_data);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double generate_3nf_matrix_blocks_time = 
		(t_end.tv_sec - t_start.tv_sec)*1e6 +
		(t_end.tv_nsec - t_start.tv_nsec)*1e-3;
	printf("Geneerate 3nf matrix blocks ends after %lg µs\n",
	       generate_3nf_matrix_blocks_time);
}

static
void process_matrix_energy_block(matrix_energy_block_t current_block,
				 transform_3nf_block_manager_t manager,
				 const char *output_path_base,
				 size_t block_index,
				 FILE *finished_block_file)
{
	printf("Launching block %lu\n",block_index);
	transformed_block_t current_transformed_block = 
		get_transformed_block(manager,current_block);
	while (has_next_energy_matrix_block(current_block))
	{
		matrix_block_setting_t current_matrix_block_settings =
			next_energy_matrix_block(current_block);
		mercury_matrix_block_t matrix_block =
			get_3nf_mercury_matrix(current_transformed_block,
					       current_matrix_block_settings);
		save_mercury_matrix_block(matrix_block,output_path_base);
		free_mercury_matrix_block(matrix_block);
	}	
	free_transformed_block(current_transformed_block);
	free_matrix_energy_block(current_block);
	mark_block_as_finished(finished_block_file,block_index);	
	printf("Finished block %lu\n",block_index);
}

static
void mark_block_as_finished(FILE* finished_block_file,size_t block_index)
{
#pragma omp critical (mark_block_as_finished)
	{
		fwrite(&block_index,sizeof(size_t),1,finished_block_file);
		fflush(finished_block_file);
	}
}

static
void load_finished_blocks(size_t **finished_blocks,
			  size_t *num_finished_blocks,
			  const char *finished_energy_blocks_path)
{
	FILE *file = fopen(finished_energy_blocks_path,"r");
	if (file == NULL)
		return; // File don't exists so all blocks should be done
	fseek(file,0,SEEK_END);
	size_t num_bytes = ftell(file);
	fseek(file,0,SEEK_SET);
	*num_finished_blocks = num_bytes/sizeof(size_t);
	size_t *loaded_blocks = (size_t*)malloc(num_bytes);	
	if (fread(loaded_blocks,sizeof(size_t),*num_finished_blocks,file) <
	    *num_finished_blocks)
		error("Could not read %lu bytes from %s.\n",
      			num_bytes,
			finished_energy_blocks_path);
	fclose(file);	
	// Sorting since the blocks may not be in order
	rsort(loaded_blocks,
	      *num_finished_blocks,
	      sizeof(size_t),
	      (__key_function_t)size_t_key);
	*finished_blocks = loaded_blocks;
}

static
uint64_t size_t_key(const size_t *element)
{
	return *element;
}
