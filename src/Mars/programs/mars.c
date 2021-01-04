#include <stdlib.h>
#include <stdio.h>
#include <combination_table/combination_table.h>
#include <settings/settings.h>
#include <bucket_sort/bucket_sort.h>
#include <radix_sort/radix_sort.h>
#include <log/log.h>
#include <global_constants/global_constants.h>

__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MARS_LOGFILE",
			 "mars.log");
}

static
calculation_blocks_t 
optimize_calculation_blocks(calculation_blocks_t calculation_blocks,
			    combination_table_t combination_table,
			    settings_t settings);

static
void save_evaluation_order_file(calculation_blocks_t calculation_blocks,
			       const char *output_path);

static
uint64_t num_multiplications(calculation_block_t *block,
			     size_t *array_sizes);

static
uint64_t smallest_block_array(calculation_block_t *block,
			      size_t *array_sizes);

static
uint64_t second_smallest_block_array(calculation_block_t *block,
				     size_t *array_sizes);


static
uint64_t third_smallest_block_array(calculation_block_t *block,
				    size_t *array_sizes);

static
uint64_t second_largest_block_array(calculation_block_t *block,
				    size_t *array_sizes);

static
uint64_t largest_block_array(calculation_block_t *block,
			     size_t *array_sizes);

static
uint64_t largest_block_array_size(calculation_block_t *block,
				  size_t *array_sizes);

static
void quint_sort(size_t array[5], size_t *keys);

static
size_t find_force_block(calculation_block_t *blocks,
			size_t num_particles,
			size_t num_blocks);

int main(int num_arguments,
	 char **argument_list)
{
	settings_t settings = parse_settings(num_arguments,argument_list);
	if (should_display_usage(settings))
	{
		display_usage(settings);
		return EXIT_FAILURE;
	}
	combination_table_t combination_table =
	       	new_combination_table(get_combination_table_path(settings),
				      get_num_protons_argument(settings),
				      get_num_neutrons_argument(settings));	
	calculation_blocks_t calculation_blocks =
	       	get_calculation_blocks(combination_table);
	calculation_blocks_t optimized_calculation_blocks =
		optimize_calculation_blocks(calculation_blocks,
					    combination_table,
					    settings);
	save_evaluation_order_file(optimized_calculation_blocks,
				  get_output_path(settings));
	// frees both combination_table and calculation_blocks
	free_combination_table(combination_table);
	free_calculation_blocks(optimized_calculation_blocks);
	free_settings(settings);
	return EXIT_SUCCESS;
}

static
calculation_blocks_t
optimize_calculation_blocks(calculation_blocks_t calculation_blocks,
			    combination_table_t combination_table,
			    settings_t settings)
{
	calculation_block_t *blocks =
		get_calculation_blocks_array(calculation_blocks);
	size_t num_calculation_blocks = 
		get_num_calculation_blocks(calculation_blocks);
	bucket_sort(blocks,
		    num_calculation_blocks,
		    sizeof(calculation_block_t),
		    3,
		    (bucket_index_function_t)
		    get_calculation_block_num_particles);
	// since all the 1nf matrix elements are 0 we skip all such blocks
	size_t two_particle_forces = 
		find_force_block(blocks,2,num_calculation_blocks);
	size_t three_particle_forces =
	       	find_force_block(blocks,3,num_calculation_blocks);
	calculation_block_t *needed_blocks = blocks+two_particle_forces;
	size_t num_needed_calculation_blocks = 
		two_nucleon_force_only_mode(settings) ?
		three_particle_forces - two_particle_forces :
		num_calculation_blocks - two_particle_forces;
	size_t *array_sizes = get_array_sizes(combination_table);

			
	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)smallest_block_array,
		array_sizes);


	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)second_smallest_block_array,
		array_sizes);

	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)third_smallest_block_array,
		array_sizes);

	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)second_largest_block_array,
		array_sizes);

	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)largest_block_array,
		array_sizes);

	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)largest_block_array_size,
		array_sizes);

	rsort_r(needed_blocks,
		num_needed_calculation_blocks,
		sizeof(calculation_block_t),
		(__key_function_r_t)num_multiplications,
		array_sizes);
	free(array_sizes);

	calculation_blocks_t optimized_calculation_blocks =
		new_calculation_blocks(needed_blocks,
				       num_needed_calculation_blocks);
	free(blocks);
	return optimized_calculation_blocks;
}

static
void save_evaluation_order_file(calculation_blocks_t calculation_blocks,
			       const char *output_path)
{
	FILE *evaluation_order_file = fopen(output_path,"w");
	reset_calculation_block_iterator(calculation_blocks);
	while (has_next_calculation_block(calculation_blocks))
	{
		calculation_block_t block = 
			next_calculation_block(calculation_blocks);
		if (block.secondary_index_list == no_index)
			fprintf(evaluation_order_file,
				"BLOCK: %lu %lu %lu %lu\n",
				block.vector_block_in,
				block.vector_block_out,
				block.primary_index_list,
				block.matrix_element_block);
		else
			fprintf(evaluation_order_file,
				"BLOCK: %lu %lu %lu %lu %lu\n",
				block.vector_block_in,
				block.vector_block_out,
				block.primary_index_list,
				block.secondary_index_list,
				block.matrix_element_block);
	}
	fclose(evaluation_order_file);
}

static
uint64_t num_multiplications(calculation_block_t *block,
			     size_t *array_sizes)
{
	return -block->num_multiplications;
}

static
uint64_t smallest_block_array(calculation_block_t *block,
			      size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return arrays[0];
}

static
uint64_t second_smallest_block_array(calculation_block_t *block,
				     size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return arrays[1];
}


static
uint64_t third_smallest_block_array(calculation_block_t *block,
				    size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return arrays[2];
}

static
uint64_t second_largest_block_array(calculation_block_t *block,
				    size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return arrays[3];
}

static
uint64_t largest_block_array(calculation_block_t *block,
			     size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return arrays[4];
}

static
uint64_t largest_block_array_size(calculation_block_t *block,
				  size_t *array_sizes)
{
	size_t arrays[5] =
	{
		block->vector_block_in,
		block->vector_block_out,
		block->primary_index_list,
		block->secondary_index_list,
		block->matrix_element_block
	};	
	quint_sort(arrays,array_sizes);
	return array_sizes[arrays[4]-1];
}

static
void quint_sort(size_t array[5], size_t *keys)
{
	// simply a buble sort
	for (size_t i = 1; i<5; i++)
		for (size_t j = 0; j<i; j++)
		{
			size_t j_size = 
				array[j] == no_index ? 
				0 : keys[array[j]-1];
			size_t jn_size =
				array[j+1] == no_index ?
				0 : keys[array[j+1]-1];
			if (j_size > jn_size)
			{
				size_t tmp = array[j];
				array[j] = array[j+1];
				array[j+1] = tmp;
			}
		}
}

static
size_t find_force_block(calculation_block_t *blocks,
			size_t num_particles,
			size_t num_blocks)
{
	for (size_t i = 0; i<num_blocks; i++)
		if (blocks[i].num_particles == num_particles)
			return i;
	return no_index;
}
