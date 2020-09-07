#include <combination_table/combination_table.h>
#include <basis_block/basis_block.h>
#include <log/log.h>
#include <array_builder/array_builder.h>
#include <string_tools/string_tools.h>
#include <error/error.h>
#include <thundertester/test.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <debug_mode/debug_mode.h>

struct _combination_table_
{
	basis_block_t *basis_blocks;
	size_t num_basis_blocks;
	size_t current_basis_block_index;
	index_list_setting_t *index_list_settings;
	size_t num_index_list_settings;
	size_t current_index_list_setting;
	matrix_block_setting_t *matrix_block_settings;
	size_t num_matrix_block_settings;
	size_t current_matrix_block_index;
	size_t length_1nf_blocks;
	size_t iterator_index_1nf_blocks;
	size_t index_to_2nf_blocks;
	size_t length_2nf_blocks;
	size_t iterator_index_2nf_blocks;
	size_t index_to_3nf_blocks;
	size_t length_3nf_blocks;
	size_t iterator_index_3nf_blocks;
	size_t *id_to_index_map;
	size_t num_ids;
};

static
void read_basis_blocks(FILE *table_file,
		       array_builder_t basis_blocks_builder,
		       array_builder_t id_to_index_map_builder,
		       size_t num_protons,
		       size_t num_neutrons);

static
void read_index_list_settings(FILE *table_file,
			      array_builder_t index_list_settings_builder,
			      array_builder_t id_to_index_map_builder);

static
void read_matrix_block_settings(FILE *table_file,
				array_builder_t matrix_block_settings_builder,
				array_builder_t id_to_index_map_builder);

static void 
sort_matrix_blocks_on_num_particles(combination_table_t table,
				    array_builder_t id_to_index_map_builder);
static
int is_title_row(const char *row);

static
void move_to_title(FILE *table_file,
		   const char *title);

static
size_t interpret_id_string(const char *id_string);

static
int compare_matrix_block_settings(matrix_block_setting_t *block_a,
				  matrix_block_setting_t *block_b);

combination_table_t new_combination_table(const char *filename,
					  size_t num_protons,
					  size_t num_neutrons)
{
	log_entry("new_combination_table(%s, %lu, %lu)",
		  filename,
		  num_protons,
		  num_neutrons);
	FILE *table_file = fopen(filename,"r");
	if (table_file == NULL)
		error("Could not open \"%s\". %s\n",
		      filename,
		      strerror(errno));
	combination_table_t table = 
		(combination_table_t)
		calloc(1,sizeof(struct _combination_table_));
	array_builder_t id_to_index_map_builder =
		new_array_builder((void**)&table->id_to_index_map,
				  &table->num_ids,
				  sizeof(size_t));
		array_builder_t basis_blocks_builder = 
		new_array_builder((void**)&table->basis_blocks,
				  &table->num_basis_blocks,
				  sizeof(basis_block_t));
	read_basis_blocks(table_file,
			  basis_blocks_builder,
			  id_to_index_map_builder,
			  num_protons,
			  num_neutrons);
	free_array_builder(basis_blocks_builder);
	array_builder_t index_list_settings_builder =
		new_array_builder((void**)
				  &table->index_list_settings,
				  &table->num_index_list_settings,
				  sizeof(index_list_setting_t));
	read_index_list_settings(table_file,
				 index_list_settings_builder,
				 id_to_index_map_builder);
	free_array_builder(index_list_settings_builder);
	array_builder_t matrix_block_settings_builder =
		new_array_builder((void**)
				  &table->matrix_block_settings,
				  &table->num_matrix_block_settings,
				  sizeof(matrix_block_setting_t));
	read_matrix_block_settings(table_file,
				   matrix_block_settings_builder,
				   id_to_index_map_builder);
	fclose(table_file);
	free_array_builder(matrix_block_settings_builder);
	sort_matrix_blocks_on_num_particles(table,id_to_index_map_builder);
	free_array_builder(id_to_index_map_builder);
	return table;
}

size_t get_full_dimension(combination_table_t combination_table)
{
	size_t full_dimension = 0;
	for (size_t i = 0; i<combination_table->num_basis_blocks; i++)
		full_dimension+=
			combination_table->basis_blocks[i].num_proton_states*
			combination_table->basis_blocks[i].num_neutron_states;	
	return full_dimension;
}

particle_type_t get_index_list_type(combination_table_t combination_table,
				    size_t index_list_id)
{
	size_t index = combination_table->id_to_index_map[index_list_id];
	return combination_table->index_list_settings[index].type;
}

basis_block_t get_basis_block(combination_table_t combination_table,
			      size_t basis_block_id)
{
	size_t index = combination_table->id_to_index_map[basis_block_id];
	return combination_table->basis_blocks[index];
}

void reset_index_list_interator(combination_table_t combination_table)
{
	combination_table->current_index_list_setting = 0;
}

int has_next_index_list_setting(combination_table_t combination_table)
{
	return combination_table->current_index_list_setting <
	       	combination_table->num_index_list_settings;
}

index_list_setting_t 
next_index_list_setting(combination_table_t combination_table)
{
	size_t index = combination_table->current_index_list_setting++;
	return combination_table->index_list_settings[index];
}

size_t get_num_basis_blocks(combination_table_t combination_table)
{
	return combination_table->num_basis_blocks;
}
void reset_basis_block_iterator(combination_table_t combination_table)
{
	combination_table->current_basis_block_index = 0;	
}

basis_block_t next_basis_block(combination_table_t table)
{
	assert(table->current_basis_block_index <
	       table->num_basis_blocks);
	return table->basis_blocks[table->current_basis_block_index++];
}

int has_next_basis_block(combination_table_t table)
{
	return table->current_basis_block_index < table->num_basis_blocks;
}

size_t get_num_matrix_blocks(combination_table_t combination_table)
{
	return combination_table->num_matrix_block_settings;
}

matrix_block_setting_t get_matrix_block_by_id(combination_table_t table,
					      size_t id)
{
	size_t index = table->id_to_index_map[id];
	assert(index < table->num_matrix_block_settings);
	return table->matrix_block_settings[index];
}

void reset_matrix_block_iterator(combination_table_t combination_table)
{
	combination_table->current_matrix_block_index = 0;
}

matrix_block_setting_t next_matrix_block(combination_table_t table)
{
	assert(table->current_matrix_block_index <
	       table->num_matrix_block_settings);
	return table->matrix_block_settings[table->
		current_matrix_block_index++];
}

int has_next_matrix_block(combination_table_t table)
{
	return table->current_matrix_block_index <
		table->num_matrix_block_settings;
}

void reset_1nf_block_iterator(combination_table_t combination_table)
{
	combination_table->iterator_index_1nf_blocks = 0;
}

matrix_block_setting_t
next_1nf_block_iterator(combination_table_t combination_table)
{
	size_t index = combination_table->iterator_index_1nf_blocks++;
	return combination_table->matrix_block_settings[index];
}

int has_next_1nf_block(combination_table_t combination_table)
{
	return combination_table->iterator_index_1nf_blocks <
	       	combination_table->length_1nf_blocks;
}

void reset_2nf_block_iterator(combination_table_t combination_table)
{
	combination_table->iterator_index_2nf_blocks = 0;
}

matrix_block_setting_t
next_2nf_block_iterator(combination_table_t combination_table)
{
	size_t index = combination_table->index_to_2nf_blocks + 
		combination_table->iterator_index_2nf_blocks++;
	return combination_table->matrix_block_settings[index];
}

int has_next_2nf_block(combination_table_t combination_table)
{
	return combination_table->iterator_index_2nf_blocks <
	       	combination_table->length_2nf_blocks;
}

void reset_3nf_block_iterator(combination_table_t combination_table)
{
	combination_table->iterator_index_3nf_blocks = 0;
}

matrix_block_setting_t
next_3nf_block_iterator(combination_table_t combination_table)
{
	size_t index = combination_table->index_to_3nf_blocks + 
		combination_table->iterator_index_3nf_blocks++;
	return combination_table->matrix_block_settings[index];
}

int has_next_3nf_block(combination_table_t combination_table)
{
	return combination_table->iterator_index_3nf_blocks <
	       	combination_table->length_3nf_blocks;
}

size_t get_num_arrays(combination_table_t combination_table)
{
	return combination_table->num_ids;
}

void free_combination_table(combination_table_t combination_table)
{
	free(combination_table->basis_blocks);
	free(combination_table->index_list_settings);
	free(combination_table->matrix_block_settings);
	free(combination_table->id_to_index_map);
	free(combination_table);
}

static
void read_basis_blocks(FILE *table_file,
		       array_builder_t basis_blocks_builder,
		       array_builder_t id_to_index_map_builder,
		       size_t num_protons,
		       size_t num_neutrons)
{
	move_to_title(table_file,"mp-states");
	char *current_row = NULL;
	size_t current_row_length = 0;
	while (!feof(table_file))
	{
		if (getline(&current_row,
			    &current_row_length,
			    table_file) < 0 ||
		    is_title_row(current_row))
			break;
		// The only rows we are interested in are the 
		// once containing ARRAYMP
		if (strstr(current_row,"ARRAYMP:") == NULL)
			continue;	
		char **words = NULL;
		size_t num_words = extract_words(&words,current_row);
		basis_block_t current_block =
		{
			.Ep = atoi(words[0]),
			.Mp = atoi(words[1]),
			.En = atoi(words[2]),
			.Mn = atoi(words[3]),
			.num_proton_states = atoll(words[4]),
			.num_neutron_states = atoll(words[5]),
			.num_protons = num_protons,
			.num_neutrons = num_neutrons,
			.block_id = interpret_id_string(words[9])
		};			
		size_t index = num_array_elements(basis_blocks_builder);
		set_array_element(id_to_index_map_builder,
				  current_block.block_id,
				  &index);
		append_array_element(basis_blocks_builder,
				     &current_block);
		if (words != NULL)
		{
			for (size_t i = 0; i<num_words; i++)
				free(words[i]);
			free(words);
		}

	}
	if (current_row != NULL)
		free(current_row);
}

static
void read_index_list_settings(FILE *table_file,
			      array_builder_t index_list_settings_builder,
			      array_builder_t id_to_index_map_builder)
{
	move_to_title(table_file,"Conn lists");
	char *current_row = NULL;
	size_t current_row_length = 0;
	while (!feof(table_file))
	{
		if (getline(&current_row,
			    &current_row_length,
			    table_file) < 0 ||
		    is_title_row(current_row))
			break;
		log_entry("current_row = %s",current_row);
		if (strstr(current_row,"ARRAY:") == NULL)
			continue;
		char **words = NULL;
		size_t num_words = extract_words(&words,current_row);
		log_entry("num_words = %lu\n",num_words);
		for (size_t i = 0; i<num_words; i++)
			log_entry("Word (%lu): %s\n",i,words[i]);
		log_entry("parse_particle_type(%s) = %d (NEUTRON:%d,PROTON:%d)",
			  words[0],
			  parse_particle_type(words[0]),
			  NEUTRON,
			  PROTON);
		index_list_setting_t current_setting =
		{
			.type = parse_particle_type(words[0]),
			.block_type = parse_block_type(words[0]),
			.energy_bra = atoi(words[1]),
			.M_bra = atoi(words[2]),
			.energy_ket = atoi(words[3]),
			.M_ket = atoi(words[4]),
			.depth = atoi(words[5]),
			.length = atoll(words[6]),
			.index_list_id = interpret_id_string(words[9])
		};
		size_t index = num_array_elements(index_list_settings_builder);
		set_array_element(id_to_index_map_builder,
				  current_setting.index_list_id,
				  &index);
		append_array_element(index_list_settings_builder,
				     &current_setting);
		if (words != NULL)
		{

			for (size_t i = 0; i<num_words; i++)
				free(words[i]);
			free(words);
		}
	}
	if (current_row != NULL)
		free(current_row);
}

static
void read_matrix_block_settings(FILE *table_file,
				array_builder_t matrix_block_settings_builder,
				array_builder_t id_to_index_map_builder)
{
	move_to_title(table_file,"Matrix-elements V (cross p-n)");
	char *current_row = NULL;
	size_t current_row_length = 0;
	while (!feof(table_file))
	{
		if (getline(&current_row,
			    &current_row_length,
			    table_file) <0 ||
		    is_title_row(current_row))
			break;
		log_entry("current_row = %s",current_row);
		if (strstr(current_row,"ARRAY:") == NULL)
			continue;
		char **words = NULL;
		size_t num_words = extract_words(&words,current_row);
		char *type_string = concatinate_strings(words[4],words[0]);
		matrix_block_setting_t current_setting =
		{
			.type = parse_block_type(type_string),
			.difference_energy_protons = atoi(words[1]),
			.difference_M_protons = atoi(words[2]),
			.depth_protons = atoi(words[3]),
			.difference_energy_neutrons = atoi(words[5]),
			.difference_M_neutrons = atoi(words[6]),
			.depth_neutrons = atoi(words[7]),
			.num_proton_combinations = atoll(words[8]),
			.num_neutron_combinations = atoll(words[9]),
			.matrix_block_id = interpret_id_string(words[13])
		};
		//size_t index =
		//       	num_array_elements(matrix_block_settings_builder);
		//set_array_element(id_to_index_map_builder,
		//		  current_setting.matrix_block_id,
		//		  &index);
		append_array_element(matrix_block_settings_builder,
				     &current_setting);
		if (words != NULL)
		{
			for (size_t i = 0; i<num_words; i++)
				free(words[i]);
			free(words);
		}
		free(type_string);
		
	}
	move_to_title(table_file,"Matrix-elements V (same p/n)");
	while (!feof(table_file))
	{
		if (getline(&current_row,
			    &current_row_length,
			    table_file) <0)
		       break;
		if (is_title_row(current_row))
			break;
		if (strstr(current_row,"ARRAY:") == NULL)
			continue;
		char **words = NULL;
		size_t num_words = extract_words(&words,current_row);
		int differnce_energy = atoi(words[1]);
		int differnce_M = atoi(words[2]);
		int depth = atoi(words[3]);
		size_t num_combinations = atoi(words[4]);
		matrix_block_setting_t current_setting;
		if (*(words[0]) == 'p')
		{
			matrix_block_setting_t proton_setting =
			{
				.type = parse_block_type(words[0]),
				.difference_energy_protons = differnce_energy,
				.difference_M_protons = differnce_M,
				.depth_protons = depth,
				.difference_energy_neutrons = 0,
				.difference_M_neutrons = 0,
				.depth_neutrons = 0,
				.num_proton_combinations = num_combinations,
				.num_neutron_combinations = 0,
				.matrix_block_id =
				       	interpret_id_string(words[8])
			};
			current_setting = proton_setting;
		}
		else
		{
			matrix_block_setting_t proton_setting =
			{
				.type = parse_block_type(words[0]),
				.difference_energy_protons = 0,
				.difference_M_protons = 0,
				.depth_protons = 0,
				.difference_energy_neutrons = differnce_energy,
				.difference_M_neutrons = differnce_M,
				.depth_neutrons = depth,
				.num_proton_combinations = 0,
				.num_neutron_combinations = num_combinations,
				.matrix_block_id =
				       	interpret_id_string(words[8])
			};
			current_setting = proton_setting;
		}
		//size_t index =
		//       	num_array_elements(matrix_block_settings_builder);
		//set_array_element(id_to_index_map_builder,
		//		  current_setting.matrix_block_id,
		//		  &index);
		append_array_element(matrix_block_settings_builder,
				     &current_setting);
		if (words != NULL)
		{
			for (size_t i = 0; i<num_words; i++)
				free(words[i]);
			free(words);
		}
		
	}
	if (current_row != NULL)
		free(current_row);
}

static void 
sort_matrix_blocks_on_num_particles(combination_table_t table,
				    array_builder_t id_to_index_map_builder)
{
	size_t bucket_sizes[3] = {0};
	for (size_t i = 0; i<table->num_matrix_block_settings; i++)
	{
		const block_type_t type = table->matrix_block_settings[i].type;
		bucket_sizes[count_particles(type)-1]++; 
	}
	size_t bucket_positions[3] = {0};
	for (size_t i = 1; i<3; i++)
		bucket_positions[i]=bucket_sizes[i-1]+bucket_positions[i-1];
	table->length_1nf_blocks = bucket_sizes[0];
	table->length_2nf_blocks = bucket_sizes[1];
	table->length_3nf_blocks = bucket_sizes[2];
	table->index_to_2nf_blocks = bucket_positions[1];
	table->index_to_3nf_blocks = bucket_positions[2];
	table->iterator_index_1nf_blocks = 0;
	table->iterator_index_2nf_blocks = 0;
	table->iterator_index_3nf_blocks = 0;
	qsort(table->matrix_block_settings,
	      table->num_matrix_block_settings,
	      sizeof(matrix_block_setting_t),
	      (__compar_fn_t)compare_matrix_block_settings);
	for (size_t index = 0; index<table->num_matrix_block_settings; index++)
	{
		matrix_block_setting_t current_block =
			table->matrix_block_settings[index];
		set_array_element(id_to_index_map_builder,
				  current_block.matrix_block_id,
				  &index);
	}

}

static
int is_title_row(const char *row)
{
	log_entry("row = %s",row);
	log_entry("strstr(\"***\",row) = %s",
		  strstr(row,"***"));  
	return strstr(row,"***") != NULL;
}

static
void move_to_title(FILE *table_file,
		   const char *title)
{
	log_entry("Moving to title %s",
		  title);
	fseek(table_file,0,SEEK_SET);
	char *current_row = NULL;
	size_t current_row_length = 0;
	while (!feof(table_file))
	{
		if (getline(&current_row,
			    &current_row_length,
			    table_file) < 0)
			break;
		if (strstr(current_row,title) != NULL)
			break;	
	}
	if (current_row != NULL)
		free(current_row);
}

static
size_t interpret_id_string(const char *id_string)
{
	log_entry("id_string = %s",id_string);
	size_t id = 0;
	while (*id_string != '=' && *id_string != 0)
		id = 10*id+(size_t)(*(id_string++)-'0');
	return id;
}

static
int compare_matrix_block_settings(matrix_block_setting_t *block_a,
				  matrix_block_setting_t *block_b)
{
	int diff = count_particles(block_a->type) -
	       	count_particles(block_b->type);
	if (diff)
		return diff;
	diff = count_neutrons(block_a->type) - count_neutrons(block_b->type);
	if (diff)
		return diff;
	diff = block_a->depth_protons - block_b->depth_protons;
	if (diff)
		return diff;
	diff = block_a->difference_energy_protons -
	       	block_b->difference_energy_protons;
	if (diff)
		return diff;
	diff = block_a->depth_neutrons - block_b->depth_neutrons;
	if (diff)
		return diff;
	diff = block_a->difference_energy_neutrons -
	       	block_b->difference_energy_neutrons;
	if (diff)
		return diff;
	diff = block_a->difference_M_protons -block_b->difference_M_protons;
	if (diff)
		return diff;
	diff = block_a->difference_M_neutrons -block_b->difference_M_neutrons;
	if (diff)
		return diff;
	return 0;
}

new_test(interpreting_id_string,
	 const char* input = "12=";
	 assert_that(interpret_id_string(input) == 12);
	);

new_test(is_title_row_should_work,
	 const char *row = "*** Matrix-elements V (same p/n) ***\n";
	 assert_that(is_title_row(row));
	);


