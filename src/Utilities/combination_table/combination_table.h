#ifndef __COMBINATION_TABLE__
#define __COMBINATION_TABLE__

#include <stdlib.h>
#include <particle_type/particle_type.h>
#include <basis_block/basis_block.h>
#include <index_list_setting/index_list_setting.h>
#include <matrix_block_setting/matrix_block_setting.h>

struct _combination_table_;
typedef struct _combination_table_ *combination_table_t;

combination_table_t new_combination_table(const char *filename,
					  size_t num_protons,
					  size_t num_neutrons);
size_t get_full_dimension(combination_table_t combination_table);

particle_type_t get_index_list_type(combination_table_t combination_table,
				    size_t index_list_id);

basis_block_t get_basis_block(combination_table_t combination_table,
			      size_t basis_block_id);

void reset_index_list_interator(combination_table_t combination_table);

int has_next_index_list_setting(combination_table_t combination_table);

index_list_setting_t 
next_index_list_setting(combination_table_t combination_table);

size_t get_num_basis_blocks(combination_table_t combination_table);

void reset_basis_block_iterator(combination_table_t combination_table);

basis_block_t next_basis_block(combination_table_t combination_table);

int has_next_basis_block(combination_table_t combination_table);

size_t get_num_matrix_blocks(combination_table_t combination_table);

matrix_block_setting_t get_matrix_block_by_id(combination_table_t table,
					      size_t id);

void reset_matrix_block_iterator(combination_table_t combination_table);

matrix_block_setting_t 
next_matrix_block(combination_table_t combination_table);

int has_next_matrix_block(combination_table_t combination_table);

void reset_1nf_block_iterator(combination_table_t combination_table);

matrix_block_setting_t
next_1nf_block_iterator(combination_table_t combination_table);

int has_next_1nf_block(combination_table_t combination_table);

void reset_2nf_block_iterator(combination_table_t combination_table);

matrix_block_setting_t
next_2nf_block_iterator(combination_table_t combination_table);

int has_next_2nf_block(combination_table_t combination_table);

void reset_3nf_block_iterator(combination_table_t combination_table);

matrix_block_setting_t
next_3nf_block_iterator(combination_table_t combination_table);

int has_next_3nf_block(combination_table_t combination_table);

size_t get_num_arrays(combination_table_t combination_table);

void free_combination_table(combination_table_t combination_table);

#endif
