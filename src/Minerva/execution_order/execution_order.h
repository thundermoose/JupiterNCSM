#ifndef __EXECUTION_ORDER__
#define __EXECUTION_ORDER__

#include <stdlib.h>
#include <combination_table/combination_table.h>

struct _execution_order_;
typedef struct _execution_order_ *execution_order_t;

struct _execution_order_iterator_;
typedef struct _execution_order_iterator_ *execution_order_iterator_t;

typedef enum
{
	unknown,
	neutron_block,
	proton_block,
	neutron_proton_block,
	unload
} instruction_type_t;

typedef struct
{
	instruction_type_t type;
	size_t vector_block_in;
	size_t vector_block_out;
	size_t matrix_element_file;
	size_t neutron_index;
	size_t proton_index;
} execution_instruction_t;

execution_order_t read_execution_order(const char *filename,
				       combination_table_t combination_table);

size_t get_num_instructions(execution_order_t execution_order);

execution_order_iterator_t 
get_execution_order_iterator(execution_order_t execution_order);

void reset_execution_order(execution_order_iterator_t iterator);

execution_instruction_t 
next_instruction(execution_order_iterator_t iterator);

int has_next_instruction(execution_order_iterator_t iterator);

void free_execution_order_iterator(execution_order_iterator_t iterator);

void free_execution_order(execution_order_t execution_order);

#endif
