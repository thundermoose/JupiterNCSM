#ifndef __EXECUTION_ORDER__
#define __EXECUTION_ORDER__

#include <stdlib.h>
#include <combination_table/combination_table.h>

struct _evaluation_order_;
typedef struct _evaluation_order_ *evaluation_order_t;

struct _evaluation_order_iterator_;
typedef struct _evaluation_order_iterator_ *evaluation_order_iterator_t;

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
	size_t instruction_index;
} evaluation_instruction_t;

evaluation_order_t read_evaluation_order(const char *filename,
				       combination_table_t combination_table);

size_t get_num_instructions(evaluation_order_t evaluation_order);

evaluation_order_iterator_t 
get_evaluation_order_iterator(evaluation_order_t evaluation_order);

void reset_evaluation_order(evaluation_order_iterator_t iterator);

evaluation_instruction_t 
next_instruction(evaluation_order_iterator_t iterator);

int has_next_instruction(evaluation_order_iterator_t iterator);

void free_evaluation_order_iterator(evaluation_order_iterator_t iterator);

void free_evaluation_order(evaluation_order_t evaluation_order);

#endif
