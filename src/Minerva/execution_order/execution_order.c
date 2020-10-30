#include <execution_order/execution_order.h>
#include <array_builder/array_builder.h>
#include <string_tools/string_tools.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

const execution_instruction_t empty_instruction = {
	.type = unknown,
	.vector_block_in = no_index,
	.vector_block_out = no_index,
	.matrix_element_file = no_index,
	.neutron_index = no_index,
	.proton_index = no_index
};

struct _execution_order_
{
	execution_instruction_t *instructions;
	size_t num_instruction;
	size_t current_instruction_index;
};

static
void parse_row(const char *row,
	       array_builder_t instructions_builder,
	       combination_table_t combination_table);

execution_order_t read_execution_order(const char *filename,
				       combination_table_t combination_table)
{
	execution_order_t execution_order =
		(execution_order_t)calloc(1,sizeof(struct _execution_order_));
	array_builder_t instructions_builder =
		new_array_builder((void**)&execution_order->instructions,
				  &execution_order->num_instruction,
				  sizeof(execution_instruction_t));
	FILE *file = fopen(filename,"r");
	char *row = NULL;
	size_t row_length = 0;
	while (!feof(file))
	{
		if (getline(&row,&row_length,file) < 0)
			break;	
		parse_row(row,
			  instructions_builder,
			  combination_table);
	}	
	fclose(file);
	if (row)
		free(row);
	free_array_builder(instructions_builder);
	for (size_t i = 0; i<execution_order->num_instruction-1; i++)
	{
		execution_instruction_t instruction =
		       	execution_order->instructions[i];
		execution_instruction_t next =
			execution_order->instructions[i+1];
		execution_instruction_t next_next = 
			i < execution_order->num_instruction-3 ?
			execution_order->instructions[i+3] :
			empty_instruction;
		if (instruction.type != unload)
			continue;
		if (next.type == unload)
			continue;
		if (instruction.vector_block_in == next.vector_block_in ||
		    instruction.vector_block_in == next_next.vector_block_in)
		{
			instruction.vector_block_in = no_index;
		}
		if (instruction.vector_block_out == next.vector_block_out ||
		    instruction.vector_block_out == next_next.vector_block_out)
		{
			instruction.vector_block_out = no_index;
		}
		if (instruction.matrix_element_file == next.matrix_element_file ||
		    instruction.matrix_element_file == next_next.matrix_element_file)
		{
			instruction.matrix_element_file = no_index;
		}
		if (instruction.neutron_index == next.neutron_index ||
		    instruction.neutron_index == next_next.neutron_index)
		{
			instruction.neutron_index = no_index;
		}
		if (instruction.proton_index == next.proton_index ||
		    instruction.proton_index == next_next.proton_index)
		{
			instruction.proton_index = no_index;
		}
		execution_order->instructions[i] = instruction;
	}
	return execution_order;
}

void reset_execution_order(execution_order_t execution_order)
{
	execution_order->current_instruction_index = 0;
}

execution_instruction_t next_instruction(execution_order_t execution_order)
{
	assert(execution_order->current_instruction_index <
	       execution_order->num_instruction);
	size_t index = execution_order->current_instruction_index++;
	return execution_order->instructions[index];
}

int has_next_instruction(execution_order_t execution_order)
{
	return execution_order->current_instruction_index < 
		execution_order->num_instruction;	
}

void free_execution_order(execution_order_t execution_order)
{
	free(execution_order->instructions);
	free(execution_order);
}

	static
void parse_row(const char *row,
	       array_builder_t instructions_builder,
	       combination_table_t combination_table)
{
	log_entry("row = %s",row);
	char *block_pointer = strstr(row,"BLOCK:");
	if (block_pointer == NULL)
		return;
	char *unload_pointer = strstr(row,"UNLOAD_");	
	execution_instruction_t current_instruction;
	current_instruction.type = unknown;
	if (unload_pointer != NULL)
		current_instruction.type = unload;
	row = block_pointer;
	char **words = NULL;
	size_t num_words = extract_words(&words,row);
	current_instruction.vector_block_in = atoll(words[1]);
	current_instruction.vector_block_out = atoll(words[2]);
	if (num_words == 6)
	{
		if (current_instruction.type != unload)
			current_instruction.type = neutron_proton_block;
		size_t index_list_one = atoll(words[3]);
		size_t index_list_two = atoll(words[4]);
		current_instruction.matrix_element_file = atoll(words[5]);
		if (get_index_list_type(combination_table,
					index_list_one) == NEUTRON)
		{
			current_instruction.neutron_index = index_list_one;
			current_instruction.proton_index = index_list_two;
		}
		else
		{
			current_instruction.neutron_index = index_list_two;
			current_instruction.proton_index = index_list_one;
		}
	}
	else
	{
		size_t index_list = atoll(words[3]);
		current_instruction.matrix_element_file = atoll(words[4]);
		if (get_index_list_type(combination_table,
					index_list) == NEUTRON)
		{
			if (current_instruction.type != unload)
				current_instruction.type = neutron_block;
			current_instruction.neutron_index = index_list;
			current_instruction.proton_index = no_index;
		}
		else
		{
			if (current_instruction.type != unload)
				current_instruction.type = proton_block;
			current_instruction.proton_index = index_list;
			current_instruction.neutron_index = no_index;
		}
	}
	append_array_element(instructions_builder,
			     &current_instruction);
	if (words != NULL)
	{
		for (size_t i = 0; i<num_words; i++)
			free(words[i]);
		free(words);
	}
}
