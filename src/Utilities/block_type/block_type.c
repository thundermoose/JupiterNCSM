#include <block_type/block_type.h>
#include <error/error.h>
#include <string.h>
#include <unit_testing/test.h>
#include <debug_mode/debug_mode.h>

#define NUM_TYPES 9

static
const char *type_names[NUM_TYPES] =
{
	"np",
	"n",
	"p",
	"nnp",
	"npp",
	"nn",
	"pp",
	"nnn",
	"ppp"
};

block_type_t parse_block_type(const char *string)
{
	for (int i = 0; i<NUM_TYPES; i++)
	{
		if (strcmp(string,type_names[i]) == 0)
			return (block_type_t)(i);
	}
	error("Unknown block type %s\n",string);
}

const char *block_type_to_string(block_type_t type)
{
	return type_names[type];
}

size_t count_neutrons(block_type_t type)
{
	switch (type)
	{
		case N_block:
		case NP_block:
		case NPP_block:
			return 1;
		case NN_block:
		case NNP_block:
			return 2;
		case NNN_block:
			return 3;	
		default:
			return 0;
	}
}

size_t count_protons(block_type_t type)
{
	switch (type)
	{
		case P_block:
		case NP_block:
		case NNP_block:
			return 1;
		case PP_block:
		case NPP_block:
			return 2;
		case PPP_block:
			return 3;	
		default:
			return 0;
	}
}

size_t count_particles(block_type_t type)
{
	switch (type)
	{
		case N_block:
		case P_block:
			return 1;
		case NN_block:
		case NP_block:
		case PP_block:
			return 2;
		case NNN_block:
		case NNP_block:
		case NPP_block:
		case PPP_block:
			return 3;
		default:
			return 0;
	}
}

new_test(test_parse_block_type_nnp,
	 assert_that(parse_block_type("nnp") == NNP_block);
	);
