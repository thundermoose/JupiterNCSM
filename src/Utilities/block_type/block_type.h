#ifndef __BLOCK_TYPE__
#define __BLOCK_TYPE__

#include <stdlib.h>

typedef enum
{
	NP_block,
	N_block,
	P_block,
	NNP_block,
	NPP_block,
	NN_block,
	PP_block,
	NNN_block,
	PPP_block	
} block_type_t;

block_type_t parse_block_type(const char *string);

const char *block_type_to_string(block_type_t type);

size_t count_neutrons(block_type_t type);

size_t count_protons(block_type_t type);

#endif
