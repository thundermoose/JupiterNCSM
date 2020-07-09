#ifndef __INDEX_LIST_SETTING__
#define __INDEX_LIST_SETTING__

#include <particle_type/particle_type.h>
#include <block_type/block_type.h>

typedef struct
{
	particle_type_t type;
	block_type_t block_type;
	int energy_bra;
	int M_bra;
	int energy_ket;
	int M_ket;
	int depth;
	size_t length;
	size_t index_list_id;
} index_list_setting_t;

#endif
