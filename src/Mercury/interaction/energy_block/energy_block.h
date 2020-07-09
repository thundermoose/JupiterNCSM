#ifndef __ENERGY_BLOCK__
#define __ENERGY_BLOCK__

#include <stdlib.h>
#include <energy_block_info/energy_block_info.h>

struct _energy_block_;
typedef struct _energy_block_ *energy_block_t;

energy_block_t open_energy_block(const char *interaction_path,
				 energy_block_info_t info);

double get_energy_block_element(energy_block_t energy_block,
				size_t bra_index,
				size_t ket_index);

void free_energy_block(energy_block_t energy_block);

#endif
