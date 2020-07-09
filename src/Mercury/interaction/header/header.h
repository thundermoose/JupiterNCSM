#ifndef __HEADER__
#define __HEADER__

#include <stdlib.h>
#include <energy_block_info/energy_block_info.h>

struct _header_;
typedef struct _header_ *header_t;

header_t read_header(const char *interaction_path);

size_t get_num_particles(const header_t header);

size_t find_block_index(const header_t header,
			int E1, int E2, int Tz, int M);

energy_block_info_t get_energy_block_info(header_t header,
					  size_t index);

void free_header(header_t header);

#endif
