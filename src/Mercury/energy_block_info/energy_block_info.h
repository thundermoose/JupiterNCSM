#ifndef __ENERGY_BLOCK_INFO__
#define __ENERGY_BLOCK_INFO__

typedef struct
{
	int dE, E1, E2, Tz, M;
	char configuration_file[128];
	char element_file[128];
	int is_old;
} energy_block_info_t;

int compare_energy_block_info(const energy_block_info_t *block_info_one,
			      const energy_block_info_t *block_info_two);

#endif
