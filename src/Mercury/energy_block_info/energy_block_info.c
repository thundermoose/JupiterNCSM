#include <energy_block_info/energy_block_info.h>
#include <debug_mode/debug_mode.h>

int compare_energy_block_info(const energy_block_info_t *block_info_one,
			      const energy_block_info_t *block_info_two)
{
	int diff = block_info_one->dE - block_info_two->dE;
	if (diff)
		return diff;
	diff = block_info_one->E1 - block_info_two->E1;
	if (diff)
		return diff;	
	diff = block_info_one->E2 - block_info_two->E2;
	if (diff)
		return diff;	
	diff = block_info_one->Tz - block_info_two->Tz;
	if (diff)
		return diff;	
	diff = block_info_one->M - block_info_two->M;
	if (diff)
		return diff;	
	return 0;
}
