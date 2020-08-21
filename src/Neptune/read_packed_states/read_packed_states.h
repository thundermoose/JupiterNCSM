#ifndef __READ_PACKED_STATES__
#define __READ_PACKED_STATES__

#include <stdlib.h>

void read_packed_states(void **states,
			size_t *num_states,
			size_t size_states,
			const char *states_filename);

#endif
