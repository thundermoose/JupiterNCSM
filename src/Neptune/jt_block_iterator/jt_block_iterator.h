#ifndef __JT_BLOCK_ITERATOR__
#define __JT_BLOCK_ITERATOR__

#include <bases/shells.h>
#include <stdint.h>

typedef struct
{
	uint64_t T:1;
	uint64_t J:32;
} jt_block_t;

typedef union
{
	size_t index;
	jt_block_t block;
} jt_block_iterator_t;

jt_block_iterator_t initial_block(quantum_number j_a,
				  quantum_number j_b,
				  quantum_number m_a,
				  quantum_number m_b,
				  quantum_number tz_a,
				  quantum_number tz_b);

int has_finished(jt_block_iterator_t block_iterator,
	       quantum_number j_a,
	       quantum_number j_b);

jt_block_t get_next_block(jt_block_iterator_t *iterator,
			  quantum_number tz_a,
			  quantum_number tz_b);
#endif
