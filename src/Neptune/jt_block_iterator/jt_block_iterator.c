#include <jt_block_iterator/jt_block_iterator.h>
#include <utils/helpful_macros.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>

jt_block_iterator_t initial_block(quantum_number j_a,
				  quantum_number j_b,
				  quantum_number m_a,
				  quantum_number m_b,
				  quantum_number tz_a,
				  quantum_number tz_b)
{
	const quantum_number J_min = max(abs(j_a-j_b),
					 abs(m_a+m_b))/2;
	const quantum_number T_min = abs(tz_a+tz_b)/2;
	log_entry("J_min = %d, T_min = %d\n",J_min,T_min);
	jt_block_iterator_t iterator;
	iterator.block.T = T_min;
	iterator.block.J = J_min;
	log_entry("block: %u %u\n",
		   iterator.block.J,iterator.block.T);
	return iterator;
}

int has_finished(jt_block_iterator_t block_iterator,
	       quantum_number j_a,
	       quantum_number j_b)
{
	const quantum_number J_max = (j_a+j_b)/2;
	return block_iterator.block.J>J_max;
}

jt_block_t get_next_block(jt_block_iterator_t *iterator,
			  quantum_number tz_a,
			  quantum_number tz_b)
{
	iterator->index+= 1+abs(tz_a+tz_b)/2;
	return iterator->block;
}
