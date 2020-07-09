#ifndef __BIT_OPERATIONS__
#define __BIT_OPERATIONS__

#include <stdlib.h>
#include <stdint.h>
#include <thundertester/test.h>

static inline
uint64_t cyclic_rshift(uint64_t number,
		const size_t shift)
{
	const size_t modular_shift = shift % 64;
	return (number>>modular_shift) | (number<<(64-modular_shift));
}

static inline
uint64_t cyclic_lshift(uint64_t number,
		const size_t shift)
{
	const size_t modular_shift = shift % 64;
	return (number<<modular_shift) | (number>>(64-modular_shift));
}

new_test(cyclic_rshift_shift_by_8,
		const uint64_t number = 0xFEDCBA9876543210;
		const uint64_t expected_result = 0x10FEDCBA98765432;
		const uint64_t result = cyclic_rshift(number,8);
		printf("%lX -> %lX = %lX\n",number,result,expected_result);
		assert_that(result == 
			expected_result));

new_test(cyclic_lshift_shift_by_8,
		assert_that(cyclic_lshift(0xFEDCBA9876543210,8) == 
			0xDCBA9876543210FE));
#endif
