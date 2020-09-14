#include <utils/range.h>
#include <unit_testing/test.h>

size_t *range(size_t start,size_t stop)
{
	const size_t length = stop-start;
	size_t *array = (size_t*)malloc(length*sizeof(size_t));
	for (size_t i = 0; i<length; i++)
		array[i] = start+i;
	return array;
}

new_test(short_range,
	 size_t *array = range(4,10);
	 for (size_t i = 0; i<6; i++)
	 	assert_that(i+4 == array[i]);
	);
