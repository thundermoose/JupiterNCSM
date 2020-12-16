#ifndef __RADIX_SORT__
#define __RADIX_SORT__

#include <stdlib.h>
#include <stdint.h>

typedef uint64_t (*__key_function_t)(const void *element);

typedef uint64_t (*__key_function_r_t)(const void *element,
				       void *data);

void rsort(void *array,
	   size_t array_length,
	   size_t element_size,
	   __key_function_t keyfunction);

void rsort_r(void *array,
	     size_t array_length,
	     size_t element_size,
	     __key_function_r_t keyfunction,
	     void *data);

#endif
