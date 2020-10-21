#ifndef __BUCKET_SORT__
#define __BUCKET_SORT__

#include <stdlib.h>

typedef size_t (*bucket_index_function_t)(void *element);

void bucket_sort(void *array,
		 size_t num_elements,
		 size_t size_element,
		 size_t num_buckets,
		 bucket_index_function_t bucket_index_function);
#endif
