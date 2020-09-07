#include <radix_sort/radix_sort.h>
#include <string.h>
#include <debug_mode/debug_mode.h>


typedef struct
{
	uint64_t key;
	size_t index;
} key_index_t;

static inline
void setup_buckets(size_t buckets[256],
		   const key_index_t *current_order,
		   const size_t array_length);

static inline
void setup_next_order(key_index_t *next_order,
		      const key_index_t *current_order,
		      size_t buckets[256],
		      const size_t array_length);

static inline
void swap(void **a,void **b);

void rsort_r(void *array,
	     size_t array_length,
	     size_t element_size,
	     __key_function_r_t keyfunction,
	     void *data)
{
	key_index_t *current_order =
	       	(key_index_t*)malloc(array_length*sizeof(key_index_t));
	key_index_t *next_order =
		(key_index_t*)malloc(array_length*sizeof(key_index_t));
	char *byte_array = (char*)array;
	uint64_t max_key = 0;
	for (size_t i = 0; i < array_length; i++)
	{
		current_order[i].index = i;
		current_order[i].key =
		       	keyfunction((void*)(byte_array+element_size*i),
				    data);
		if (current_order[i].key > max_key)
			max_key = current_order[i].key;
	}
	while (max_key > 0)
	{
		size_t buckets[16];
		setup_buckets(buckets,current_order,array_length);
		setup_next_order(next_order,current_order,buckets,array_length);
		swap((void**)&next_order,(void**)&current_order);
		max_key >>=4;
	}
	char *sorted_array = (char*)malloc(array_length*element_size);	
	for (size_t i = 0; i < array_length; i++)
		memcpy(sorted_array+i*element_size,
		       byte_array+current_order[i].index*element_size,
		       element_size);
	memcpy(array,
	       sorted_array,
	       array_length*element_size);		
	free(sorted_array);
	free(current_order);
	free(next_order);
}

static inline
void setup_buckets(size_t buckets[16],
		   const key_index_t *current_order,
		   const size_t array_length)
{
	size_t bucket_sizes[16] = {0};
	for (size_t i = 0; i<array_length; i++)
	{
		bucket_sizes[current_order[i].key & 15]++;
	}
	buckets[0] = 0;
	for (size_t i = 1; i<16; i++)
		buckets[i]=bucket_sizes[i-1] + buckets[i-1];
}

static inline
void setup_next_order(key_index_t *next_order,
		      const key_index_t *current_order,
		      size_t buckets[16],
		      const size_t array_length)
{
	for (size_t i = 0; i<array_length; i++)
	{
		key_index_t current = current_order[i];
		size_t j = buckets[current.key & 15]++;
		current.key >>= 4;
		next_order[j] = current;		
	}
}

static inline
void swap(void **a,void **b)
{
	void *tmp = *a;
	*a = *b;
	*b = tmp;
}
