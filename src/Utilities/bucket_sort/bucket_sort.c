#include <bucket_sort/bucket_sort.h>
#include <string.h>

void bucket_sort(void *array,
		 size_t num_elements,
		 size_t size_element,
		 size_t num_buckets,
		 bucket_index_function_t bucket_index_function)
{
	size_t *buckets = (size_t*)calloc(num_buckets,sizeof(size_t));
	char *array_bytes = (char*)array;
	for (size_t i = 0; i<num_elements; i++)
		buckets[bucket_index_function((void*)(array_bytes+
						      i*size_element))-1]++;
	size_t position = 0;
	for (size_t i = 0; i<num_buckets; i++)
	{
		size_t current_position = position;
		position+=buckets[i];
		buckets[i] = current_position;
	}
	char *sorted_array = malloc(num_elements*size_element);
	for (size_t i = 0; i < num_elements; i++)
	{
		void *current_element = (void*)(array_bytes+i*size_element);
		size_t target_index =
			buckets[bucket_index_function(current_element)-1]++;
		void *target_element =
		       	(void*)(sorted_array+target_index*size_element);
		memcpy(target_element,
		       current_element,
		       size_element);
	}
	memcpy(array,
	       sorted_array,
	       num_elements*size_element);
	free(buckets);
	free(sorted_array);
}
