#ifndef __ARRAY_BUILDER__
#define __ARRAY_BUILDER__

#include <stdlib.h>

/* The main handle type of the array_builder.
 */
struct _array_builder_;
typedef struct _array_builder_ *array_builder_t;

/* To define custom resize policies.
 */
typedef size_t (*resize_policy_t)(size_t current_length,
				  size_t needed_expansion);

/* To set up an array builder provide a pointer to your array pointer
 * and one to the length of the array.
 */
array_builder_t new_array_builder(void **array,
				  size_t *array_length,
				  size_t element_size);

/* This is the default resize policy, the output is equal to
 * 2*current_length + needed_length which means that the array increases 
 * exponentially, which is good for fast construction of small arrays, but
 * could cause problems if the final array is in the same size as max ram 
 * memory.
 */
size_t resize_exponentialy(size_t current_length,
			   size_t needed_length);

/* To set a new resize policy.
 */
void set_resize_policy(array_builder_t builder,
		       resize_policy_t policy);

/* To get the current size of the array, the size of the array that it will
 * be trimmed down to if free_array_builder is called.
 */
size_t num_array_elements(const array_builder_t builder);

/* To add a new element to the end of the current array.
 */
void append_array_element(array_builder_t builder,
			  void *element);

/* To set an element with a specific index. If the index is larger than
 * the current size, the array is rescaled.
 */
void set_array_element(array_builder_t builder,
		       size_t index,
		       void *element);
/* Freeing the array builder, this will also trim the array, and
 * free up unused memory.
 */
void free_array_builder(array_builder_t builder);
#endif
