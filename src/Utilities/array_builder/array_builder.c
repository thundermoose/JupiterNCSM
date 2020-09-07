#include <array_builder/array_builder.h>
#include <string.h>
#include <thundertester/test.h>
#include <debug_mode/debug_mode.h>

struct _array_builder_
{
	void **array;
	size_t *array_length;
	size_t allocated_length;
	size_t element_size;
	resize_policy_t resize_policy;	
};

static
void expand_array(array_builder_t builder,
		  size_t needed_expansion);

static
void trim_array(array_builder_t builder);

array_builder_t new_array_builder(void **array,
				  size_t *array_length,
				  size_t element_size)
{
	array_builder_t builder =
	       	(array_builder_t)malloc(sizeof(struct _array_builder_));
	builder->array = array;
	builder->array_length = array_length;
	builder->element_size = element_size;
	builder->allocated_length = *array_length;
	builder->resize_policy = resize_exponentialy;
	return builder;
}

size_t resize_exponentialy(size_t current_length,
			   size_t needed_expansion)
{
	return 2*current_length + needed_expansion;
}

void set_resize_policy(array_builder_t builder,
		       resize_policy_t policy)
{
	builder->resize_policy = policy;
}

size_t num_array_elements(array_builder_t builder)
{
	return *(builder->array_length);
}

void append_array_element(array_builder_t builder,
			  void *element)
{
	expand_array(builder,1);
	char *array = *builder->array;
	array+=builder->element_size*(*builder->array_length);
	memcpy(array,element,builder->element_size);
	(*builder->array_length)++;
}

void set_array_element(array_builder_t builder,
		       size_t index,
		       void *element)
{
	if (index >= (*builder->array_length))
	{
		expand_array(builder,index-(*builder->array_length)+1);
		*builder->array_length = index+1;
	}
	char *array = *builder->array;
	array+=builder->element_size*index;	
	memcpy(array,element,builder->element_size);
}

void free_array_builder(array_builder_t builder)
{
	trim_array(builder);
	free(builder);
}

static
void expand_array(array_builder_t builder,
		  size_t needed_expansion)
{
	if ((*builder->array_length)+needed_expansion >
	    builder->allocated_length)
	{
		builder->allocated_length =
		       	builder->resize_policy(builder->allocated_length,
					       needed_expansion);
		*builder->array =
		       	realloc(*(builder->array),
				builder->element_size*
				builder->allocated_length);
	}
}

static
void trim_array(array_builder_t builder)
{
	if ((*builder->array_length) < builder->allocated_length)
	{
		builder->allocated_length = *builder->array_length;
		*builder->array =
		       	realloc(*builder->array,
				builder->element_size*
				builder->allocated_length);
	}
}

new_test(appending_twenty_elements_to_array,
	 size_t length_array = 0;
	 size_t *array = NULL;
	 array_builder_t builder = new_array_builder((void**)&array,
						     &length_array,
						     sizeof(size_t));
	 for (size_t i = 0; i<20; i++)
	 	append_array_element(builder,&i);
	 free_array_builder(builder);
	 assert_that(length_array == 20);
	 for (size_t i = 0; i<length_array; i++)
	 	assert_that(array[i] == i);
	 free(array);
	);
