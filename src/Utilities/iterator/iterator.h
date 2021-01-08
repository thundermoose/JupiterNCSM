#ifndef __ITERATOR__
#define __ITERATOR__

#include <stdlib.h>

struct _iterator_;
typedef struct _iterator_ *iterator_t;

iterator_t new_iterator(void *first_element,
			size_t num_elements,
			size_t size_element);

void initialize(iterator_t iterator,
		void *first_element);

int has_next_element(iterator_t iterator);

void next_element(iterator_t iterator,
		  void *current_element);

void free_iterator(iterator_t iterator);

#endif
