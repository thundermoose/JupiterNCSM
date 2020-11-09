#include <iterator/iterator.h>
#include <omp.h>

struct _iterator_
{
	void *first_element;
	size_t num_elements;
	size_t size_element;
	size_t current_element;
	omp_nest_lock_t next_element_lock;
};

iterator_t new_iterator(void *first_element,
			size_t num_elements,
			size_t size_element)
{
	iterator_t iterator = (iterator_t)malloc(sizeof(struct _iterator_));
	iterator->first_element = first_element;
	iterator->num_elements = num_elements;
	iterator->size_element = size_element;
	iterator->current_element = 0;
	omp_init_nest_lock(&iterator->next_element_lock);
	return iterator;
}

void initialize(iterator_t iterator,
		void *first_element)
{
	omp_set_nest_lock(&iterator->next_element_lock);
	iterator->current_element = 0;
	memcpy(first_element,
	       iterator->first_element,
	       iterator->size_element);
	omp_unset_nest_lock(&iterator->next_element_lock);
}

int has_next_element(iterator_t iterator)
{
	omp_set_nest_lock(&iterator->next_element_lock);
	return iterator->current_element<iterator->num_elements;
}

void next_element(iterator_t iterator,
		  void *current_element)
{
	char *element = (char*)(iterator->first_element)+
		iterator->current_element*iterator->size_element;
	memcpy(current_element,
	       element,
	       iterator->size_element);
	omp_unset_nest_lock(&iterator->next_element_lock);
}

void free_iterator(iterator_t iterator)
{
	omp_destroy_nest_lock(&iterator->next_element_lock);
	free(iterator);
}
