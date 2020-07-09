#include <matrix_builder/matrix_builder.h>
#include <utils/index_hash.h>

struct _sparse_matrix_builder_
{
	Sparse_Matrix *matrix;
	size_t allocated_num_elements;
	index_hash_t index_hash_map;
};

static
int is_full(sparse_matrix_builder_t matrix_builder);

static 
void expand_sparese_matrix(sparse_matrix_builder_t matrix_builder);

static
void trim_sparse_matrix(sparse_matrix_builder_t matrix_builder);

static
void reallocate_sparse_matrix(sparse_matrix_builder_t matrix_builder);

sparse_matrix_builder_t new_sparse_matrix_builder(Sparse_Matrix *matrix)
{
	sparse_matrix_builder_t matrix_builder =
		(sparse_matrix_builder_t)
		malloc(sizeof(struct _sparse_matrix_builder_));		
	matrix_builder->matrix = matrix;
	matrix_builder->allocated_num_elements = matrix->num_elements;
	const size_t num_bins = 2*matrix->m*matrix->n;
	matrix_builder->index_hash_map = new_index_hash(num_bins);
	return matrix_builder;
}

void set_sparse_matrix_element(sparse_matrix_builder_t matrix_builder,
			     size_t m_index,
			     size_t n_index,
			     double element)
{
	size_t element_index = get_index(matrix_builder->index_hash_map,
					 m_index,n_index);
	if (element_index != no_index)
	{
		matrix_builder->matrix->elements[element_index] = element;
		return;
	}
	if (is_full(matrix_builder))
		expand_sparese_matrix(matrix_builder);
	element_index = matrix_builder->matrix->num_elements++;
	set_index(matrix_builder->index_hash_map,
		  m_index,n_index,element_index);
	matrix_builder->matrix->m_index[element_index] = m_index;
	matrix_builder->matrix->n_index[element_index] = n_index;
	matrix_builder->matrix->elements[element_index] = element;
}

void free_sparse_matrix_builder(sparse_matrix_builder_t matrix_builder)
{
	trim_sparse_matrix(matrix_builder);
	free_index_hash(matrix_builder->index_hash_map);
	free(matrix_builder);
}

static
int is_full(sparse_matrix_builder_t matrix_builder)
{
	return matrix_builder->allocated_num_elements ==
		matrix_builder->matrix->num_elements;
}

static 
void expand_sparese_matrix(sparse_matrix_builder_t matrix_builder)
{
	matrix_builder->allocated_num_elements +=
		matrix_builder->allocated_num_elements +1;
	reallocate_sparse_matrix(matrix_builder);
}

static
void trim_sparse_matrix(sparse_matrix_builder_t matrix_builder)
{
	matrix_builder->allocated_num_elements =
		matrix_builder->matrix->num_elements;
	reallocate_sparse_matrix(matrix_builder);
}

static
void reallocate_sparse_matrix(sparse_matrix_builder_t matrix_builder)
{
	matrix_builder->matrix->m_index = 
		(size_t*)realloc(matrix_builder->matrix->m_index,
				 matrix_builder->allocated_num_elements*
				 sizeof(size_t));
	matrix_builder->matrix->n_index = 
		(size_t*)realloc(matrix_builder->matrix->n_index,
				 matrix_builder->allocated_num_elements*
				 sizeof(size_t));
	matrix_builder->matrix->elements = 
		(double*)realloc(matrix_builder->matrix->elements,
				 matrix_builder->allocated_num_elements*
				 sizeof(double));
}
