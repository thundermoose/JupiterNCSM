#ifndef __MATRIX_BUILDER__
#define __MATRIX_BUILDER__

#include <stdlib.h>
#include <matrix_transform/matrix_transform.h>

struct _sparse_matrix_builder_;
typedef struct _sparse_matrix_builder_ *sparse_matrix_builder_t;

sparse_matrix_builder_t new_sparse_matrix_builder(Sparse_Matrix *matrix);

void set_sparse_matrix_element(sparse_matrix_builder_t matrix_builder,
			     size_t m_index,
			     size_t n_index,
			     double element);

void free_sparse_matrix_builder(sparse_matrix_builder_t matrix_builder);
#endif
