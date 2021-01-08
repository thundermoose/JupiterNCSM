#ifndef __NORM_MATRIX__
#define __NORM_MATRIX__

#include <settings/settings.h>
#include <vector/vector.h>

void create_norm_matrix(const char *norm_matrix_path,
			vector_t *training_vectors,
			size_t num_training_vectors);
#endif
