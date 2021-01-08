#ifndef __NUMPY_MATRIX__
#define __NUMPY_MATRIX__

#include <stdlib.h>

void save_as_numpy_matrix(const char *file_path,
			  double *elements,
			  size_t num_rows,
			  size_t num_columns);

#endif
