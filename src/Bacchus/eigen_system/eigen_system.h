#ifndef __EIGEN_SYSTEM__
#define __EIGEN_SYSTEM__

#include <stdlib.h>

struct _eigen_system_;
typedef struct _eigen_system_ *eigen_system_t;

eigen_system_t new_empty_eigensystem(const size_t num_eigen_values);

size_t get_num_eigen_values(eigen_system_t eigen_system);

double get_eigen_value(eigen_system_t eigen_system,
		size_t eigen_value_index);

double* get_eigen_values(eigen_system_t eigen_system);

void set_eigen_values(eigen_system_t eigen_system,
	double *eigen_values);

void print_eigen_system(const eigen_system_t eigen_system);

void free_eigen_system(eigen_system_t eigen_system);
#endif
