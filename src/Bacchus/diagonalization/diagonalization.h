#ifndef __DIAGONALIZATION__
#define __DIAGONALIZATION__

#include <eigen_system/eigen_system.h>
#include <matrix/matrix.h>
#include <stdlib.h>

eigen_system_t diagonalize_tridiagonal_matrix(
	const double *diagonal_elements,
	const double *off_diagonal_elements,
	size_t dimension);

eigen_system_t diagonalize_symmetric_matrix(
	matrix_t matrix);
#endif
