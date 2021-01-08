#ifndef __DIAGONALIZATION__
#define __DIAGONALIZATION__

#include <eigensystem/eigensystem.h>
#include <matrix/matrix.h>
#include <stdlib.h>

eigensystem_t diagonalize_tridiagonal_matrix(
	const double *diagonal_elements,
	const double *off_diagonal_elements,
	size_t dimension);

eigensystem_t diagonalize_symmetric_matrix(
	matrix_t matrix);
#endif
