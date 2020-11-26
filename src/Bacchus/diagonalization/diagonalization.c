#include <diagonalization/diagonalization.h>
#include <log/log.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

extern void dsteqr_(char *compz,
		    int *matrix_side,
		    double *diagonal,
		    double *off_diagonal,
		    double *eigen_vectors,
		    int *leading_dimension,
		    double *work_array,
		    int *info);

eigen_system_t diagonalize_tridiagonal_matrix(
	const double *diagonal_elements,
	const double *off_diagonal_elements,
	const size_t dimension)
{
	char compz = 'V';
	int matrix_side = (int)dimension; 

	eigen_system_t eigen_system = new_empty_eigensystem(dimension);

	double *eigen_values = get_eigen_values(eigen_system);
	memcpy(eigen_values,
		diagonal_elements,
		sizeof(double)*dimension);
	// dsteqr destroys the off diagonal elements, I don't want that
	double *off_diagonal = (double*)calloc(dimension-1,sizeof(double));
	memcpy(off_diagonal,
		off_diagonal_elements,
		sizeof(double)*(dimension-1));		
	
	double *eigen_vectors = 
		(double*)calloc(dimension*dimension,
				sizeof(double));		
	for (size_t i = 0; i<dimension; i++)
		eigen_vectors[i*dimension+i] = 1;
	double *work_array =
		(double*)calloc(2*dimension-2,
				sizeof(double));
	int info = 0;
	dsteqr_(&compz,
		&matrix_side,
		eigen_values,
		off_diagonal,
		eigen_vectors,
		&matrix_side,
		work_array,
		&info);
	log_entry("info = %d",info);
	assert(info == 0);
	set_eigen_values(eigen_system,
			eigen_values);
	set_raw_eigen_vectors(eigen_system,
			  eigen_vectors);
	free(eigen_values);
	free(work_array);
	free(eigen_vectors);	
	free(off_diagonal);
	return eigen_system;
}

extern void dsyev_(char *jobz,
		char *uplo,
		int *side,
		double *matrix,
		int *lda,
		double *eigen_values,
		double *work,
		int *lwork,
		int *info);	

eigen_system_t diagonalize_symmetric_matrix(
	matrix_t matrix)
{
	assert(get_num_rows(matrix) == get_num_columns(matrix));
	int side = (int)get_num_rows(matrix);
	double *matrix_elements = get_matrix_elements(matrix);
	
	eigen_system_t eigen_system = new_empty_eigensystem(side);
	double *eigen_values = get_eigen_values(eigen_system);
	int lwork = 3*side-1;
	double *work = (double*)calloc(lwork,
			sizeof(double));
	int info = 0;
	dsyev_("N","U",
		&side,
		matrix_elements,
		&side,
		eigen_values,
		work,
		&lwork,
		&info);
	log_entry("info = %d",info);
	assert(info == 0);
	set_eigen_values(eigen_system,eigen_values);
	free(eigen_values);
	free(matrix_elements);
	free(work);
	return eigen_system;	
}
