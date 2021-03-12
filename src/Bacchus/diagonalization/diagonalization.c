#include <diagonalization/diagonalization.h>
#include <log/log.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

extern void dsteqr_(char *compz,
		    int *matrix_side,
		    double *diagonal,
		    double *off_diagonal,
		    double *eigenvectors,
		    int *leading_dimension,
		    double *work_array,
		    int *info);

eigensystem_t diagonalize_tridiagonal_matrix(
	const double *diagonal_elements,
	const double *off_diagonal_elements,
	const size_t dimension)
{
	char compz = 'V';
	int matrix_side = (int)dimension; 

	eigensystem_t eigensystem = new_empty_eigensystem(dimension);

	double *eigenvalues = get_eigenvalues(eigensystem);
	memcpy(eigenvalues,
		diagonal_elements,
		sizeof(double)*dimension);
	// dsteqr destroys the off diagonal elements, I don't want that
	double *off_diagonal = (double*)calloc(dimension-1,sizeof(double));
	memcpy(off_diagonal,
		off_diagonal_elements,
		sizeof(double)*(dimension-1));		
	
	double *eigenvectors = 
		(double*)calloc(dimension*dimension,
				sizeof(double));		
	for (size_t i = 0; i<dimension; i++)
		eigenvectors[i*dimension+i] = 1;
	double *work_array =
		(double*)calloc(2*dimension-2,
				sizeof(double));
	int info = 0;
	dsteqr_(&compz,
		&matrix_side,
		eigenvalues,
		off_diagonal,
		eigenvectors,
		&matrix_side,
		work_array,
		&info);
	log_entry("info = %d",info);
	assert(info == 0);
	for (size_t i = 0; i<dimension; i++)
		if (eigenvectors[i*dimension]<0)
			for (size_t j = 0; j<dimension; j++)
				eigenvectors[i*dimension+j]*=-1;
	set_eigenvalues(eigensystem,
			eigenvalues);
	set_raw_eigenvectors(eigensystem,
			  eigenvectors);
	free(eigenvalues);
	free(work_array);
	free(eigenvectors);	
	free(off_diagonal);
	return eigensystem;
}

extern void dsyev_(char *jobz,
		char *uplo,
		int *side,
		double *matrix,
		int *lda,
		double *eigenvalues,
		double *work,
		int *lwork,
		int *info);	

eigensystem_t diagonalize_symmetric_matrix(
	matrix_t matrix)
{
	assert(get_num_rows(matrix) == get_num_columns(matrix));
	int side = (int)get_num_rows(matrix);
	double *matrix_elements = get_matrix_elements(matrix);
	
	eigensystem_t eigensystem = new_empty_eigensystem(side);
	double *eigenvalues = get_eigenvalues(eigensystem);
	int lwork = 3*side-1;
	double *work = (double*)calloc(lwork,
			sizeof(double));
	int info = 0;
	dsyev_("N","U",
		&side,
		matrix_elements,
		&side,
		eigenvalues,
		work,
		&lwork,
		&info);
	log_entry("info = %d",info);
	assert(info == 0);
	set_eigenvalues(eigensystem,eigenvalues);
	free(eigenvalues);
	free(matrix_elements);
	free(work);
	return eigensystem;	
}
