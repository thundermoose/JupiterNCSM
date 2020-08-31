#ifndef __MATRIX_TRANSFORM__
#define __MATRIX_TRANSFORM__
#include <stdlib.h>
#include <stdio.h>
#include <bases/m_scheme_3p_basis.h>
#include <bases/jjj_coupled_3p.h>
#include <clebsch_gordan/clebsch_gordan.h>


typedef struct
{
	double *elements;
	size_t m,n;
} Dens_Matrix;

#define ELEM(densmat,i,j) densmat->elements[i*densmat->n+j]

Dens_Matrix* new_zero_matrix(size_t m,
			     size_t n);

Dens_Matrix* new_identity_matrix(size_t m,
				 size_t n);

void accumulate_matrix(Dens_Matrix mat_acc,
		       Dens_Matrix mat_term);

void to_python_matrix(FILE* file,
		      const char* name,
		      Dens_Matrix mat);
void create_matrix_plot(Dens_Matrix m,
			const char *title);

void print_matrix(Dens_Matrix mat);

double mean_square_difference(Dens_Matrix *matrix_a,
			      Dens_Matrix *matrix_b);

double get_dens_matrix_element(Dens_Matrix *matrix,
			       size_t i,
			       size_t j);

void free_dens_matrix(Dens_Matrix* matrix);

typedef struct _sparse_matrix_{
	size_t m,n;
	size_t num_elements;
	size_t *m_index;
	size_t *n_index;
	double *elements;
} Sparse_Matrix;

// Transform generators

Sparse_Matrix *new_zero_sparse_matrix(size_t m,size_t n);

Sparse_Matrix *identity_transform(size_t n);

Sparse_Matrix *jjj_transform(M_Scheme_3p_Basis* ms3b,
			     JJJ_Basis* jjjbasis,
			     Clebsch_Gordan_Data* cgd);
// Checks if a sparse matrix is unitary
int is_unitary(Sparse_Matrix *spm);

void print_sparse_matrix(Sparse_Matrix mat);
void transpose_sparse_matrix(Sparse_Matrix *matrix);

void sparse_to_python_matrix(FILE* file,
			     const char* name,
			     Sparse_Matrix mat);

Dens_Matrix *sparse_to_dens_matrix(Sparse_Matrix *mat);


void free_sparse_matrix(Sparse_Matrix *matrix);
// Transform applicator 
Dens_Matrix *transform_matrix(Sparse_Matrix *bra_transform,
			      Dens_Matrix *matrix,
			      Sparse_Matrix *ket_transform);



#endif
