#include "matrix_transform.h"
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <utils/helpful_macros.h>
#include <utils/permutation_tools.h>
#include <utils/assertion.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>

Dens_Matrix* new_zero_matrix(size_t m,
			     size_t n){
	Dens_Matrix* out_matrix = (Dens_Matrix*)malloc(sizeof(Dens_Matrix));
	out_matrix->m = m;
	out_matrix->n = n;
	out_matrix->elements = (double*)calloc(m*n,sizeof(double));
	return out_matrix;
}

Dens_Matrix* new_identity_matrix(size_t m,
				 size_t n)
{
	Dens_Matrix* matrix = new_zero_matrix(m,n);
	for (size_t i = 0; i<min(m,n); i++)
		matrix->elements[i*(n+1)] = 1;
	return matrix;
}

void accumulate_matrix(Dens_Matrix mat_acc,
		       Dens_Matrix mat_term){
	if (__builtin_expect(mat_acc.m != mat_term.m ||
			     mat_acc.n != mat_term.n,0)){
		fprintf(stderr,"Matrices need to have the same size to be accumulated\n");
		exit(1);
	}

	size_t i,j;
	for (i = 0; i<mat_acc.m; i++){
		for (j = 0; j<mat_acc.n; j++){
#pragma omp atomic
			mat_acc.elements[i*mat_acc.n+j] += mat_term.elements[i*mat_acc.n+j];
		}
	}

}

void to_python_matrix(FILE* file,
		      const char* name,
		      Dens_Matrix mat){
	size_t i,j;
	fprintf(file,"%s = np.array([",name);
	for (i = 0; i<mat.m; i++){
		fprintf(file,"%s[",i>0 ? " ," : "");
		for (j = 0; j<mat.n; j++){
			fprintf(file,"%s%lg",j>0 ? " ," : "",mat.elements[i*mat.n+j]);
		}
		fprintf(file,"]");
	}
	fprintf(file,"])\n");
}

void create_matrix_plot(Dens_Matrix m,
			const char *title)
{
	FILE* python = popen("python","w");
	fprintf(python,
		"import numpy as np\n"
		"import matplotlib.pyplot as plt\n");
	to_python_matrix(python,
			 "m",
			 m);
	fprintf(python,
		"plt.matshow(m)\n"
		"plt.title('%s')\n"
		"plt.show()\n"
		"quit()\n",
		title);
	pclose(python);
}

void print_matrix(Dens_Matrix mat){
	size_t i,j;
	for (i = 0; i<mat.m; i++){
		for (j = 0; j<mat.n; j++){
			printf(" %5.2lf",mat.elements[i*mat.n+j]);
		}
		printf("\n");
	}
}

double mean_square_difference(Dens_Matrix *matrix_a,
			      Dens_Matrix *matrix_b)
{
	double accumulator = 0;
	for (size_t i = 0; i<matrix_a->m*matrix_a->n; i++)
	{
		double difference =
			(matrix_a->elements[i]-matrix_b->elements[i]);
		accumulator += difference*difference;
	}
	return accumulator/(matrix_a->m*matrix_a->n);
}

double get_dens_matrix_element(Dens_Matrix *matrix,
			       size_t i,
			       size_t j)
{
	return matrix->elements[i*matrix->n+j];
}

void free_dens_matrix(Dens_Matrix* matrix){
	free(matrix->elements);
	free(matrix);
}

Sparse_Matrix *new_zero_sparse_matrix(size_t m,size_t n)
{
	Sparse_Matrix *matrix =
		(Sparse_Matrix*)calloc(1,sizeof(Sparse_Matrix));
	matrix->m = m;
	matrix->n = n;
	return matrix;
}

Sparse_Matrix *identity_transform(size_t n){
	Sparse_Matrix *out_matrix = (Sparse_Matrix*)malloc(sizeof(Sparse_Matrix));
	out_matrix->n = n;
	out_matrix->m = n;
	out_matrix->num_elements = n;
	out_matrix->elements = (double*)malloc(n*sizeof(double));
	out_matrix->m_index = (size_t*)malloc(n*sizeof(size_t));
	out_matrix->n_index = (size_t*)malloc(n*sizeof(size_t));
	size_t i;
	//#pragma parallel for private(i)
	for (i = 0; i<n; i++){
		out_matrix->m_index[i] = i;
		out_matrix->n_index[i] = i;
		out_matrix->elements[i] = 1.0;
	}
	return out_matrix;
}


/* This function creates a transformation matrix from
 * a m-scheme basis to a j-scheme basis. It assumes
 * that the shells that the particles occupies are
 * ordered such that a<=b, and c is free.
 */
Sparse_Matrix *jjj_transform(M_Scheme_3p_Basis* ms3b,
			     JJJ_Basis* jjjbasis,
			     Clebsch_Gordan_Data* cgd){
	// here we set up the basic structure of the
	// output matrix,
	Sparse_Matrix *out_matrix = (Sparse_Matrix*)malloc(sizeof(Sparse_Matrix));
	out_matrix->n = ms3b->dimension;
	out_matrix->m = jjjbasis->dimension;
	out_matrix->num_elements = 0;
	size_t max_num_elements = 0;
	out_matrix->m_index = NULL;
	out_matrix->n_index = NULL;
	out_matrix->elements = NULL;
	// These loops are for looping over
	size_t i;
	size_t j;
	//printf("Inside jjj_transform\n");
	for (i = 0; i<out_matrix->m; i++){
		JJJ_State m_state = jjjbasis->states[i];
		for (j = 0; j<out_matrix->n; j++){
			M_Scheme_3p_State n_state = ms3b->states[j];
			int sign = 1;
			quantum_number ma,mb,mc;
			ma = ms3b->sp_states->sp_states[n_state.a].m;
			mb = ms3b->sp_states->sp_states[n_state.b].m;
			mc = ms3b->sp_states->sp_states[n_state.c].m;


			if (m_state.a == ms3b->sp_states->sp_states[n_state.a].shell &&
			    m_state.b == ms3b->sp_states->sp_states[n_state.b].shell &&
			    m_state.c == ms3b->sp_states->sp_states[n_state.c].shell)
			{
				// Nothing should happen
			}
			else if (m_state.a == ms3b->sp_states->sp_states[n_state.b].shell &&
				 m_state.b == ms3b->sp_states->sp_states[n_state.a].shell &&
				 m_state.c == ms3b->sp_states->sp_states[n_state.c].shell)
			{
				SWAP(ma,mb);
				sign=-1;
			}
			else if (m_state.a == ms3b->sp_states->sp_states[n_state.a].shell &&
				 m_state.b == ms3b->sp_states->sp_states[n_state.c].shell &&
				 m_state.c == ms3b->sp_states->sp_states[n_state.b].shell)
			{
				SWAP(mb,mc);
				sign=-1;
			}
			else if (m_state.a == ms3b->sp_states->sp_states[n_state.b].shell &&
				 m_state.b == ms3b->sp_states->sp_states[n_state.c].shell &&
				 m_state.c == ms3b->sp_states->sp_states[n_state.a].shell)
			{
				SWAP(ma,mb);
				SWAP(mb,mc);
				sign=1;
			}
			else if (m_state.a == ms3b->sp_states->sp_states[n_state.c].shell &&
				 m_state.b == ms3b->sp_states->sp_states[n_state.a].shell &&
				 m_state.c == ms3b->sp_states->sp_states[n_state.b].shell)
			{
				SWAP(mb,mc);
				SWAP(ma,mb);
				sign=1;
			}
			else if (m_state.a == ms3b->sp_states->sp_states[n_state.c].shell &&
				 m_state.b == ms3b->sp_states->sp_states[n_state.b].shell &&
				 m_state.c == ms3b->sp_states->sp_states[n_state.a].shell)
			{

				SWAP(mb,mc);
				SWAP(ma,mb);
				SWAP(mb,mc);
				sign=-1;
			}
			else
			{
				continue;
			}
			ASSERT(m_state.a<=m_state.b,
			       exit(1),
			       "m_state.a > m_state.b\n");

			log_entry("computing: <((%ld %ld)%d %ld)%d|(%ld %d) (%ld %d) (%ld %d)>\n",
				   m_state.a,m_state.b,m_state.j_ab,m_state.c,m_state.j_abc,
				   ms3b->sp_states->sp_states[n_state.a].shell,
				   ms3b->sp_states->sp_states[n_state.a].m,
				   ms3b->sp_states->sp_states[n_state.b].shell,
				   ms3b->sp_states->sp_states[n_state.b].m,
				   ms3b->sp_states->sp_states[n_state.c].shell,
				   ms3b->sp_states->sp_states[n_state.c].m);


			quantum_number ja,jb,jc;
			ja = jjjbasis->shells->shells[m_state.a].j;
			jb = jjjbasis->shells->shells[m_state.b].j;
			jc = jjjbasis->shells->shells[m_state.c].j;
			quantum_number M_tot = ma+mb+mc;
			if (m_state.j_abc<abs(M_tot)){
				log_entry("e = %4.3lf\n",0.0);
				continue;
			}

			double cg1,cg2;
			cg1 = clebsch_gordan(ja,jb,m_state.j_ab,
					     ma,mb,ma+mb,
					     cgd);
			log_entry("<%d %d,%d %d|%d %d> = %lf\n",ja,ma,jb,mb,m_state.j_ab,ma+mb,cg1);
			cg2 = clebsch_gordan(m_state.j_ab,jc,m_state.j_abc,
					     ma+mb,mc,M_tot,
					     cgd);
			log_entry("<%d %d,%d %d|%d %d> = %lf\n",m_state.j_ab,ma+mb,jc,mc,m_state.j_abc,ma+mb+mc,cg2);
			if (max_num_elements == out_matrix->num_elements){
				max_num_elements = max_num_elements*2 + 1;
				out_matrix->m_index = (size_t*)realloc(out_matrix->m_index,
								       sizeof(size_t)*max_num_elements);
				out_matrix->n_index = (size_t*)realloc(out_matrix->n_index,
								       sizeof(size_t)*max_num_elements);
				out_matrix->elements = (double*)realloc(out_matrix->elements,
									sizeof(double)*max_num_elements);
			}
			out_matrix->m_index[out_matrix->num_elements] = i;
			out_matrix->n_index[out_matrix->num_elements] = j;
			out_matrix->elements[out_matrix->num_elements++] = cg1*cg2*sign;
			log_entry("cg1 * cg2 = %4.3lf * %4.3lf = %4.3lf \n",cg1,cg2,cg1*cg2);
		}
		log_entry("\n");
	}
	out_matrix->m_index = (size_t*)realloc(out_matrix->m_index,
					       sizeof(size_t)*out_matrix->num_elements);
	out_matrix->n_index = (size_t*)realloc(out_matrix->n_index,
					       sizeof(size_t)*out_matrix->num_elements);
	out_matrix->elements = (double*)realloc(out_matrix->elements,
						sizeof(double)*out_matrix->num_elements);
	return out_matrix;
}


int is_unitary(Sparse_Matrix *spm)
{
	// Declaring general purpose index variables
	size_t i,j;
	// Makeing sure we got correct input
	ASSERT(spm != NULL,
	       exit(1),
	       "Can't check if null pointer spm is unitary\n");

	// An empty sparse matrix is a matrix filed with zeros
	// and therefore can not be unitary
	if (spm->num_elements == 0)
	{
		return 0;
	}

	// checks if M^TM == I (I is nxn)
	Dens_Matrix *res = new_zero_matrix(spm->n,
					   spm->n);
	ASSERT(res != NULL,
	       exit(1),
	       "Couldn't create %d x %d zero dens matrix needed for unitarity test\n",
	       spm->n,spm->n);
	// now perform the sparse sparse matrix multiplication
	for (i = 0; i<spm->num_elements; i++)
	{
		for (j= 0; j<spm->num_elements; j++)
		{
			if (spm->m_index[i] == spm->m_index[j])
			{
				log_entry("M^TM_{%ld %ld} += %1.17lg {%ld}\n",
					   spm->n_index[i],spm->n_index[j],
					   spm->elements[i]*spm->elements[j],
					   spm->m_index[i]);
				ELEM(res,spm->n_index[i],spm->n_index[j])+=spm->elements[i]*spm->elements[j];
			}
		}
	}

	// check if res is a unitary matrix
	for (i = 0; i<res->m; i++)
	{
		for (j = 0; j<res->n; j++)
		{
			if (fabs(ELEM(res,i,j) - (i == j ? 1.0 : 0.0))>1e-10)
			{
				fprintf(stderr,"M^TM_{%ld %ld} is %lg\n",
					i,j,ELEM(res,i,j));
				return 0;
			}
		}
	}
	free_dens_matrix(res);
	// if we get this far M^TM = I

	// now we check if MM^T = I (I is a m x m matrix)
	res = new_zero_matrix(spm->m,spm->m);
	ASSERT(res != NULL,
	       exit(1),
	       "Couldn't create %ld x %ld zero dens matrix needed for unitarity test\n",
	       spm->m,spm->m);
	// now perform the sparse sparse matrix multiplication
	for (i = 0; i<spm->num_elements; i++)
	{
		for (j= 0; j<spm->num_elements; j++)
		{
			if (spm->n_index[i] == spm->n_index[j])
			{
				ELEM(res,spm->m_index[i],spm->m_index[j])+=spm->elements[i]*spm->elements[j];
			}
		}
	}

	// check if res is a unitary matrix
	for (i = 0; i<res->m; i++)
	{
		for (j = 0; j<res->n; j++)
		{
			if (fabs(ELEM(res,i,j) - (i == j ? 1.0 : 0.0))>1e-10)
			{
				fprintf(stderr,"MM^T_{%ld %ld} is %lg\n",
					i,j,ELEM(res,i,j));
				return 0;
			}
		}
	}
	free_dens_matrix(res);
	// if we get this far we know for certain that spm (aka M) is unitary
	return 1;
}


void transpose_sparse_matrix(Sparse_Matrix *matrix)
{
	size_t *tmp_index = matrix->m_index;
	matrix->m_index = matrix->n_index;
	matrix->n_index = tmp_index;
}

void print_sparse_matrix(Sparse_Matrix mat)
{
	Dens_Matrix *dm = new_zero_matrix(mat.m,mat.n);
	size_t i;
	for (i = 0; i<mat.num_elements; i++)
	{
		ELEM(dm,mat.m_index[i],mat.n_index[i]) = mat.elements[i];
	}
	print_matrix(*dm);
	free_dens_matrix(dm);
}

void sparse_to_python_matrix(FILE* file,
			     const char* name,
			     Sparse_Matrix mat)
{
	fprintf(file,
		"import numpy as np\n");
	fprintf(file,
		"%s = np.zeros((%ld, %ld))\n",
		name,mat.m,mat.n);
	size_t i;
	for (i = 0; i<mat.num_elements; i++)
	{
		fprintf(file,
			"%s[%ld,%ld] = %lg\n",
			name,
			mat.m_index[i],
			mat.n_index[i],
			mat.elements[i]);
	}
}

Dens_Matrix *sparse_to_dens_matrix(Sparse_Matrix *mat)
{
	Dens_Matrix *out = new_zero_matrix(mat->m,
					   mat->n);
	size_t i;
	for (i = 0; i<mat->num_elements; i++)
	{
		ELEM(out,mat->m_index[i],mat->n_index[i]) = mat->elements[i];
	}
	return out;
}

void free_sparse_matrix(Sparse_Matrix *matrix)
{
	free(matrix->m_index);
	free(matrix->n_index);
	free(matrix->elements);
	free(matrix);
}

Dens_Matrix *transform_matrix(Sparse_Matrix *bra_transform,
			      Dens_Matrix *matrix,
			      Sparse_Matrix *ket_transform)
{
	if (matrix->n != ket_transform->m)
	{
		fprintf(stderr,"The matrix to transform has to have n equal to m of the ket_transform, n = %ld, m = %ld\n",
			matrix->n,ket_transform->m);
		exit(1);
	}

	if (matrix->m != bra_transform->m)
	{
		fprintf(stderr,"The matrix to transform has to have m equal to m of the bra_transform, m_d = %ld, m_b = %ld\n",
			matrix->m,bra_transform->m);
		exit(1);
	}

	Dens_Matrix *intermediate = new_zero_matrix(matrix->m,
						    ket_transform->n);

	// D*S
	size_t i,t,j,k;
	//#pragma omp parallel for private(i,t,j,k)
	for (i=0; i<matrix->m; i++)
	{
		double* row_out = intermediate->elements+i*ket_transform->n;
		double* row_in = matrix->elements+i*ket_transform->m;
		for (t = 0; t<ket_transform->num_elements; t++)
		{
			k = ket_transform->m_index[t];
			j = ket_transform->n_index[t];
			row_out[j]+=ket_transform->elements[t]*row_in[k];
		}

	}

	Dens_Matrix* out_matrix = new_zero_matrix(bra_transform->n,
						  ket_transform->n);
	// S^t*D'
	//#pragma parallel for private(t,i,j,k)
	for (t = 0; t<bra_transform->num_elements; t++)
	{
		i = bra_transform->n_index[t];
		k = bra_transform->m_index[t];
		double e = bra_transform->elements[t];
		double *row_in = intermediate->elements+k*intermediate->n;
		double *row_out = out_matrix->elements+i*intermediate->n;
		for (j = 0; j<intermediate->n; j++)
		{
			row_out[j]+=row_in[j]*e;
		}
	}
	free_dens_matrix(intermediate);
	return out_matrix;
}
