#include <string.h>
#include <errno.h>
#include <math.h>
#include <jt_transformation/jt_transformation.h>
#include <matrix_builder/matrix_builder.h>
#include <jt_block_iterator/jt_block_iterator.h>
#include <utils/assertion.h>
#include <utils/debug_messages.h>
#include <utils/automatic_test.h>
#include <utils/permutation_tools.h>
#include <thundertester/test.h>

typedef struct
{
	Shells *m_shells;
	SP_States *sp_states;
	strue_shell_index *true_shell_map_jt_to_m;
	Clebsch_Gordan_Data* clebsch_gordan_data;
} jt_workspace_t;

static
double get_jt_to_m_element(m_scheme_2p_state_t m_scheme_state,
			   jt_state_t jt_state,
			   jt_workspace_t *workspace);

static
double normalization(Shell shell_a,
		     Shell shell_b,
		     jt_state_t jt_state,
		     jt_workspace_t *workspace);

static
int is_odd(int number);

Sparse_Matrix* new_jt_transformation_3p(M_Scheme_3p_Basis* m_basis,
					JT_Basis* jt_basis,
					Clebsch_Gordan_Data* cgd)
{
	Sparse_Matrix *final_transform =
		(Sparse_Matrix*)calloc(1,sizeof(Sparse_Matrix));
	ASSERT(final_transform != NULL,
	       return NULL,
	       "calloc could not allocated enough memory\n"
	       "%s\n",strerror(errno));

	ASSERT(m_basis != NULL,
	       return NULL,
	       "null pointer to m-scheme elements given\n");
	ASSERT(jt_basis != NULL,
	       return NULL,
	       "null pointer to jt-scheme elements given\n");
	ASSERT(cgd != NULL,
	       return NULL,
	       "null pointer to clebsch_gordan_data given\n");

	// setting up the basic properties of the m- to j-scheme
	// transformation

	// we have the j-scheme states on the bra side 
	final_transform->m = jt_basis->dimension;

	// we have the m-shceme states on the ket side
	final_transform->n = m_basis->dimension;

	// we know that at most there can be m*n matrix elements
	size_t alloc_num_elements = final_transform->m*final_transform->n;
	// allocating all arrays
	final_transform->m_index = (size_t*)calloc(alloc_num_elements,
						   sizeof(size_t*));
	final_transform->n_index = (size_t*)calloc(alloc_num_elements,
						   sizeof(size_t*));
	final_transform->elements = (double*)calloc(alloc_num_elements,
						    sizeof(double*));

	// keep a pointer to the single particle basis
	// and one to the shells to simplify expressions
	// further down
	SP_States *sp_basis = m_basis->sp_states;
	Shells* shells = sp_basis->shells;


	size_t m_i; // index on the bra side
	size_t n_i; // index on the ket side
	for (m_i = 0;
	     m_i < jt_basis->dimension;
	     m_i++)
	{
		// get the current jt-scheme state
		JT_State jt_state = jt_basis->states[m_i];

		// to simplify future expressions
		quantum_number j_a = shells->true_shells[jt_state.a].j;
		quantum_number j_b = shells->true_shells[jt_state.b].j;
		quantum_number j_c = shells->true_shells[jt_state.c].j;
		quantum_number j_ab = jt_state.jab;
		quantum_number j_abc = jt_state.jabc;

		quantum_number t_ab = jt_state.tab;
		quantum_number t_abc = jt_state.tabc;

		for (n_i = 0;
		     n_i < m_basis->dimension;
		     n_i++)
		{
			// get the current m-scheme state
			M_Scheme_3p_State m_state = m_basis->states[n_i];

			quantum_number m_a = sp_basis->sp_states[m_state.a].m;
			quantum_number m_b = sp_basis->sp_states[m_state.b].m;
			quantum_number m_c = sp_basis->sp_states[m_state.c].m;
			quantum_number m_ab = m_a+m_b;
			quantum_number m_abc = m_ab+m_c;
			DEBUG_MESS("j_ab: %d m_ab: %d j_abc: %d m_abc: %d\n",
				   j_ab,m_ab,j_abc,m_abc);
			// checking the triangular condition
			if (abs(m_abc)>j_abc)
			{
				continue;
			}
			quantum_number tz_a = shells->shells[sp_basis->sp_states[m_state.a].shell].tz;
			quantum_number tz_b = shells->shells[sp_basis->sp_states[m_state.b].shell].tz;
			quantum_number tz_c = shells->shells[sp_basis->sp_states[m_state.c].shell].tz;
			quantum_number tz_ab = tz_a+tz_b;
			quantum_number tz_abc = tz_ab+tz_c;
			// to accumulate all contributions
			double curr_element = 0.0;
			size_t num_contributions = 0;
			// For each posible permutation of
			// the particles in the m_state
			// compute a contribution
			size_t p; // if p is even we have a positive permutation
			// if p is odd we have a negative permutation
			for (p = 0; p<6; p++)
			{
				// check if the permutation is correct
				if (jt_state.a == shells->shells[sp_basis->sp_states[m_state.a].shell].tse &&
				    jt_state.b == shells->shells[sp_basis->sp_states[m_state.b].shell].tse &&
				    jt_state.c == shells->shells[sp_basis->sp_states[m_state.c].shell].tse)
				{
					DEBUG_MESS("m: %ld, n: %ld, p:%ld\n",m_i,n_i,p);
					m_a = sp_basis->sp_states[m_state.a].m;
					m_b = sp_basis->sp_states[m_state.b].m;
					m_c = sp_basis->sp_states[m_state.c].m;
					m_ab = m_a+m_b;
					m_abc = m_ab+m_c;

					tz_a = shells->shells[sp_basis->sp_states[m_state.a].shell].tz;
					tz_b = shells->shells[sp_basis->sp_states[m_state.b].shell].tz;
					tz_c = shells->shells[sp_basis->sp_states[m_state.c].shell].tz;

					tz_ab = tz_a+tz_b;  
					tz_abc = tz_ab+tz_c;



					// compute the two clebsch-gordan coefficients needed
					double cg1 = clebsch_gordan(j_a,j_b,j_ab,
								    m_a,m_b,m_ab,
								    cgd);
					double cg2 = clebsch_gordan(j_ab,j_c,j_abc,
								    m_ab,m_c,m_abc,
								    cgd);
					double cg3 = clebsch_gordan(1,    1,    t_ab,
								    tz_a, tz_b, tz_ab,
								    cgd);
					double cg4 = clebsch_gordan(t_ab,  1,    t_abc,
								    tz_ab, tz_c, tz_abc,
								    cgd);
					curr_element+=cg1*cg2*cg3*cg4*(p&1 ? -1.0 : 1.0);
					num_contributions++;
					break;
				}
				// perform a swap to get to the next
				// permutation of the m-schem state
				if (p&1)
				{
					SWAP(m_state.a,m_state.b);
				}
				else
				{
					SWAP(m_state.b,m_state.c);
				}
			}
			if (num_contributions > 1)
				curr_element/=sqrt(num_contributions);
			DEBUG_MESS("curr_element is %lg\n",curr_element);
			if (fabs(curr_element)>1e-10)
			{
				// Add the current element to our final matrix
				final_transform->m_index[final_transform->num_elements] = m_i;
				final_transform->n_index[final_transform->num_elements] = n_i;
				final_transform->elements[final_transform->num_elements] = curr_element;
				final_transform->num_elements++;
			}
		}
	}

	// our initial guess of the number of non-zero elements are probably wrong
	if (final_transform->num_elements < final_transform->m*final_transform->n)
	{
		alloc_num_elements = final_transform->num_elements;
		final_transform->m_index = (size_t*)realloc(final_transform->m_index,
							    alloc_num_elements*
							    sizeof(size_t*));
		final_transform->n_index = (size_t*)realloc(final_transform->n_index,
							    alloc_num_elements*
							    sizeof(size_t*));
		final_transform->elements = (double*)realloc(final_transform->elements,
							     alloc_num_elements*
							     sizeof(double*));
	}

	return final_transform;
}

Sparse_Matrix* new_jt_transformation(m_scheme_2p_basis_t m_scheme_2p_basis,
				     jt_basis_t jt_basis,
				     Clebsch_Gordan_Data *clebsch_gordan_data)
{
	Shells *jt_shells = get_jt_basis_shells(jt_basis);
	Shells *m_shells = get_m_scheme_shells(m_scheme_2p_basis);
	jt_workspace_t workspace =
	{
		.m_shells = m_shells,
		.sp_states = get_m_scheme_sp_states(m_scheme_2p_basis),
		.true_shell_map_jt_to_m = matching_true_shells(jt_shells,
							       m_shells),
		.clebsch_gordan_data = clebsch_gordan_data

	};
	for(size_t i = 0; i<m_shells->num_of_true_shells; i++)
	{
		DEBUG_MESS("%lu->%lu\n",
			   i,workspace.true_shell_map_jt_to_m[i]);
	}
	const size_t dimension_m_scheme =
		get_m_scheme_2p_dimension(m_scheme_2p_basis);
	const size_t dimension_jt_scheme = get_dimension(jt_basis);
	Sparse_Matrix *transfom = new_zero_sparse_matrix(dimension_jt_scheme,
							 dimension_m_scheme);
	size_t num_allocated_elements = dimension_jt_scheme*dimension_m_scheme;
	transfom->m_index =
		(size_t*)malloc(num_allocated_elements*sizeof(size_t));
	transfom->n_index =
		(size_t*)malloc(num_allocated_elements*sizeof(size_t));
	transfom->elements =
		(double*)malloc(num_allocated_elements*sizeof(double));
	size_t element_index = 0;
	for (size_t i = 0; i < dimension_m_scheme; i++)
	{
		m_scheme_2p_state_t m_scheme_state = 
			get_m_scheme_2p_state(m_scheme_2p_basis,
					      i);
		for (size_t j = 0; j < dimension_jt_scheme; j++)
		{
			jt_state_t jt_state = get_jt_state(jt_basis,j);
			double element = get_jt_to_m_element(m_scheme_state,
							     jt_state,
							     &workspace);			
			if (fabs(element)<1e-10)
				continue;
			transfom->m_index[element_index] = j;
			transfom->n_index[element_index] = i;
			transfom->elements[element_index] = element;
			element_index++;
		}
	}	
	transfom->m_index = (size_t*)realloc(transfom->m_index,
					     element_index*sizeof(size_t));
	transfom->n_index = (size_t*)realloc(transfom->n_index,
					     element_index*sizeof(size_t));
	transfom->elements = (double*)realloc(transfom->elements,
					      element_index*sizeof(double));
	transfom->num_elements = element_index;
	free(workspace.true_shell_map_jt_to_m);
	return transfom;

}

double get_jt_to_m_element(m_scheme_2p_state_t m_scheme_state,
			   jt_state_t jt_state,
			   jt_workspace_t *workspace)
{
	jt_state.a = workspace->true_shell_map_jt_to_m[jt_state.a];
	jt_state.b = workspace->true_shell_map_jt_to_m[jt_state.b];
	SP_State state_a = workspace->sp_states->sp_states[m_scheme_state.a];
	SP_State state_b = workspace->sp_states->sp_states[m_scheme_state.b];
	Shell shell_a = workspace->m_shells->shells[state_a.shell];
	Shell shell_b = workspace->m_shells->shells[state_b.shell];
	int swaped = 0;
	if (jt_state.a == shell_b.tse && jt_state.b == shell_a.tse &&
	   jt_state.a != jt_state.b)
	{
		swaped = 1;
		SP_State tmp_state = state_a;
		state_a = state_b;
		state_b = tmp_state;
		Shell tmp_shell = shell_a;
		shell_a = shell_b;
		shell_b = tmp_shell;
	}	
	int jt_space_symmetric = is_odd((shell_a.j+
					 shell_b.j)/2 - 
					jt_state.J -jt_state.T);
	double normalization_constant = normalization(shell_a,shell_b,
						      jt_state,
						      workspace);
	normalization_constant = jt_space_symmetric && swaped ?
	       	-normalization_constant : normalization_constant;
	if (fabs(normalization_constant) <1e-10)
		return 0.0;	
	double cg1 = clebsch_gordan(shell_a.j,shell_b.j,2*jt_state.J,
				    state_a.m,state_b.m,state_a.m+state_b.m,
				    workspace->clebsch_gordan_data);
	double cg2 = clebsch_gordan(1,1,2*jt_state.T,
				    shell_a.tz,shell_b.tz,
				    shell_a.tz+shell_b.tz,
				    workspace->clebsch_gordan_data);
	return normalization_constant*cg1*cg2;	
}

double normalization(Shell shell_a,
		     Shell shell_b,
		     jt_state_t jt_state,
		     jt_workspace_t *workspace)
{
	int m_shell_a = shell_a.tse;
	int m_shell_b = shell_b.tse;
	if ((jt_state.a != m_shell_a && jt_state.b != m_shell_b &&
	     jt_state.a != m_shell_b && jt_state.b != m_shell_a) ||
	    (jt_state.a != m_shell_a && jt_state.b == m_shell_b) ||
	    (jt_state.a != m_shell_b && jt_state.b == m_shell_a) ||
	    (jt_state.b != m_shell_a && jt_state.a == m_shell_b) ||
	    (jt_state.b != m_shell_b && jt_state.a == m_shell_a))
	{
		return 0.0;
	}
	else if (m_shell_a != m_shell_b)
	{
		return 1;
	}
	else if (is_odd(jt_state.J+jt_state.T) && jt_state.a == jt_state.b)
	{
		return sqrt(2.0);
	}
	else
	{
		return 0.0;
	}
}

int is_odd(int number)
{
	return number % 2 == 1;
}


new_test(jt_transform_3p_unitarity_test,
	 {
	 const quantum_number nmax = 4;
	 // Setup the shells and the m-scheme basis
	 Shells *shells = new_shells(nmax);
	 SP_States *sp_states = new_sp_states(shells);
	 M_Scheme_3p_Basis *m_scheme_2p_basis = new_m_scheme_3p_basis3(nmax,1,1,sp_states);
	 /*ASSERT(m_scheme_2p_basis->dimension == 1,
	   TEST_FAILED("m_scheme_2p_basis = %ld, but should be 1\n",
	   m_scheme_2p_basis->dimension),
	   "");*/
	 JT_Basis* jt_basis = new_jt_basis_from_m_scheme(m_scheme_2p_basis);

	 assert_that(jt_basis != NULL);

	 assert_that(!contains_duplicates(jt_basis));
	 DEBUG_MESS("jt_basis->dimension = %ld\n",jt_basis->dimension);



	 Clebsch_Gordan_Data *cgd = initiate_clebsch_gordan(3*(2*nmax+1));

	 Sparse_Matrix* transform =
		 new_jt_transformation_3p(m_scheme_2p_basis,
					  jt_basis,cgd);
	 free_clebsch_gordan(cgd);

	 DEBUG_MESS("transform->size = (%ld %ld)\n",
		    transform->m,transform->n);

	 double *acc = (double*)calloc(transform->n,
				       sizeof(double));
	 size_t i;
	 for (i = 0; i<transform->num_elements; i++)
	 {
		 acc[transform->n_index[i]]+=transform->elements[i]*transform->elements[i];
	 }
	 //DEBUG_MESS("acc = %lg\n",acc);
	 print_jt_basis_3p(jt_basis);
	 list_m_scheme_3p_basis(m_scheme_2p_basis);
	 for (i = 0; i<transform->n; i++)
	 {
		 assert_that(fabs(acc[i]-1)<1e-10);
	 }
	 free(acc);
	 free_sparse_matrix(transform);
	 free_jt_basis_3p(jt_basis);
	 free_m_scheme_3p_basis(m_scheme_2p_basis);
	 free_sp_states(sp_states);
	 free_shells(shells);
	 });

new_test(jt_transformation_unitarity,
	 {
	 const quantum_number e_max1 = 2;
	 const quantum_number e_max2 = 2;
	 const double tollerance = 1e-6;
	 const char *transform_filename =
	 get_test_file_path("transform.py");
	 jt_basis_t jt_basis = new_antoine_basis(e_max1,
						 e_max2);
	 print_jt_basis(jt_basis);
	 m_scheme_2p_basis_t m_scheme_2p_basis = 
	 new_m_scheme_2p_basis_fixed_isospin(e_max1,
					     e_max2,0);
	 print_m_scheme_2p_basis(m_scheme_2p_basis);
	 assert_that(get_dimension(jt_basis) ==
		     get_m_scheme_2p_dimension(m_scheme_2p_basis));
	 Clebsch_Gordan_Data *cgd = initiate_clebsch_gordan(2*(2*e_max2+1));
	 Sparse_Matrix *transform = new_jt_transformation(m_scheme_2p_basis,
							  jt_basis,
							  cgd);
	 free_clebsch_gordan(cgd);
	 //	transpose_sparse_matrix(transform);
	 Dens_Matrix *identity_matrix =
		 new_identity_matrix(get_dimension(jt_basis),
				     get_dimension(jt_basis));
	 Dens_Matrix *transformed_matrix = transform_matrix(transform,
							    identity_matrix,
							    transform);
	 FILE *transform_file = fopen(transform_filename,"w");
	 if (transform_file == NULL)
	 {
		 printf("filename = %s\n",transform_filename);
		 printf("ERROR_MESSAGE: %s\n",strerror(errno));
	 }
	 assert_that(transform_file != NULL);
	 sparse_to_python_matrix(transform_file,"transform",*transform);
	 print_matrix(*transformed_matrix);
	 assert_that(mean_square_difference(identity_matrix,
					    transformed_matrix) < tollerance);	
	 }
);
