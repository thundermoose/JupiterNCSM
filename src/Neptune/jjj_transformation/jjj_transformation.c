#include <string.h>
#include <errno.h>
#include <math.h>
#include "jjj_transformation.h"
#include "../utils/assertion.h"
#include "../utils/debug_messages.h"
#include "../utils/automatic_test.h"
#include "../utils/permutation_tools.h"
#include <debug_mode/debug_mode.h>


Sparse_Matrix* new_jjj_transformation(M_Scheme_3p_Basis* m_basis,
				      JJJ_Basis* j_basis,
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
  ASSERT(j_basis != NULL,
	 return NULL,
	 "null pointer to j-scheme elements given\n");
  ASSERT(cgd != NULL,
	 return NULL,
	 "null pointer to clebsch_gordan_data given\n");

  // setting up the basic properties of the m- to j-scheme
  // transformation
  
  // we have the j-scheme states on the bra side 
  final_transform->m = j_basis->dimension;

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
       m_i < j_basis->dimension;
       m_i++)
    {
      // get the current j-scheme state
      JJJ_State j_state = j_basis->states[m_i];

      // to simplify future expressions
      quantum_number j_a = shells->shells[j_state.a].j;
      quantum_number j_b = shells->shells[j_state.b].j;
      quantum_number j_c = shells->shells[j_state.c].j;
      quantum_number j_ab = j_state.j_ab;
      quantum_number j_abc = j_state.j_abc;
      
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
	      if (j_state.a == sp_basis->sp_states[m_state.a].shell &&
		  j_state.b == sp_basis->sp_states[m_state.b].shell &&
		  j_state.c == sp_basis->sp_states[m_state.c].shell)
		{
		  DEBUG_MESS("m: %ld, n: %ld, p:%ld\n",m_i,n_i,p);
		  m_a = sp_basis->sp_states[m_state.a].m;
		  m_b = sp_basis->sp_states[m_state.b].m;
		  m_c = sp_basis->sp_states[m_state.c].m;
		  m_ab = m_a+m_b;
		  m_abc = m_ab+m_c;
		  // compute the two clebsch-gordan coefficients needed
		  double cg1 = clebsch_gordan(j_a,j_b,j_ab,
					      m_a,m_b,m_ab,
					      cgd);
		  double cg2 = clebsch_gordan(j_ab,j_c,j_abc,
					      m_ab,m_c,m_abc,
					      cgd);
		  
		  curr_element+=cg1*cg2*(p&1 ? -1.0 : 1.0);
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

BEGIN_TEST(simple_jjj_transformation)
// Setting up the m-scheme basis, consisting of one
// Slater Determinant
M_Scheme_3p_Basis ms3b;
ms3b.dimension = 1;
M_Scheme_3p_State phi = {0,1,2};
ms3b.states = &phi;
Shells *shells = new_shells(0);
ms3b.sp_states = new_sp_states(shells);
ms3b.e_max = 0;

list_m_scheme_3p_basis(&ms3b);

// sett up the j-scheme basis
JJJ_Basis jjj_basis; 

JJJ_State psi1 = {0, 1, 0, // a, b, c
		  0, 1, // j_ab, j_abc,
		  -1, // tz
		  0}; // parity
JJJ_State psi2 = {0, 1, 0, // a, b, c
		  2, 1, // j_ab, j_abc,
		  -1, // tz
		  0}; // parity
JJJ_State psi3 = {0, 1, 0, // a, b, c
		  2, 3, // j_ab, j_abc,
		  -1, // tz
		  0}; // parity

jjj_basis.dimension = 3;
jjj_basis.states = (JJJ_State*)calloc(3,sizeof(JJJ_State));
jjj_basis.states[0] = psi1;
jjj_basis.states[1] = psi2;
jjj_basis.states[2] = psi3;
jjj_basis.shells = shells;

list_jjj_states(&jjj_basis);

// Checking that we get the j-scheme basis we expect
ASSERT(jjj_basis.dimension == 3,
       TEST_FAILED("jjj_basis.dimension = %ld, but should be 3\n",
		   jjj_basis.dimension),
       "");

ASSERT(jjj_basis.states[0].a == 0 &&
       jjj_basis.states[0].b == 1 &&
       jjj_basis.states[0].c == 0 &&
       jjj_basis.states[0].j_ab == 0 &&
       jjj_basis.states[0].j_abc == 1,
       TEST_FAILED("The first element in jjj_basis should be |((0 1)0 0) 1>\n"
		   "but is |((%ld %ld)%d %ld)%d>\n",
		   jjj_basis.states[0].a,
		   jjj_basis.states[0].b,
		   jjj_basis.states[0].j_ab,
		   jjj_basis.states[0].c,
		   jjj_basis.states[0].j_abc),
       "");

ASSERT(jjj_basis.states[1].a == 0 &&
       jjj_basis.states[1].b == 1 &&
       jjj_basis.states[1].c == 0 &&
       jjj_basis.states[1].j_ab == 2 &&
       jjj_basis.states[1].j_abc == 1,
       TEST_FAILED("The second element in jjj_basis should be |((0 1)2 0) 1>\n"
		   "but is |((%ld %ld)%d %ld)%d>\n",
		   jjj_basis.states[1].a,
		   jjj_basis.states[1].b,
		   jjj_basis.states[1].j_ab,
		   jjj_basis.states[1].c,
		   jjj_basis.states[1].j_abc),
       "");

ASSERT(jjj_basis.states[2].a == 0 &&
       jjj_basis.states[2].b == 1 &&
       jjj_basis.states[2].c == 0 &&
       jjj_basis.states[2].j_ab == 2 &&
       jjj_basis.states[2].j_abc == 3,
       TEST_FAILED("The third element in jjj_basis should be |((0 1)2 0) 3>\n"
		   "but is |((%ld %ld)%d %ld)%d>\n",
		   jjj_basis.states[1].a,
		   jjj_basis.states[1].b,
		   jjj_basis.states[1].j_ab,
		   jjj_basis.states[1].c,
		   jjj_basis.states[1].j_abc),
       "");

// Set up the clebsch-gordan computation scheme
Clebsch_Gordan_Data* cgd = initiate_clebsch_gordan(3); // 3 is the maximum j that we have
// Generates the transformation
Sparse_Matrix* transform = new_jjj_transformation(&ms3b,
						  &jjj_basis,
						  cgd);
free_clebsch_gordan(cgd);
ASSERT(transform != NULL,
       TEST_FAILED("new_jjj_transformation gave a null pointer\n"),
       "");

ASSERT(transform->m == 3,
       TEST_FAILED("transform->m = %ld, but should equal 3\n",transform->m),
       "");

ASSERT(transform->n == 1,
       TEST_FAILED("transform->n = %ld, but should equal 1\n",transform->n),
       "");

ASSERT(transform->num_elements == 2, // There are two non zero matrix elements
       TEST_FAILED("transform->num_elements is %ld but should be 2\n",
		   transform->num_elements),
       "");

ASSERT(transform->m_index != NULL,
       TEST_FAILED("transform->m_index is NULL\n"),
       "");
ASSERT(transform->n_index != NULL,
       TEST_FAILED("transform->n_index is NULL\n"),
       "");
ASSERT(transform->elements != NULL,
       TEST_FAILED("transform->elements is NULL\n"),
       "");

// Print the non zeros matrix elements
printf("The matrix elements:\n");
size_t i;
for (i = 0; i<transform->num_elements; i++)
  {
    printf("M(%ld,%ld) = %1.17lg\n",
	   transform->m_index[i],
	   transform->n_index[i],
	   transform->elements[i]);
  }


// Checking if the matrix column is unitary
double accumulator = 0.0;
for (i = 0; i<transform->num_elements; i++)
  {
    accumulator+=transform->elements[i]*transform->elements[i];
  }

ASSERT(fabs(accumulator-1.0)<1e-10,
       TEST_FAILED("The single column is not unitary\n"),
       "");
free_sparse_matrix(transform);

free(jjj_basis.states);

free_sp_states(ms3b.sp_states);
free_shells(shells);

END_TEST


// Testing transforming between m-scheme with Nmax = 0 to corresponding j-scheme
BEGIN_TEST(m_to_j_nmax0_transfom);
Shells *shells = new_shells(0);
SP_States *sp_basis = new_sp_states(shells);
M_Scheme_3p_Basis *m_basis =
  new_m_scheme_3p_basis_no_m_rest(0,sp_basis);

list_m_scheme_3p_basis(m_basis);

JJJ_Basis *j_basis =
  new_jjj_basis_m_scheme(m_basis);

list_jjj_states(j_basis);
// Setting up
Clebsch_Gordan_Data* cgd = initiate_clebsch_gordan(3);
// Generates the transformation
Sparse_Matrix* transform = new_jjj_transformation(m_basis,
						  j_basis,
						  cgd);
free_clebsch_gordan(cgd);
printf("The matrix elements:\n");


size_t i;
for (i = 0; i<transform->num_elements; i++)
  {
    printf("M(%ld,%ld) = %1.17lg\n",
	   transform->m_index[i],
	   transform->n_index[i],
	   transform->elements[i]);
  }

printf("In matrix form:\n");
print_sparse_matrix(*transform);

free_jjj_basis(j_basis);
free_sparse_matrix(transform);
free_m_scheme_3p_basis(m_basis);
free_sp_states(sp_basis);
free_shells(shells);

END_TEST
     
