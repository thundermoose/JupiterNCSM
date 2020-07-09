#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include "jjj_coupled_3p.h"
#include <string_tools/string_tools.h>
#include <utils/permutation_tools.h>
#include <utils/automatic_test.h>
#include <utils/assertion.h>
#include <utils/debug_messages.h>


#define MAX(a,b) ((a)<(b) ? (b) : (a))


/* Compares two JJJ_States
 * it starts by comparing Jabc
 * if Jabc is equal it compars Tz
 * if Tz are equal it compairs
 * parity
 * and lastly it compares 
 * total harmonic oscillator energy
 */
int comp_jjj_states(const void *a,
		    const void *b,
		    void *sh)
{
  Shells* shells = (Shells*)sh;
  JJJ_State* state_a = (JJJ_State*)a;
  JJJ_State* state_b = (JJJ_State*)b;
  
  if (state_a->j_abc-state_b->j_abc)
    {
      return state_a->j_abc-state_b->j_abc;
    }

  if (state_a->j_ab-state_b->j_ab)
    {
      return state_a->j_ab-state_b->j_ab;
    }
  quantum_number Tza =
    (shells->shells[state_a->a].tz+
     shells->shells[state_a->b].tz+
     shells->shells[state_a->c].tz);
  quantum_number Tzb =
    (shells->shells[state_b->a].tz+
     shells->shells[state_b->b].tz+
     shells->shells[state_b->c].tz);
  if (Tza-Tzb)
    {
      return Tza-Tzb;
    }
  quantum_number pa =
    (shells->shells[state_a->a].l+
     shells->shells[state_a->b].l+
     shells->shells[state_a->c].l)&1;
  quantum_number pb =
    (shells->shells[state_b->a].l+
     shells->shells[state_b->b].l+
     shells->shells[state_b->c].l)&1;
  if (pa-pb)
    {
      return pa-pb;
    }
  quantum_number Ea =
    (shells->shells[state_a->a].e+
     shells->shells[state_a->b].e+
     shells->shells[state_a->c].e);
  quantum_number Eb =
    (shells->shells[state_b->a].e+
     shells->shells[state_b->b].e+
     shells->shells[state_b->c].e);
  //if (Ea-Eb)
  return Ea-Eb;
}

/* The following function generates an
 * List of all possible jjj-coupled
 * three-body states
 */
JJJ_Basis* new_jjj_basis(quantum_number e_max,
			 Shells* shells)
{
  JJJ_Basis* jjj_basis = (JJJ_Basis*)malloc(sizeof(JJJ_Basis));
  jjj_basis->shells = shells;
  // We do not yet know
  // the final number
  // of coupled three-particle
  // states
  jjj_basis->states = NULL;
  jjj_basis->dimension = 0;
  size_t max_dimension = 0;
  JJJ_State jjj_state;
  // The first single particle-state can reside in
  // any shell
  for (jjj_state.a=0;
       jjj_state.a<shells->num_of_shells;
       jjj_state.a++)
    {
      Shell sh_a = shells->shells[jjj_state.a];
      if (sh_a.e>e_max)
	break;
      // We follow the convention by
      // Gustav Jansen that the shell
      // of the second single-particle state
      // can only reside in shells above of
      // equal to the first shell
      for (jjj_state.b=jjj_state.a;
	   jjj_state.b<shells->num_of_shells;
	   jjj_state.b++)
	{
	  Shell sh_b = shells->shells[jjj_state.b];
	  if (sh_a.e+sh_b.e>e_max)
	    break;
	  quantum_number j_ab_min = abs(sh_a.j-sh_b.j);
	  quantum_number j_ab_max = sh_a.j+sh_b.j;
	  // The third single-particle state
	  // are allowed to be in any shall
	  for (jjj_state.c = 0;
	       jjj_state.c<=shells->num_of_shells;
	       jjj_state.c++)
	    {
	      Shell sh_c = shells->shells[jjj_state.c];
	      if (sh_a.e+sh_b.e+sh_c.e>e_max)
		break;
	      jjj_state.parity = (sh_a.l+sh_b.l+sh_c.l)&1;
	
	      // Using the triangle inequality
	      // to couple the two first
	      // single-particle states
	      for (jjj_state.j_ab = j_ab_min;
		   jjj_state.j_ab <=j_ab_max;
		   jjj_state.j_ab+=2)
		{
		  quantum_number j_abc_min = abs(sh_c.j-jjj_state.j_ab);
		  quantum_number j_abc_max = sh_c.j+jjj_state.j_ab;
		  // Using the triangle inequality
		  // to couple the third
		  // single-particle state
		  // with the total angular momentum
		  // of the two first particle
		  for (jjj_state.j_abc = j_abc_min;
		       jjj_state.j_abc <=j_abc_max;
		       jjj_state.j_abc+=2)
		    {
		      // If the previous guess of the
		      // final size was to small
		      // we resize the array
		      if (jjj_basis->dimension == max_dimension)
			{
			  max_dimension = 2*max_dimension+1;
			  jjj_basis->states = (JJJ_State*)realloc(jjj_basis->states,
								  sizeof(JJJ_State)*max_dimension);
			}
		      jjj_basis->states[jjj_basis->dimension++]=jjj_state;
		    }
		}
	    }
	}
    }
  jjj_basis->states = (JJJ_State*)realloc(jjj_basis->states,
					  sizeof(JJJ_State)*max_dimension);

  // The order in which we have generated the states
  // is not the most optimal order to do the
  // decoupling in, and therefore we here sort the
  // list in to a more optimal order
  qsort_r(jjj_basis->states,jjj_basis->dimension,sizeof(JJJ_State),
	  comp_jjj_states,shells);
  return jjj_basis;
}

/* The following routine 
 * looks for the index of a given state
 * if it do not find the state it returns -1
 * if it do find the state it returns its index
 */
ssize_t find_jjj_state(JJJ_State* states,
		       size_t dim,
		       JJJ_State look_for){

  /* This could perhaps be done more optimally
   * but for the moment I use a linear search
   * to find the state
   */
  size_t i;
  for (i = 0; i<dim; i++){
    if (states[i].a == look_for.a &&
	states[i].b == look_for.b &&
	states[i].c == look_for.c &&
	states[i].j_ab == look_for.j_ab &&
	states[i].j_abc == look_for.j_abc){
      return i;
    }
  }
  return -1;
}

/* The following routine
 * appends a jjj-coupled state
 * to the end of an array of
 * jjj-coupled states
 */
JJJ_State* append_state(JJJ_State *states,
			size_t *dim,
			size_t *dim_max,
			JJJ_State cand){
  // First it is necessary to check
  // if the state do not already exist
  ssize_t i = find_jjj_state(states,*dim,cand);
  if (i<0){
    // If there are no more space
    // to append the state into
    // we have to allocate more
    if (*dim == *dim_max){
      *dim_max = (*dim_max)*2 +1;
      states =
	(JJJ_State*)realloc(states,
			    sizeof(JJJ_State)*(*dim_max));
    }
    states[(*dim)++] = cand;
  }
  return states;
}


/* Given a m-scheme three-particle basis
 * the following function generates 
 * the corresponding jjj-coupled basis
 */
JJJ_Basis* new_jjj_basis_m_scheme(M_Scheme_3p_Basis* ms_basis){
  JJJ_Basis* jjj_basis = (JJJ_Basis*)malloc(sizeof(JJJ_Basis));
  //jjj_basis->e_max = ms_basis->e_max;
  jjj_basis->shells = ms_basis->sp_states->shells;
  jjj_basis->states = NULL;
  jjj_basis->dimension = 0;
  size_t max_dim = 0;
  size_t i;
#define SHELL(a) ms_basis->sp_states->sp_states[a].shell

#define M(a) ms_basis->sp_states->sp_states[a].m
  
#define ORB_ANG(a) ms_basis->sp_states->shells->shells[a].l

#define ANG(a) ms_basis->sp_states->shells->shells[a].j

#define ISO_SPIN(a) ms_basis->sp_states->shells->shells[a].tz

  for (i = 0; i<ms_basis->dimension; i++)
    {
      M_Scheme_3p_State curr = ms_basis->states[i];
      // first sort the curr.abc such that a.shell <= b.shell<= c.shell
      curr = sort_on_shells(curr,
			    ms_basis->sp_states);
      // Now we have the possibilities that
      // a.shell == b.shell == c.shell,
      // a.shell == b.shell < c.shell,
      // a.shell < b.shell == c.shell,
      // a.shell < b.shell < c.shell
      // which we have to consider

      
      // A) a == b == c => a b c && a c b && b c a
#define CASE_A SHELL(curr.a) == SHELL(curr.b) && SHELL(curr.b) == SHELL(curr.c)
      // B) a == b < c =>  a c b && b c a
#define CASE_B SHELL(curr.a) == SHELL(curr.b) && SHELL(curr.b) < SHELL(curr.c)
      // C) a < b == c => a b c && a c b
#define CASE_C SHELL(curr.a) < SHELL(curr.b) && SHELL(curr.b) == SHELL(curr.c)
      // D) a < b < c => a b c
#define CASE_D SHELL(curr.a) < SHELL(curr.b) && SHELL(curr.b) < SHELL(curr.c)
      M_Scheme_3p_State cands[3];
      size_t num;
      if (CASE_A)
	{
	  cands[0] = curr;
	  cands[1] = curr;
	  cands[2] = curr;
	  SWAP(cands[1].b,cands[1].c);
	  CYCLIC_SHIFT_R(cands[2].c,cands[2].a,cands[2].b);
	  num = 3;
	}
      else if (CASE_B)
	{
	  cands[0] = curr;
	  cands[1] = curr;
	  SWAP(cands[0].b,cands[0].c);
	  CYCLIC_SHIFT_L(cands[1].a,cands[1].b,cands[1].c);
	  num = 2;
	}
      else if (CASE_C)
	{
	  cands[0] = curr;
	  cands[1] = curr;
	  SWAP(cands[1].b,cands[1].c);
	  num = 2;
	}
      else // CASE_D
	{
	  cands[0] = curr;
	  num = 1;
	}

      
      size_t k;
      for (k = 0; k<num; k++)
	{
	  // Here we loop over all possible coupling orders and work out the corresponding j-scheme elements
	  JJJ_State current =
	    {SHELL(cands[k].a),
	     SHELL(cands[k].b),
	     SHELL(cands[k].c),
	     0,0,0,0};
	  // computing parity
	  current.parity =
	    (ORB_ANG(current.a) +
	     ORB_ANG(current.b) +
	     ORB_ANG(current.c))&1;

	  // computing total Tz
	  current.tz =
	    ISO_SPIN(current.a)+
	    ISO_SPIN(current.b)+
	    ISO_SPIN(current.c);

	  // computing  m_ab and m_abc

	  quantum_number m_ab =
	    M(cands[k].a) + M(cands[k].b);
	  DEBUG_MESS("m_ab = %d\n",m_ab);
	  
	  quantum_number m_abc =
	    m_ab + M(cands[k].c);
	  DEBUG_MESS("m_abc = %d\n",m_abc);
	  // computing the j_ab loop limits using the triangular conditions

	  quantum_number j_ab_min =
	    MAX(abs(ANG(current.a)-ANG(current.b)),
		abs(m_ab));

	  quantum_number j_ab_max =
	    ANG(current.a)+ANG(current.b);

	  // To take the pauli principle in to account
	  quantum_number j_ab_step =
	    current.a == current.b ? 4 : 2;

	  for (current.j_ab  = j_ab_min;
	       current.j_ab <= j_ab_max;
	       current.j_ab += j_ab_step)
	    {
	      // compute the j_abc loop limits using the triangular conditions
	      quantum_number j_abc_min =
		MAX(abs(current.j_ab-ANG(current.c)),
		    abs(m_abc));

	      quantum_number j_abc_max =
		current.j_ab+ANG(current.c);

	      quantum_number j_abc_step = 2;
	      
	      for (current.j_abc  = j_abc_min;
		   current.j_abc <= j_abc_max;
		   current.j_abc += j_abc_step)
		{
		  // determin if the current state already exists

		  // if not, add it to the basis, increase the allocation size if needed

		  jjj_basis->states =  append_state(jjj_basis->states,
						    &jjj_basis->dimension,
						    &max_dim,
						    current);
		}
	    }
	}
    }
  // As usual the guessed size of the basis, migth not be correct
  // so we need to trim it to the correct size
  jjj_basis->states = (JJJ_State*)realloc(jjj_basis->states,
					  sizeof(JJJ_State)*jjj_basis->dimension);
  //printf("There are %ld, jjj_states\n",jjj_basis->dimension);
  qsort_r(jjj_basis->states,jjj_basis->dimension,sizeof(JJJ_State),
	  comp_jjj_states,jjj_basis->shells);
  return jjj_basis;

}

JJJ_Basis* new_jjj_basis_from_ascii(FILE* config_file,
				    Shells* shells)
{
  size_t n = 0;
  char* row = NULL;
  JJJ_Basis* out = NULL;
  size_t allocated_size = 0;
  // We allow for any number of blank
  // rows in the begining of the file
  do
    {
      if (getline(&row,&n,
		  config_file)<0)
	{
	  if (row != NULL)
	    free(row);
	    
	  return NULL;
	}
    }
  while (strcmp(row,"\n")==0);
  // Verify that it is the right type of block
  if (strcmp(row,
	     "---three-particle-states---\n") == 0)
    {
      out = (JJJ_Basis*)malloc(sizeof(JJJ_Basis));
      out->dimension = 0;
      out->shells = shells;
      out->states = NULL;
    }
  else
    {
      free(row);
      return NULL;
    }
  // skip any aditional comments or blanck rows
  do
    {
      if (getline(&row,&n,
		  config_file)<0)
	{
	  free(row);
	  free(out);
	  return NULL;
	}
    }
  while (strcmp(row,"\n")==0 ||
	 begins_with(row,"---"));

  size_t position;
  do
    {
      JJJ_State jjj_state;
      if (sscanf(row,
		 "%*d %ld %ld %ld %d %d",
		 &jjj_state.a,
		 &jjj_state.b,
		 &jjj_state.c,
		 &jjj_state.j_ab,
		 &jjj_state.j_abc)<5 &&
	  strcmp(row,"\n")!=0)
	{
	  fprintf(stderr,
		  "wrong line format in \"%s\"\n",
		  row);
	  free(row);
	  free(out);
	  exit(1);
	}
      else if (strcmp(row,"\n") != 0)
	{
	  jjj_state.tz =
	    shells->shells[jjj_state.a].tz+
	    shells->shells[jjj_state.b].tz+
	    shells->shells[jjj_state.c].tz;
	  jjj_state.parity =
	    (shells->shells[jjj_state.a].l+
	     shells->shells[jjj_state.b].l+
	     shells->shells[jjj_state.c].l)&1;

	  if (allocated_size == out->dimension)
	    {
	      allocated_size=allocated_size*2+1;
	      out->states =
		(JJJ_State*)realloc(out->states,
				    sizeof(JJJ_State)*
				    allocated_size);
	    }
	  out->states[out->dimension++] = jjj_state;
	}

      position = ftell(config_file);
      if (getline(&row,&n,
		  config_file)<0)
	{
	  break;
	}
    }
  while (!begins_with(row,"---"));
  if (begins_with(row,"---"))
    {
      fseek(config_file,
	    position,
	    SEEK_SET);
    }

  out->states =
    (JJJ_State*)realloc(out->states,
			sizeof(JJJ_State)*
			out->dimension);
  return out;
}


JJJ_Basis* get_block(JJJ_Basis* in_jjj_basis,
		     quantum_number j_abc,
		     quantum_number tz,
		     quantum_number parity)
{
  JJJ_Basis* jjj_basis =
    (JJJ_Basis*)malloc(sizeof(JJJ_Basis));
  jjj_basis->shells = in_jjj_basis->shells;
  ssize_t low_index = -1;
  size_t high_index = 0;
  size_t i;
  for (i = 0; i<in_jjj_basis->dimension; i++)
    {
      if (in_jjj_basis->states[i].j_abc == j_abc &&
	  in_jjj_basis->states[i].tz == tz &&
	  in_jjj_basis->states[i].parity == parity)
	{
	  low_index = i;
	  break;
	}
    }
  if (low_index < 0)
    {
      free(jjj_basis);
      return NULL;
    }
  for (i = low_index+1;
       i<in_jjj_basis->dimension;
       i++)
    {
      if (in_jjj_basis->states[i].j_abc != j_abc ||
	  in_jjj_basis->states[i].tz != tz ||
	  in_jjj_basis->states[i].parity != parity)
	{
	  high_index = i;
	  goto end_of_search; // Used only for a "for if else" construction
	}
    }
  high_index = in_jjj_basis->dimension;
 end_of_search: // Used only for a "for if else" construction
  jjj_basis->dimension = high_index-low_index;
  jjj_basis->states = (JJJ_State*)malloc(sizeof(JJJ_State)*jjj_basis->dimension);

  memcpy(jjj_basis->states,
	 in_jjj_basis->states+low_index,
	 jjj_basis->dimension*sizeof(JJJ_State));
  return jjj_basis;
}

void list_jjj_states(JJJ_Basis* jjj_basis)
{
  size_t i;
  printf("jjj_states:\n");
  printf("(#):\ta\tb\tc\t|\tj_ab\tj_abc\ttz\tparity\n");
  for (i = 0;
       i<jjj_basis->dimension;
       i++)
    {
      printf("(%ld):\t%ld\t%ld\t%ld\t|\t%d\t%d\t%d\t%d\n",
	     i,
	     jjj_basis->states[i].a,
	     jjj_basis->states[i].b,
	     jjj_basis->states[i].c,
	     jjj_basis->states[i].j_ab,
	     jjj_basis->states[i].j_abc,
	     jjj_basis->states[i].tz,
	     jjj_basis->states[i].parity);
    }
}

// compares two state,
// it assumes that they
// are equal if all fields are
// equal down to a permuation
// of the a and b fields
//inline
int fast_cmp(JJJ_State a,
	     JJJ_State b)
{
  if ((a.a == b.a && a.b == b.b)&&
      a.c == b.c && a.j_ab == b.j_ab &&
      a.j_abc == b.j_abc && a.tz == b.tz &&
      a.parity == b.parity)
    return 1;
  else if ((a.b == b.a && a.a == b.b)&&
	   a.c == b.c && a.j_ab == b.j_ab &&
	   a.j_abc == b.j_abc && a.tz == b.tz &&
	   a.parity == b.parity)
    return -1;
  else
    return 0;
}

ssize_t find_trans_jjj_state(JJJ_Basis* looking_in,
			     JJJ_State looking_for,
			     ssize_t *shell_translation)
{
  // translate state if any of states has non existing shells
  // there are no corresponding state
  if (shell_translation[looking_for.a] < 0)
    return -1;
  looking_for.a =
    shell_translation[looking_for.a];
  if (shell_translation[looking_for.b] < 0)
    return -1;
  looking_for.b =
    shell_translation[looking_for.b];
  if (shell_translation[looking_for.c] < 0)
    return -1;
  looking_for.c =
    shell_translation[looking_for.c];  
  // now we know that the state is in the correct format
  // down to a permutation of a and b
  size_t i;
  for (i = 0;
       i<looking_in->dimension;
       i++)
    {
      int phase = 0;
      if ((phase = fast_cmp(looking_for,
			    looking_in->states[i])))
	{
	  return (i+1)*phase;
	}

    }
  return 0;
}

ssize_t* matching_jjj_states(JJJ_Basis *jjj_basis_a,
			     JJJ_Basis *jjj_basis_b,
			     sshell_index* a_to_b_shells)
{
  // we need to find all JJJ_State:s in jjj_basis_a that
  // has a corresponding JJJ_State in jjj_basis_b
  // and if so what index it has
  // for JJJ_State:s in jjj_basis_a that do not
  // have a corresponding JJJ_State in jjj_basis_b
  // we assign -1. Further the shells used for a and
  // b might not be the same therefore we use the
  // translation table a_to_b_shells when comparing

  // To do this we first need an array to store the result in
  // must have the same length as jjj_basis_a->states
  ssize_t* corresponding =
    (ssize_t*)malloc(sizeof(ssize_t)*jjj_basis_a->dimension);

  // now we go through each and every state in jjj_basis_a
  // to find the corresponding state in jjj_basis_b if such exist
  size_t i;
  for (i = 0;
       i < jjj_basis_a->dimension;
       i++)
    {
      // search for corresponding state in jjj_basis_b
      corresponding[i] =
	find_trans_jjj_state(jjj_basis_b, // looking in
			     jjj_basis_a->states[i], // looking for
			     a_to_b_shells); //shell translation
    }
  return corresponding;
}


void free_jjj_basis(JJJ_Basis* jjj_basis)
{
  free(jjj_basis->states);
  free(jjj_basis);
}


// Tests


BEGIN_TEST(trivial_test_new_jjj_basis_m_scheme)
{
  Shells *shells = new_shells(0);
  SP_States *sp_states = new_sp_states(shells);
  M_Scheme_3p_Basis in_basis;
  in_basis.dimension = 1;
  M_Scheme_3p_State m_state[1] = {{0,1,2}};
  in_basis.states = m_state;
  in_basis.sp_states = sp_states;
  in_basis.e_max = 0;

  JJJ_Basis* out_basis = new_jjj_basis_m_scheme(&in_basis);
  list_jjj_states(out_basis);
  ASSERT(out_basis->dimension == 3,
	 TEST_FAILED("out_basis->dimension = %ld, but should equal 3\n",
		     out_basis->dimension),
	 "");
  JJJ_State expected[3] =
    {{0,1,0,0,1,-1,0},
     {0,1,0,2,1,-1,0},
     {0,1,0,2,3,-1,0}};

  size_t i;
  for (i = 0; i<3; i++)
    {
      ASSERT(fast_cmp(expected[i],out_basis->states[i]) == 1,
	     TEST_FAILED("state[%ld] is wrong\n",i),
	     "");
    }
  (void)(expected);
  
  free_jjj_basis(out_basis);
  free_sp_states(sp_states);
  free_shells(shells);
}
END_TEST

BEGIN_TEST(nmax_0_test)
{
  Shells *shells =
    new_shells(0);
  SP_States *sp_states =
    new_sp_states(shells);
  M_Scheme_3p_Basis* m_basis =
    new_m_scheme_3p_basis_no_m_rest(0,sp_states);

  DEBUG_MESS("m-scheme basis:\n");
  list_m_scheme_3p_basis(m_basis);

  JJJ_Basis *j_basis =
    new_jjj_basis_m_scheme(m_basis);

  DEBUG_MESS("j-scheme basis:\n");
  list_jjj_states(j_basis);
  free_jjj_basis(j_basis);
  free_m_scheme_3p_basis(m_basis);
  free_sp_states(sp_states);
  free_shells(shells);
}
END_TEST
