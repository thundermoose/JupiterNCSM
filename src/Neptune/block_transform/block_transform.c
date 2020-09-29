#include <stdlib.h>
#include <math.h>
#include <block_transform/block_transform.h>
#include <bases/jjj_coupled_3p.h>
#include <bases/jt_basis_3p.h>
#include <utils/assertion.h>
#include <utils/debug_messages.h>
#include <utils/helpful_macros.h>
#include <jjj_transformation/jjj_transformation.h>
#include <jt_transformation/jt_transformation.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>

#ifdef DEBUG
static inline
void log_matrix(Dens_Matrix *matrix)
{
	log_entry("matrix (%p):",matrix);
	for (size_t i = 0; i<matrix->m; i++)
	{
		for (size_t j = 0; j<matrix->n; j++)
		{
			log_entry("(%lu %lu) = %g",
				  i,j,
				  matrix->elements[matrix->n*i+j]);
		}
	}
}
#endif

// The following four functions
// compute the extremal values
// of Jabc and Tz that exists
// within a given basis

quantum_number find_min_jabc(M_Scheme_3p_Basis* basis)
{
	// Since we know that J>=0
	// we use -1 to represent
	// that this min_jabc is not
	// set yet
	quantum_number min_jabc= -1;
	size_t i;
	for (i = 0; i<basis->dimension; i++)
	{
		// For each state s we use the triangular
		// condition for angular momentum coupling
		// twice to compute the smallest possible J_abc
		// that state s can correspond to
		M_Scheme_3p_State s = basis->states[i];
		SP_State a = basis->sp_states->sp_states[s.a];
		SP_State b = basis->sp_states->sp_states[s.b];
		SP_State c = basis->sp_states->sp_states[s.c];
		quantum_number candidate_jabc =
			max(abs(abs(basis->sp_states->shells->shells[a.shell].j-
							basis->sp_states->shells->shells[b.shell].j)-
						basis->sp_states->shells->shells[c.shell].j),
					abs(a.m+b.m+c.m));

		// if the computed J is smaller than
		// the currently smallest (or if min_jabc == -1)
		// we know that min_jabc is not the smallest
		if (candidate_jabc<min_jabc ||
				min_jabc == -1)
		{
			min_jabc=candidate_jabc;
		}
	}
	return min_jabc;
}

quantum_number find_max_jabc(M_Scheme_3p_Basis* basis)
{
	quantum_number max_jabc= -1;
	size_t i;
	for (i = 0; i<basis->dimension; i++)
	{
		// Exactly the same principle as for find_min_jabc
		// except that we use the triangular condition
		// to compute the maximum jabc and that < is
		// replaced by >
		M_Scheme_3p_State s = basis->states[i];
		SP_State a = basis->sp_states->sp_states[s.a];
		SP_State b = basis->sp_states->sp_states[s.b];
		SP_State c = basis->sp_states->sp_states[s.c];
		quantum_number candidate_jabc =
			basis->sp_states->shells->shells[a.shell].j+
			basis->sp_states->shells->shells[b.shell].j+
			basis->sp_states->shells->shells[c.shell].j;
		if (candidate_jabc>max_jabc || max_jabc == -1)
		{
			max_jabc=candidate_jabc;
		}
	}
	return max_jabc;
}


quantum_number find_min_tz(M_Scheme_3p_Basis* basis)
{
	quantum_number min_tz= -7;
	size_t i;
	for (i = 0; i<basis->dimension; i++)
	{
		// This follows a similar principle as
		// for find_min_Jabc and find_max_Jabc
		// except that we compute the states Tz
		M_Scheme_3p_State s = basis->states[i];
		SP_State a = basis->sp_states->sp_states[s.a];
		SP_State b = basis->sp_states->sp_states[s.b];
		SP_State c = basis->sp_states->sp_states[s.c];
		quantum_number candidate_tz =
			basis->sp_states->shells->shells[a.shell].tz+
			basis->sp_states->shells->shells[b.shell].tz+
			basis->sp_states->shells->shells[c.shell].tz;
		if (candidate_tz<min_tz || min_tz == -7)
		{
			min_tz=candidate_tz;
		}
	}
	return min_tz;
}

quantum_number find_max_tz(M_Scheme_3p_Basis* basis)
{
	quantum_number max_tz= -1;
	size_t i;
	for (i = 0; i<basis->dimension; i++)
	{
		// See find_min_tz for explanation
		M_Scheme_3p_State s = basis->states[i];
		SP_State a = basis->sp_states->sp_states[s.a];
		SP_State b = basis->sp_states->sp_states[s.b];
		SP_State c = basis->sp_states->sp_states[s.c];
		quantum_number candidate_tz =
			basis->sp_states->shells->shells[a.shell].tz+
			basis->sp_states->shells->shells[b.shell].tz+
			basis->sp_states->shells->shells[c.shell].tz;
		if (candidate_tz>max_tz || max_tz == -1)
		{
			max_tz=candidate_tz;
		}
	}
	return max_tz;
}


Dens_Matrix* compute_jjj_block(M_Scheme_3p_Basis* bra_basis,
		M_Scheme_3p_Basis* ket_basis,
		Data_File* datafile,
		Clebsch_Gordan_Data* cgd)
{


	log_entry("\n\n========================New Block=============================\n\n");
	log_entry("Test\n");

	// This matrix is the one that will
	// be returned, we will build it up
	// through out the method, by adding
	// the result from each Jabc,Tz and
	// parity contribution.
	Dens_Matrix* acc =
		new_zero_matrix(bra_basis->dimension,
				ket_basis->dimension);


	// We here computes lists of the fully
	// j-coupled bases, we need one for
	// the bra and one for the ket side
	// since we might be computing an
	// of diagonal matrix-block
	JJJ_Basis *bra_jjj_basis =
		new_jjj_basis_m_scheme(bra_basis);

	ASSERT(is_subset_of_basis(bra_jjj_basis,
				datafile),
			exit(1),
			"bra_jjj_basis is not a subset of datafile->configurations\n");

	// Debug code
	log_entry("bra_m_scheme_basis:\n");
	DEBUG_CALL(list_m_scheme_3p_basis(bra_basis));
	log_entry("bra_jjj_basis:\n");
	DEBUG_CALL(list_jjj_states(bra_jjj_basis));

	JJJ_Basis *ket_jjj_basis =
		new_jjj_basis_m_scheme(ket_basis);

	ASSERT(is_subset_of_basis(ket_jjj_basis,
				datafile),
			exit(1),
			"ket_jjj_basis is not a subset of datafile->configurations\n");

	// Debug code
	log_entry("ket_m_scheme_basis:\n");
	DEBUG_CALL(list_m_scheme_3p_basis(ket_basis));
	log_entry("ket_jjj_basis:\n");
	DEBUG_CALL(list_jjj_states(ket_jjj_basis));

	// Determining the limits for the loops
	quantum_number min_bra_jabc = find_min_jabc(bra_basis);
	quantum_number min_ket_jabc = find_min_jabc(ket_basis);
	quantum_number max_bra_jabc = find_max_jabc(bra_basis);
	quantum_number max_ket_jabc = find_max_jabc(ket_basis);
	quantum_number min_jabc = max(min_bra_jabc,
			min_ket_jabc);
	quantum_number max_jabc = min(max_bra_jabc,
			max_ket_jabc);

	quantum_number min_bra_tz = find_min_tz(bra_basis);
	quantum_number min_ket_tz = find_min_tz(ket_basis);
	quantum_number max_bra_tz = find_max_tz(bra_basis);
	quantum_number max_ket_tz = find_max_tz(ket_basis);
	quantum_number min_tz = max(min_bra_tz,
			min_ket_tz);
	quantum_number max_tz = min(max_bra_tz,
			max_ket_tz);
	// respecting conserved observables
	// if we do not have a over lap of the J (or Tz)
	// intervals that are allowed by the triangular
	// conditions we can not have any states that
	// have equal J (or Tz)
	// this can happen for example
	// if min_bra_jabc > max_ket_jabc
	if (min_jabc>max_jabc ||
			min_tz>max_tz)
	{
		return acc; 
	}

	quantum_number jabc,tz,parity;
	// We devide the transformation in to blocks
	// since the full matrix is block diagonal in
	// theses blocks. 
	for (jabc  = min_jabc;
			jabc <= max_jabc;
			jabc += 2) // dJ = 1 and since we multiply
		// it with 2 jabc should increase
		// with 2
	{
		for (tz=min_tz; tz<=max_tz; tz+=2)
		{
			for (parity = 0; parity <= 1; parity++)
			{

				// We compute the portions
				// of the bra and ket basis
				// that corresponds to the current
				// block, if any
				// if no such portion are found
				// this block is skiped
				JJJ_Basis* cut_bra_basis =
					get_block(bra_jjj_basis,
							jabc,tz,parity);
				if (cut_bra_basis == NULL)
				{
					continue;
				}

				log_entry("cut_bra_basis:\n");
				DEBUG_CALL(list_jjj_states(cut_bra_basis));

				JJJ_Basis* cut_ket_basis =
					ket_jjj_basis == bra_jjj_basis ?
					cut_bra_basis :
					get_block(ket_jjj_basis,
							jabc,tz,parity);
				if (cut_ket_basis == NULL)
				{
					free_jjj_basis(cut_bra_basis);
					continue;
				}

				log_entry("cut_ket_basis:\n");
				DEBUG_CALL(list_jjj_states(cut_ket_basis));


				log_entry("\n\n\nObtaining elements for: %d %d %d\n",
						jabc,tz,parity);

				// Now we know what matrix elements
				// we need to read from the datafile
				Dens_Matrix* jjj_matrix =
					get_matrix(datafile,
							cut_bra_basis,
							cut_ket_basis);
				// We do also know enough to compute
				// the transformations from j-scheme
				// to m-scheme
				Sparse_Matrix* bra_transform =
					new_jjj_transformation(bra_basis,
							cut_bra_basis,
							cgd);
				/*
				   ASSERT(is_unitary(bra_transform),
				   print_sparse_matrix(*bra_transform);
				   exit(1),
				   "bra_transform for %d %d %d is not unitary\n",
				   jabc,tz,parity);*/
				Sparse_Matrix* ket_transform =
					new_jjj_transformation(ket_basis,
							cut_ket_basis,
							cgd);
				/*
				   ASSERT(is_unitary(ket_transform),
				   print_sparse_matrix(*ket_transform);
				   exit(1),
				   "ket_transform for %d %d %d is not unitary\n",
				   jabc,tz,parity);*/
				// Now we have everything to apply
				// the transformations
				Dens_Matrix* m_matrix =
					transform_matrix(bra_transform,
							jjj_matrix,
							ket_transform);
				DEBUG_CALL({
						char title[256];
						sprintf(title,
								"JJJ_Block: %d %d %d",
								jabc,tz,parity);
						create_matrix_plot(*jjj_matrix,
								title);
						});	 

#pragma omp critical(to_stdout)
				{
					DEBUG_CALL({
							printf("Transforms\n");
							print_sparse_matrix(*bra_transform);
							printf("X\n");
							print_matrix(*jjj_matrix);
							printf("X\n");
							print_sparse_matrix(*ket_transform);
							printf("|\nV\n");
							print_matrix(*m_matrix);
							});
				}


				// Now we can add the transformed
				// matrix (the one in m-scheme)
				// to the output matrix
				accumulate_matrix(*acc,*m_matrix);

				// We can now free everything we created
				// for this block
				free_sparse_matrix(bra_transform);
				free_sparse_matrix(ket_transform);
				free_dens_matrix(jjj_matrix);
				free_jjj_basis(cut_bra_basis);
				if (ket_jjj_basis != bra_jjj_basis)
				{
					free_jjj_basis(cut_ket_basis);
				}
				free_dens_matrix(m_matrix);
			}
		}
	}
	// These bases are not usefull anymore
	free_jjj_basis(bra_jjj_basis);
	free_jjj_basis(ket_jjj_basis);
	return acc;
}

Dens_Matrix* compute_jt_block(m_scheme_2p_basis_t bra_basis,
			      m_scheme_2p_basis_t ket_basis,
			      quantum_number Tz,
			      quantum_number M,
			      quantum_number J_max,
			      antoine_2nf_file_t data_file,
			      Clebsch_Gordan_Data* cgd)
{
	log_entry("Tz = %d, M = %d, J_max = %d", Tz,M,J_max);
	jt_basis_t jt_basis = get_jt_basis(data_file);

	Dens_Matrix* accumulator =
		new_zero_matrix(get_m_scheme_2p_dimension(bra_basis),
				get_m_scheme_2p_dimension(ket_basis));
	for (quantum_number J = abs(M)/2; J<=J_max/2; J++)
		for (quantum_number T = abs(Tz)/2; T<=1; T++)
		{
			jt_basis_t jt_block_basis =
				get_jt_block_basis(jt_basis,
						   J,T);
			if (jt_block_basis == NULL)
				continue;
			log_entry("J = %d, T = %d",J,T);
			Dens_Matrix* jt_potential =
				get_antoine_matrix(data_file,
						   jt_block_basis,
						   jt_block_basis,
						   Tz);
			log_entry("jt_potential:");
#ifdef DEBUG
			log_matrix(jt_potential);
#endif
			Sparse_Matrix *bra_transform = 
				new_jt_transformation(bra_basis,
						      jt_block_basis,
						      cgd);
			Sparse_Matrix *ket_transform = 
				new_jt_transformation(ket_basis,
						      jt_block_basis,
						      cgd);
			Dens_Matrix* transformed_matrix =
				transform_matrix(bra_transform,
						 jt_potential,
						 ket_transform);
			log_entry("transformed_matrix:");
#ifdef DEBUG
			log_matrix(transformed_matrix);
#endif
			accumulate_matrix(*accumulator,
					  *transformed_matrix);
			free_sparse_matrix(bra_transform);
			free_sparse_matrix(ket_transform);
			free_dens_matrix(jt_potential);
			free_dens_matrix(transformed_matrix);
			free_jt_basis(jt_block_basis);
		}
	return accumulator;
}

Dens_Matrix* compute_3jt_block(M_Scheme_3p_Basis* bra_basis,
		M_Scheme_3p_Basis* ket_basis,
		Data_File* datafile,
		Clebsch_Gordan_Data* cgd)
{
	// Same principle as compute_jjj_block
	// but we do it for JT

	// creating a accumulator matrix
	Dens_Matrix* acc =
		new_zero_matrix(bra_basis->dimension,
				ket_basis->dimension);

	// Set up the full JT-bases
	// corresponding to the
	// input m-scheme bases
	JT_Basis *bra_jt_basis =
		new_jt_basis_from_m_scheme(bra_basis);
	JT_Basis *ket_jt_basis =
		new_jt_basis_from_m_scheme(ket_basis);

	// Computing the limits for the loops
	// now also including T but not tz
	quantum_number min_bra_jabc = find_min_jabc(bra_basis);
	quantum_number min_ket_jabc = find_min_jabc(ket_basis);
	quantum_number max_bra_jabc = find_max_jabc(bra_basis);
	quantum_number max_ket_jabc = find_max_jabc(ket_basis);
	quantum_number min_jabc = max(min_bra_jabc,
			min_ket_jabc);
	quantum_number max_jabc = min(max_bra_jabc,
			max_ket_jabc);
	if (min_jabc>max_jabc)
	{
		// jabc is a conserved
		// quantity
		return acc;
	}

	// All particles have t = 1/2,
	// therefore T for 3 particles
	// can either be 1/2 or 3/2
	// and as usual, T and J are
	// internaly multiplied with 2
	// to fit in integers
	quantum_number min_tabc = 1;
	quantum_number max_tabc = 3;

	quantum_number jabc,tabc;

	// looping of the allowed JT-blocks
	for (jabc  = min_jabc;
			jabc <= max_jabc;
			jabc += 2) 
	{
		for (tabc  = min_tabc;
				tabc <= max_tabc;
				tabc += 2)
		{
			// compute the portions of
			// the JT-bases that
			// correspond to the current
			// JT-block if no such
			// portion exists
			// go to the next block

			JT_Basis* cut_bra_jt_basis =
				get_jt_block(bra_jt_basis,
						jabc,
						tabc);
			JT_Basis* cut_ket_jt_basis =
				get_jt_block(ket_jt_basis,
						jabc,
						tabc);

			// now we can compute
			// the transformations

			/* place code here */

			// and load the corresponding
			// JT-matrix elements

			Dens_Matrix* jjj_matrix =
				get_matrix(datafile,
						cut_bra_jt_basis,
						cut_ket_jt_basis);
			intentionaly_unused(jjj_matrix);
			// Everything need to
			// do the transformation
			// exists now

			/* place code here */

			// Add the result to the
			// accumulator matrix

			/* place code here */

			// Everything JT-block specific
			// can now be freed

			free_jt_basis_3p(cut_ket_jt_basis);
			free_jt_basis_3p(cut_bra_jt_basis);

		}
	}

	// Every help structure created
	// can now be freed
	free_jt_basis_3p(ket_jt_basis);
	free_jt_basis_3p(bra_jt_basis);
	return acc;
}
