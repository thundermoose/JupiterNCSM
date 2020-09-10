#include <stdlib.h>
#include <omp.h>
#include <transform_scheduller/transform_scheduller_2p.h>
#include <block_transform/block_transform.h>
#include <utils/debug_messages.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>

#ifndef NDEBUG
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

	static
void task_handler(Data_Block current_block,
		  quantum_number J_max,
		  size_t block_number,
		  Clebsch_Gordan_Data *cgd,
		  out_file_t out_file,
		  antoine_2nf_file_t data_file,
		  m_scheme_2p_basis_t basis)
{
	m_scheme_2p_basis_t bra_basis
		= generate_2p_block(basis,
				    current_block.Tz,
				    current_block.M,
				    current_block.E1);
	if (bra_basis == NULL)
	{
		return;
	}
	m_scheme_2p_basis_t ket_basis =
		current_block.E1 == current_block.E2 ?
		bra_basis :
		generate_2p_block(basis,
				  current_block.Tz,
				  current_block.M,
				  current_block.E2);
	if (ket_basis == NULL)
	{
		free_m_scheme_2p_basis(bra_basis);
		return;
	}
	log_entry("block: %d %d %d %d",
		   current_block.Tz,
		   current_block.M,
		   current_block.E1,
		   current_block.E2);
	log_entry("bra_dimension: %lu",
		   get_m_scheme_2p_dimension(bra_basis));
	log_entry("ket_dimension: %lu",
		   get_m_scheme_2p_dimension(ket_basis));
	log_entry("bra_basis:");
	log_m_scheme_2p_basis(bra_basis);
	log_entry("ket_basis:");
	log_m_scheme_2p_basis(ket_basis);
	Dens_Matrix* matrix_block =
		compute_jt_block(bra_basis,
				 ket_basis,
				 current_block.Tz,
				 current_block.M,
				 J_max,
				 data_file,
				 cgd);
#ifndef NDEBUG
	log_matrix(matrix_block);
#endif
	size_t *bra_indices = m_scheme_2p_corresponding_indices(bra_basis);
	size_t *ket_indices = m_scheme_2p_corresponding_indices(ket_basis);
	write_to_block(out_file,
		       block_number,
		       *matrix_block,
		       bra_indices,
		       ket_indices);


	free_dens_matrix(matrix_block);
	if (current_block.E1 != current_block.E2)
		free_m_scheme_2p_basis(ket_basis);
	free_m_scheme_2p_basis(bra_basis);
}

	static
int has_next_block(Data_Block *data_block,
		   quantum_number e_max,
		   quantum_number J_max)
{
	if (data_block->E2 < e_max-(data_block->E1 & 1))
		data_block->E2+=2;
	else if (data_block->E1 < e_max)
	{
		data_block->E1++;
		data_block->E2 = data_block->E1;
	}
	else if (data_block->M<J_max)
	{
		data_block->M+=2;
		data_block->E1 = 0;
		data_block->E2 = 0;
	}
	else if (data_block->Tz<2)
	{
		data_block->Tz+=2;
		data_block->M = -J_max;
		data_block->E1 = 0;
		data_block->E2 = 0;
	}
	else
	{
		return 0;
	}
	return 1;
}

void transform_2p_data(antoine_2nf_file_t data_file,
		       out_file_t out_file,
		       m_scheme_2p_basis_t basis)
{
	quantum_number e_max1 = get_m_scheme_2p_e_max1(basis);
	quantum_number e_max2 = get_m_scheme_2p_e_max2(basis);
	quantum_number J_max = 2*(e_max1*2+1);
	Data_Block current_block = 
	{
		.Tz = -2,
		.M = -J_max,
		.E1 = 0,
		.E2 = 0
	};
	Clebsch_Gordan_Data* cgd;
#pragma omp parallel
	{
		cgd = initiate_clebsch_gordan(J_max);
		#pragma omp single
		{
			do
			{
				log_entry("current_block: %d %d %d %d",
					   current_block.Tz,
					   current_block.M,
					   current_block.E1,
					   current_block.E2);
				size_t block_number =
					add_block(out_file,
						  current_block);

				#pragma omp task firstprivate(current_block,cgd,block_number)
				{
					task_handler(current_block,
						     J_max,
						     block_number,
						     cgd,
						     out_file,
						     data_file,
						     basis);
				}
			}
			while(has_next_block(&current_block,
					     e_max2,
					     J_max));
		}
		free_clebsch_gordan(cgd);
	}
}
