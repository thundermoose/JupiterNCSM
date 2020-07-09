#include <stdlib.h>
#include <omp.h>
#include <transform_scheduller/transform_scheduller_3p.h>
#include <block_transform/block_transform.h>
#include <utils/range.h>

	static
void task_handler(Data_Block current_block,
		  size_t block_number,
		  Clebsch_Gordan_Data *cgd,
		  out_file_t out_file,
		  Data_File* data_file,
		  M_Scheme_3p_Basis* basis)
{
	size_t bra_offset = 0;
	M_Scheme_3p_Basis* bra_basis
		= generate_block(basis,
				 current_block.Tz,
				 current_block.M,
				 current_block.E1,
				 &bra_offset);
	if (bra_basis == NULL)
	{
		return;
	}
	size_t ket_offset = bra_offset;
	M_Scheme_3p_Basis* ket_basis =
		current_block.E1 == current_block.E2 ?
		bra_basis :
		generate_block(basis,
			       current_block.Tz,
			       current_block.M,
			       current_block.E2,
			       &ket_offset);
	if (ket_basis == NULL)
	{
		free_m_scheme_3p_basis(bra_basis);
		return;
	}
	Dens_Matrix* matrix_block =
		compute_jjj_block(bra_basis,
				  ket_basis,
				  data_file,
				  cgd);

	size_t *bra_indices =
		range(bra_offset,bra_offset+bra_basis->dimension);
	size_t *ket_indices =
		range(ket_offset,ket_offset+ket_basis->dimension);
	write_to_block(out_file,
		       block_number,
		       *matrix_block,
		       bra_indices,
		       ket_indices);
	free(bra_indices);
	free(ket_indices);


	free_dens_matrix(matrix_block);
	if (current_block.E1 != current_block.E2)
		free_m_scheme_3p_basis(ket_basis);
	free_m_scheme_3p_basis(bra_basis);
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
	else if (data_block->Tz<3)
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

void transform_3p_data(Data_File* data_file,
		       out_file_t out_file,
		       M_Scheme_3p_Basis* basis)
{
	quantum_number e_max1 = basis->sp_states->shells->e_max;
	quantum_number e_max2 = basis->e_max;
	quantum_number J_max = 3*(e_max1*2+1);
	Data_Block current_block = 
	{
		.Tz = -3,
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
			while(has_next_block(&current_block,
					     e_max2,
					     J_max))
			{
				size_t block_number =
					add_block(out_file,
						  current_block);

#pragma omp task firstprivate(current_block,cgd,block_number)
				{
					task_handler(current_block,
						     block_number,
						     cgd,
						     out_file,
						     data_file,
						     basis);
				}
			}
		}
		free_clebsch_gordan(cgd);
	}
}
