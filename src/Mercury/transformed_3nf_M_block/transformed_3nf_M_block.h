#ifndef __TRANSFORMED_3NF_M_BLOCK__
#define __TRANSFORMED_3NF_M_BLOCK__

#include <bases/m_scheme_3p_basis.h>
#include <matrix_transform/matrix_transform.h>

typedef struct
{
	M_Scheme_3p_Basis *ket_basis;
	M_Scheme_3p_Basis *bra_basis;
	Dens_Matrix *matrix;
} transformed_3nf_M_block_t;

#endif
