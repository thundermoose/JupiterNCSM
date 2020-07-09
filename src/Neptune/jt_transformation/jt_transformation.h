#ifndef __JT_TRANSFORMATION__
#define __JT_TRANSFORMATION__

#include <matrix_transform/matrix_transform.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <bases/m_scheme_2p_basis.h>
#include <bases/jt_basis.h>
#include <bases/m_scheme_3p_basis.h>
#include <bases/jt_basis_3p.h>

/* This function generates a transformation matrix between 
 * the input M_Scheme_3p_Basis and the JJJ_Basis.
 */
Sparse_Matrix* new_jt_transformation_3p(M_Scheme_3p_Basis* ms3b,
					JT_Basis* jt_basis,
					Clebsch_Gordan_Data* cgd);

Sparse_Matrix* new_jt_transformation(m_scheme_2p_basis_t m_scheme_2p_basis,
				     jt_basis_t jt_basis,
				     Clebsch_Gordan_Data *clebsch_gordan_data);
#endif
