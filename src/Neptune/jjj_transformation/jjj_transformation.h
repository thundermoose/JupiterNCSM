#ifndef __JJJ_TRANSFORMATION__
#define __JJJ_TRANSFORMATION__

#include <matrix_transform/matrix_transform.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <bases/m_scheme_3p_basis.h>
#include <bases/jjj_coupled_3p.h>

/* This function generates a transformation matrix between 
 * the input M_Scheme_3p_Basis and the JJJ_Basis.
 */
Sparse_Matrix* new_jjj_transformation(M_Scheme_3p_Basis* ms3b,
				      JJJ_Basis* jjj_basis,
				      Clebsch_Gordan_Data* cgd);


#endif
