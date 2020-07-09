#ifndef __BLOCK_TRANSFORM__
#define __BLOCK_TRANSFORM__

#include <bases/m_scheme_2p_basis.h>
#include <bases/m_scheme_3p_basis.h>
#include <matrix_transform/matrix_transform.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <input/read_3nf_file.h>
#include <input/read_2nf_antoine_format.h>

Dens_Matrix* compute_jjj_block(M_Scheme_3p_Basis* bra_basis,
			       M_Scheme_3p_Basis* ket_basis,
			       Data_File* datafile,
			       Clebsch_Gordan_Data* cgd);

Dens_Matrix* compute_jt_block(m_scheme_2p_basis_t bra_basis,
			      m_scheme_2p_basis_t ket_basis,
			      quantum_number Tz,
			      quantum_number M,
			      quantum_number J_max,
			      antoine_2nf_file_t data_file,
			      Clebsch_Gordan_Data *cgd);
#endif
