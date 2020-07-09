#ifndef __TRANSFORM_SCHEDULLER_3P__
#define __TRANSFORM_SCHEDULLER_3P__
#include <input/read_2nf_antoine_format.h>
#include <bases/m_scheme_2p_basis.h>
#include <output/out_file.h>

void transform_2p_data(antoine_2nf_file_t data,
		       out_file_t outfile,
		       m_scheme_2p_basis_t basis);

#endif
