#ifndef __TRANSFORM_SCHEDULLER_3P__
#define __TRANSFORM_SCHEDULLER_3P__
#include <input/read_3nf_file.h>
#include <bases/m_scheme_3p_basis.h>
#include <output/out_file.h>

void transform_3p_data(Data_File* data,
		       out_file_t outfile,
		       M_Scheme_3p_Basis* basis);

#endif
