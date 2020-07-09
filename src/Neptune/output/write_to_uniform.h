#ifndef __WRITE_TO_UNIFORM__
#define __WRITE_TO_UNIFORM__

#include <stdlib.h>
#include <stdio.h>
#include "data_block.h"
#include <bases/m_scheme_3p_basis.h>
#include <bases/m_scheme_2p_basis.h>
#include <matrix_transform/matrix_transform.h>

typedef struct _uniform_out_file_
{
	char *directory_name;
	FILE *header_file;
	FILE *basis_file;
	Data_Block* data_blocks;
	block_status_t *block_status;
	size_t current_number_blocks;
	size_t max_number_blocks;
	double *phase;
} Uniform_Out_File;

Uniform_Out_File *create_new_uniform_2p_out_file(const char *directory_name,
						 m_scheme_2p_basis_t ms_basis);

Uniform_Out_File *create_new_uniform_3p_out_file(const char *directory_name,
						 M_Scheme_3p_Basis *ms_basis);

size_t add_uniform_block(Uniform_Out_File *file,
			 Data_Block new_block);

void write_to_uniform_block(Uniform_Out_File *file,
			    size_t block_number,
			    Dens_Matrix elements,
			    size_t *m_indices,
			    size_t *n_indices);

void close_uniform_out_file(Uniform_Out_File* file);
#endif
