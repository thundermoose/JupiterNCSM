#ifndef __OUT_FILE__
#define __OUT_FILE__

#include <stdlib.h>
#include <output/data_block.h>
#include <bases/m_scheme_3p_basis.h>
#include <matrix_transform/matrix_transform.h>

struct _out_file_;

typedef struct _out_file_ *out_file_t;

typedef enum 
{
	unified_2p_file,
	unified_3p_file,
	hdf5_file,
	unknown_format
} file_type_t;

file_type_t file_type_from_string(const char *file_type);

out_file_t create_new_out_file(const char *file_name,
				void *ms_basis,
				file_type_t file_type);

size_t add_block(out_file_t file,
		Data_Block data_blocks);

void write_to_block(out_file_t file,
		size_t block_number,
		Dens_Matrix elements,
		size_t *bra_indices,
		size_t *ket_indices);

void close_out_file(out_file_t file);

#endif
