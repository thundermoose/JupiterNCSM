#include <output/out_file.h>
#include <output/write_to_uniform.h>
#include <output/write_to_hdf5.h>
#include <debug_mode/debug_mode.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

struct _out_file_
{
	void *out_file_handler;
	file_type_t file_type;
};

file_type_t file_type_from_string(const char *file_type)
{
	if (strcmp("hdf5",file_type) == 0)
		return hdf5_file;
	else if (strcmp("unified-2p",file_type) == 0)
		return unified_2p_file;
	else if (strcmp("unified-3p",file_type) == 0)
		return unified_3p_file;
	else
		return unknown_format;
}

out_file_t create_new_out_file(const char *file_name,
				void *ms_basis,
				file_type_t file_type)
{
	out_file_t out_file = (out_file_t)malloc(sizeof(struct _out_file_));
	out_file->file_type = file_type;
	switch (file_type)
	{
		case unified_2p_file:
			out_file->out_file_handler =
				create_new_uniform_2p_out_file(file_name,
						(m_scheme_2p_basis_t)ms_basis);
			break;
		case unified_3p_file:
			out_file->out_file_handler =
				create_new_uniform_3p_out_file(file_name,
						(M_Scheme_3p_Basis*)ms_basis);
			break;
		case hdf5_file:
			out_file->out_file_handler = 
				create_new_hdf5_out_file(file_name,
						(M_Scheme_3p_Basis*)ms_basis);
			break;
		default:
			fprintf(stderr,"Unknown file format\n");
			exit(EXIT_FAILURE);
	}
	return out_file;
}

size_t add_block(out_file_t file,
		Data_Block data_blocks)
{
	switch (file->file_type)
	{
		case unified_2p_file:
		case unified_3p_file:
			return add_uniform_block(
					(Uniform_Out_File*)
					file->out_file_handler,
					data_blocks);
			break;
		case hdf5_file:
			return add_hdf5_block(
					(HDF5_Out_File*)
					file->out_file_handler,
					data_blocks.Tz,
					data_blocks.M,
					data_blocks.E1,
					data_blocks.E2);
			break;
		default:
			fprintf(stderr,"Unknown file format\n");
			exit(EXIT_FAILURE);
	}
}

void write_to_block(out_file_t file,
		size_t block_number,
		Dens_Matrix elements,
		size_t *bra_indices,
		size_t *ket_indices)
{
	switch (file->file_type)
	{
		case unified_2p_file:
		case unified_3p_file:
			write_to_uniform_block(
					(Uniform_Out_File*)
					file->out_file_handler,
					block_number,
					elements,
					bra_indices,
					ket_indices);
			break;
		case hdf5_file:
			write_to_hdf5_block(
					(HDF5_Out_File*)
					file->out_file_handler,
					block_number,
					elements,
					bra_indices,
					ket_indices);
			break;
		default:
			fprintf(stderr,"Unknown file format\n");
			exit(EXIT_FAILURE);
	}
}

void close_out_file(out_file_t file)
{
	switch (file->file_type)
	{
		case unified_2p_file:
		case unified_3p_file:
			close_uniform_out_file(
					(Uniform_Out_File*)
					file->out_file_handler);
			break;
		case hdf5_file:
			close_hdf5_out_file(
					(HDF5_Out_File*)
					file->out_file_handler);
			break;
		default:
			fprintf(stderr,"Unknown file format\n");
			exit(EXIT_FAILURE);
	}
	free(file);
}
