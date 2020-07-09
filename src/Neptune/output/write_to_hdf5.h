#ifndef __WRITE_TO_HDF5__
#define __WRITE_TO_HDF5__
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <output/data_block.h>
#include <bases/m_scheme_3p_basis.h>
#include <matrix_transform/matrix_transform.h>



typedef struct _hdf5_out_file_
{
  hid_t file;
  Data_Block* data_blocks;
  size_t current_number_blocks;
  size_t max_number_blocks;
} HDF5_Out_File;



HDF5_Out_File* create_new_hdf5_out_file(const char* file_name,
			      M_Scheme_3p_Basis *ms_basis);

size_t add_hdf5_block(HDF5_Out_File* file,
		 quantum_number Tz,
		 quantum_number M,
		 quantum_number E1,
		 quantum_number E2);



void write_to_hdf5_block(HDF5_Out_File* file,
		    size_t block_number,
		    Dens_Matrix elements,
		    size_t *m_indices,
		    size_t *n_indices);

void close_hdf5_out_file(HDF5_Out_File* file);

#endif
