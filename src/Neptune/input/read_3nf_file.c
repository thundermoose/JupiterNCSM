#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <input/read_3nf_file.h>
#include <string_tools/string_tools.h>
#include <input/read_3nf_hdf5_file.h>
#include <input/read_3nf_ascii_format.h>
#include <debug_mode/debug_mode.h>
#include <log/log.h>

Data_File* open_data_file(const char* file_name)
{
    // Tries to open the file as a hdf5 file
    HDF5_Data *hdf5_data =
	open_hdf5_data(file_name);
    if (hdf5_data)
	{
	    Data_File* opened_data =
		(Data_File*)malloc(sizeof(Data_File));
	    opened_data->signature = HDF5;
	    opened_data->data_pointer =
		(void*)hdf5_data;
	    opened_data->e_max = 0; // change this to something more appropriate
	    return opened_data;
	}
  
    // Tries to open the file as our ascii format
    ASCII_Data* ascii_data =
	open_ascii_data(file_name);
    if (ascii_data)
	{
	    Data_File* opened_data =
		(Data_File*)malloc(sizeof(Data_File));
	    opened_data->signature = ASCII;
	    opened_data->data_pointer =
		(void*)ascii_data;
	    opened_data->e_max = 0;
	    return opened_data;
	}

    // If we come this far we could not open the data files
    fprintf(stderr,"Could not open the file, not a recognized format\n");
    return NULL;
}

void set_max_loaded_memory(Data_File* data_file,size_t max_loaded_memory)
{
	if (data_file->signature != HDF5)
		return;
	set_hdf5_max_loaded_memory((HDF5_Data*)(data_file->data_pointer),
				   max_loaded_memory);
}

weight_t identify_weight(const char *weight)
{
#define match(W)				\
	if (strcmp(weight,#W)==0)		\
		return W
	match(CE);
	match(CD);
	match(C1);
	match(C3);
	match(C4);
	return U;
}

void set_weight(Data_File* file,
		weight_t weight,
		double value)
{
	if (file->signature != HDF5)
		return;
	HDF5_Data *data_file = file->data_pointer;
        data_file->weights[weight] = value;
}

Dens_Matrix* get_matrix(Data_File* data_file,
		        void* m_basis,
		        void* n_basis)
{
    switch (data_file->signature)
	{
	case HDF5:
	    return get_matrix_hdf5((HDF5_Data*)data_file->data_pointer,
				   (JJJ_Basis*)m_basis,
				   (JJJ_Basis*)n_basis);
	    break;
	case ASCII:
	    return get_matrix_ascii((ASCII_Data*)data_file->data_pointer,
				    m_basis,
				    n_basis);
	    break;
	default:
	    fprintf(stderr,"Could not get matrix from file, unknown format\n");
	}
    return NULL;
}


int is_subset_of_basis(JJJ_Basis* j_scheme,
		       Data_File* data_file)
{
    switch (data_file->signature)
	{
	case HDF5:
	    return is_subset_of_basis_hdf5(j_scheme,
					   data_file->data_pointer);
	    break;
	case ASCII:
	    log_entry("subset_of_basis_ascii is not yet implemented\n");
	    return 1;
	    break;
	default:
	    return 0;
	}

}

void free_data_file(Data_File* data_file)
{
    switch (data_file->signature)
	{
	case HDF5:
	    free_hdf5_data((HDF5_Data*)data_file->data_pointer);
	    break;
	case ASCII:
	    free_ascii_data((ASCII_Data*)data_file->data_pointer);
	    break;
	}
    free(data_file);
}
