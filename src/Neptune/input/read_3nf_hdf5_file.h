#ifndef __READ_3NF_HDF5_FILE__
#define __READ_3NF_HDF5_FILE__

#include <hdf5.h>
#include <bases/jjj_coupled_3p.h>
#include <bases/shells.h>
#include <matrix_transform/matrix_transform.h>
#include <input/read_3nf_file.h>
#include <utils/index_hash.h>
#include <omp.h>

typedef struct _block_configuration_
{
	int i,j;
} Block_Configuration;


typedef struct _block_
{
	int channel_number;
	index_hash_t configuration_to_index;
	double *matrix_elements;
	size_t num_matrix_elements;
	size_t min_conf_number;
	size_t max_conf_number;
	int score; // Number of threads that uses this block
	int in_use;
	omp_lock_t in_use_lock;
} Block;

typedef struct _channel_
{
	int J_abc, Tz, Parity, J_ab;
	/* J_ab is not used, 
	   but has to be included
	   due to the format of the
	   file*/
} Channel;

typedef struct _configuration_
{
	int a,b,c,J_ab;
	/* The abc are actually shell indices, however
	   in the hdf5 file format they are stored as int and not uint64*/
} Configuration;

typedef struct _hdf5_data_
{
	hid_t file_handle;
	Channel* channels;
	size_t num_channels;
	int e_max;
	Configuration* configurations;
	size_t num_configurations;
	size_t max_num_configurations;
	Block** open_blocks;
	omp_lock_t *block_locks;
	size_t loaded_memory;
	size_t max_loaded_memory;
	double weights[5];
	Shells* shells;
	// Workspace for load_channel
	double *out_array;
	size_t allocated_out_array_size;
	Block_Configuration *block_configurations;
	size_t allocated_configurations_size;
} HDF5_Data;

#define LEC_CE(data) data->weights[0]
#define LEC_CD(data) data->weights[1]
#define LEC_C1(data) data->weights[2]
#define LEC_C3(data) data->weights[3]
#define LEC_C4(data) data->weights[4]

HDF5_Data* open_hdf5_data(const char* file_name);

Dens_Matrix* get_matrix_hdf5(HDF5_Data* data_file,
			     JJJ_Basis* m_basis,
			     JJJ_Basis* n_basis);

int is_subset_of_basis_hdf5(JJJ_Basis* j_scheme,
			    HDF5_Data* datafile);


void free_hdf5_data(HDF5_Data* data_file);

#endif
