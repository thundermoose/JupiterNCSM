#include <string.h>
#include <math.h>
#include <output/write_to_hdf5.h>
#include <utils/questions.h>
#include <debug_mode/debug_mode.h>

void save_sp_basis(hid_t file_handle,
		   SP_States *sp_states)

{
	const hsize_t sp_basis_size[2] = {sp_states->dimension,5};
	hid_t sp_basis = H5Dcreate2(file_handle,
				    "Basis",
				    H5T_NATIVE_INT,
				    H5Screate_simple(2,
						     sp_basis_size,
						     NULL),
				    H5P_DEFAULT,
				    H5P_DEFAULT,
				    H5P_DEFAULT);
	Shells *shells = sp_states->shells;
	int* sp_array = (int*)malloc(sizeof(int)*sp_states->dimension*5);
	size_t i;
	for (i = 0; i<sp_states->dimension; i++)
	{
		sp_array[5*i]= shells->shells[sp_states->sp_states[i].shell].n;
		sp_array[5*i+1]= shells->shells[sp_states->sp_states[i].shell].l;
		sp_array[5*i+2]= shells->shells[sp_states->sp_states[i].shell].j;
		sp_array[5*i+3]= shells->shells[sp_states->sp_states[i].shell].tz;
		sp_array[5*i+4]= sp_states->sp_states[i].m;
	}
	H5Dwrite(sp_basis,
		 H5T_NATIVE_INT,
		 H5S_ALL,
		 H5S_ALL,
		 H5P_DEFAULT,
		 sp_array);
	H5Dclose(sp_basis);
	free(sp_array);
}

void save_configurations(hid_t file_handle,
			 M_Scheme_3p_Basis *ms_basis)
{
	const hsize_t config_size[2] = {ms_basis->dimension,3};
	hid_t configs = H5Dcreate2(file_handle,
				   "Configurations",
				   H5T_NATIVE_INT,
				   H5Screate_simple(2,
						    config_size,
						    NULL),
				   H5P_DEFAULT,
				   H5P_DEFAULT,
				   H5P_DEFAULT);
	int* array = (int*)malloc(sizeof(int)*3*ms_basis->dimension);
	size_t i;
	for (i = 0; i<ms_basis->dimension; i++)
	{
		array[3*i] = ms_basis->states[i].a;
		array[3*i+1] = ms_basis->states[i].b;
		array[3*i+2] = ms_basis->states[i].c;
	}
	H5Dwrite(configs,
		 H5T_NATIVE_INT,
		 H5S_ALL,
		 H5S_ALL,
		 H5P_DEFAULT,
		 array);
	H5Dclose(configs);
	free(array);
}

HDF5_Out_File* create_new_hdf5_out_file(const char* file_name,
					M_Scheme_3p_Basis *ms_basis)
{
	hid_t file_handle = H5Fcreate(file_name,H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
	if (file_handle<0)
	{
		char question[256];
		sprintf(question,
			"File \"%s\" already exists, do you want to overwrite it?",
			file_name);
		if (yes_no_question(question)>0)
		{
			file_handle = H5Fcreate(file_name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
		}
		else
		{
			exit(1);
		}


	}
	save_sp_basis(file_handle,
		      ms_basis->sp_states);
	save_configurations(file_handle,
			    ms_basis);
	HDF5_Out_File* out = (HDF5_Out_File*)malloc(sizeof(HDF5_Out_File));
	out->file = file_handle;
	out->data_blocks = NULL;
	out->current_number_blocks = 0;
	out->max_number_blocks = 0;
	return out;
}




size_t add_hdf5_block(HDF5_Out_File* file,
		      quantum_number Tz,
		      quantum_number M,
		      quantum_number E1,
		      quantum_number E2)
{
#pragma omp critical(output)
	{
		if (file->current_number_blocks == file->max_number_blocks)
		{
			file->max_number_blocks = 2*file->max_number_blocks+1;
			file->data_blocks =
				(Data_Block*)realloc(file->data_blocks,
						     sizeof(Data_Block)*file->max_number_blocks);
		}
		file->data_blocks[file->current_number_blocks].Tz = Tz;
		file->data_blocks[file->current_number_blocks].M = M;
		file->data_blocks[file->current_number_blocks].E1 = E1;
		file->data_blocks[file->current_number_blocks].E2 = E2;
		file->current_number_blocks++;
	}
	return file->current_number_blocks;
}



void write_to_hdf5_block(HDF5_Out_File* file,
			 size_t block_number,
			 Dens_Matrix elements,
			 size_t *m_indices,
			 size_t *n_indices)
{
#pragma omp critical(hdf5)
	{


		// determine if the block is triangular
		int is_triangular =
			file->data_blocks[block_number-1].E1 == file->data_blocks[block_number-1].E2;

		// generate config array and element array
		size_t max_num_elements;
		if (is_triangular){
			max_num_elements = (elements.m*(elements.m+1))/2;
		}else{
			max_num_elements = elements.m*elements.n;
		}
		int *configs = (int*)malloc(sizeof(int)*2*max_num_elements);
		double *mat_elements = (double*)malloc(sizeof(double)*max_num_elements);
		double treshold = 1e-20;
		size_t num_elements = 0;
		size_t i,j;
		if (is_triangular)
		{
			for (i = 0; i<elements.m; i++)
			{

				for (j = i; j<elements.n; j++)
				{
					if (fabs(elements.elements[i*elements.n+j])>treshold)
					{
						mat_elements[num_elements] = elements.elements[i*elements.m+j];
						configs[2*num_elements] = m_indices[i];
						configs[2*num_elements+1] = n_indices[j];
						num_elements++;
					}
				}
			}
		}
		else
		{
			for (i = 0; i<elements.m*elements.n; i++)
			{
				if (fabs(elements.elements[i])>treshold)
				{
					mat_elements[num_elements] = elements.elements[i];
					configs[2*num_elements] = m_indices[i/elements.n];
					configs[2*num_elements+1] = n_indices[i%elements.n];
					num_elements++;
				}
			}
		}

		if (num_elements != 0)
		{

			// Create group
			char name[256];
			sprintf(name,
				"DataSplit_%05ld",
				block_number);
			hid_t group = H5Gcreate2(file->file,
						 name,
						 H5P_DEFAULT,
						 H5P_DEFAULT,
						 H5P_DEFAULT);

			// Store configurations
			const hsize_t config_dimension[2]={num_elements,2};
			hid_t config_data =
				H5Dcreate2(group,
					   "Configurations",
					   H5T_NATIVE_INT,
					   H5Screate_simple(2,
							    config_dimension,
							    config_dimension),
					   H5P_DEFAULT,
					   H5P_DEFAULT,
					   H5P_DEFAULT);

			H5Dwrite(config_data,
				 H5T_NATIVE_INT,
				 H5S_ALL,
				 H5S_ALL,
				 H5P_DEFAULT,
				 configs);

			H5Dclose(config_data);
			free(configs);


			// Store elements

			const hsize_t element_dimension[1] = {num_elements};

			hid_t element_data =
				H5Dcreate2(group,
					   "Elements",
					   H5T_NATIVE_DOUBLE,
					   H5Screate_simple(1,
							    element_dimension,
							    element_dimension),
					   H5P_DEFAULT,
					   H5P_DEFAULT,
					   H5P_DEFAULT);

			H5Dwrite(element_data,
				 H5T_NATIVE_DOUBLE,
				 H5S_ALL,
				 H5S_ALL,
				 H5P_DEFAULT,
				 mat_elements);

			H5Dclose(element_data);
			H5Gclose(group);
		}
		free(mat_elements);
	}
}

	static
void remove_unused_channels(HDF5_Out_File* file)
{
	// This method is going through the blocks
	// and checks if such a block exists in
	// the hdf5 file
	size_t block_index; 
	size_t used_block_index; 
	for (block_index = used_block_index = 0;
	     block_index < file->current_number_blocks;
	     block_index++)
	{
		// check if the block exists
		char block_name[256];
		sprintf(block_name,
			"DataSplit_%05ld",
			block_index+1);
		if (H5Lexists(file->file,
			      block_name,
			      H5Pcreate(H5P_LINK_ACCESS))<=0)
		{
			// link do not exists
			continue;
		}
		// link do exists

		// move the block to the correct position
		file->data_blocks[used_block_index] = file->data_blocks[block_index];

		char new_block_name[256];
		sprintf(new_block_name,
			"DataSplit_%05ld",
			used_block_index+1);
		printf("%s becomes %s\n",block_name,new_block_name);
		if (H5Lmove(file->file,
			    block_name,
			    file->file,
			    new_block_name,
			    H5Pcreate(H5P_LINK_CREATE),
			    H5Pcreate(H5P_LINK_ACCESS))<0)
		{
			fprintf(stderr,"Could not move link \"%s\" to \"%s\"\n",
				block_name,new_block_name);
			exit(1);
		}

		used_block_index++;
	}
	file->current_number_blocks = file->max_number_blocks = used_block_index;
	file->data_blocks = (Data_Block*)realloc(file->data_blocks,
						 sizeof(Data_Block)*
						 file->max_number_blocks);
}

	static
void save_channels(HDF5_Out_File* file)
{
	const hsize_t channels_dimension[2] =
	{file->current_number_blocks,4};
	hid_t channels = H5Dcreate2(file->file,
				    "Channels",
				    H5T_NATIVE_INT,
				    H5Screate_simple(2,
						     channels_dimension,
						     channels_dimension),
				    H5P_DEFAULT,
				    H5P_DEFAULT,
				    H5P_DEFAULT);

	H5Dwrite(channels,
		 H5T_NATIVE_INT,
		 H5S_ALL,
		 H5S_ALL,
		 H5P_DEFAULT,
		 file->data_blocks);


	H5Fclose(file->file);
}

void close_hdf5_out_file(HDF5_Out_File* file)
{
	remove_unused_channels(file);
	save_channels(file);

	if (file->data_blocks != NULL)
	{
		free(file->data_blocks);
	}
	free(file);
}
