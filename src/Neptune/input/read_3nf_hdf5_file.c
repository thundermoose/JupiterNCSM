#include <string.h>
#include <omp.h>
#include "read_3nf_hdf5_file.h"
#include <debug_mode/debug_mode.h>
#include <log/log.h>
#include <time.h>

static
Block *get_hdf5_block(HDF5_Data* data_file, int channel_number);

static
void release_block(Block *block);

static
size_t get_size_of_block(HDF5_Data* data_file, int channel_number);

static
void unload_low_score_channels(HDF5_Data *data_file,
			       size_t num_bytes_to_unload);

static
int compare_hdf5_blocks(Block **block_a, Block **block_b);

static
int is_block_in_use(Block *block);

static
int compare_configurations(Configuration *a,
			   Configuration *b);

static
uint64_t hash_configuration(Configuration *a);


int comp_confs(Configuration conf_a,
	       Configuration conf_b)
{
	return conf_a.a == conf_b.a &&
		conf_a.b == conf_b.b &&
		conf_a.c == conf_b.c &&
		conf_a.J_ab == conf_b.J_ab;
}


HDF5_Data* open_hdf5_data(const char* file_name)
{
	HDF5_Data* data_file =
		(HDF5_Data*)malloc(sizeof(HDF5_Data));
	data_file->file_handle = H5Fopen(file_name,
					 H5F_ACC_RDONLY,
					 H5P_DEFAULT);
	if (data_file->file_handle<0)
	{
		return NULL;
	}
	// Read channels
	hid_t channels =
		H5Dopen2(data_file->file_handle,
			 "/Channels",H5P_DEFAULT);

	data_file->num_channels =
		H5Dget_storage_size(channels)/sizeof(Channel);

	data_file->channels =
		(Channel*)malloc(sizeof(Channel)*
				 data_file->num_channels);

	H5Dread(channels,
		H5T_NATIVE_INT, H5S_ALL,
		H5S_ALL, H5P_DEFAULT,
		data_file->channels);
	H5Dclose(channels);
	// Read configurations
	hid_t configs =
		H5Dopen2(data_file->file_handle,
			 "/Configurations",H5P_DEFAULT);
	hid_t attrib =
		H5Aopen(configs,"cut_e3",H5P_DEFAULT);
	H5Aread(attrib, H5T_NATIVE_INT,&data_file->e_max);
	H5Aclose(attrib);

	data_file->num_configurations =
		H5Dget_storage_size(configs)/sizeof(Configuration);

	data_file->configurations =
		(Configuration*)malloc(sizeof(Configuration)*data_file->num_configurations);
	H5Dread(configs,
		H5T_NATIVE_INT, H5S_ALL,
		H5S_ALL,H5P_DEFAULT,
		data_file->configurations);
	H5Dclose(configs);
	// Initiate the block reading system, but not reading any blocks
	data_file->open_blocks = (Block**)calloc(data_file->num_channels,
						 sizeof(Block*));
	data_file->block_locks = (omp_lock_t*)malloc(data_file->num_channels*
						     sizeof(omp_lock_t));
	for (size_t i = 0; i<data_file->num_channels; i++)
		omp_init_lock(&data_file->block_locks[i]);
	data_file->loaded_memory = 0;
	data_file->max_loaded_memory = 80*((size_t)(1)<<30);
	data_file->weights[0] = -0.03957471;// CE
	data_file->weights[1] = 0.81680589;// CD
	data_file->weights[2] = -1.12152120;// C1
	data_file->weights[3] = -3.92500586;// C3
	data_file->weights[4] = 3.76568716;// C4
	data_file->out_array = NULL;
	data_file->allocated_out_array_size = 0;
	data_file->block_configurations = NULL;
	data_file->allocated_configurations_size = 0;
	return data_file;
}

void set_hdf5_max_loaded_memory(HDF5_Data *data_file,size_t max_loaded_memory)
{
	data_file->max_loaded_memory = max_loaded_memory;
}
			  

int get_config_number(HDF5_Data* data_file,
		      JJJ_State a,
		      hash_map_t state_to_index)
{
	Configuration ac = {a.a+1,a.b+1,a.c+1,a.j_ab};
	size_t index = 0;
	if (!hash_map_get(state_to_index, &ac,&index))
		return -1;
	return index;
}

int get_channel_number(HDF5_Data* data_file,
		       int J_abc,int Tz,int parity)
{
	int i;
	for (i = 0; i<data_file->num_channels; i++)
	{
		if (data_file->channels[i].J_abc == J_abc &&
		    data_file->channels[i].Tz == Tz &&
		    data_file->channels[i].Parity == parity)
		{
			return i;
		}
	}
	return -1;
}

int is_channel_loaded(HDF5_Data* data_file,
		      int channel_number)
{
	return data_file->open_blocks[channel_number] != NULL;
}



Block* load_channel(HDF5_Data* data_file,
		    int channel_number)
{
	struct timespec t_start,t_end;
	clock_gettime(CLOCK_REALTIME,&t_start);
	char channel_name[64];
	sprintf(channel_name, "Dataset_split_%.5d",channel_number+1);
	Block* block = (Block*)malloc(sizeof(Block));
	block->score = 0;
	block->channel_number = channel_number;
	Block_Configuration *configurations = NULL;
	double *out_array = NULL;
#pragma omp critical(hdf5)
	{
	hid_t group = H5Gopen(data_file->file_handle,
			      channel_name,H5P_DEFAULT);
	if (group<0)
	{
		printf("Tried and faild to open \"%s\"\n",
		       channel_name);
	}
	hid_t conf = H5Dopen2(group,"Configurations",H5P_DEFAULT);
	block->num_matrix_elements =
		H5Dget_storage_size(conf)/sizeof(Block_Configuration);
	configurations = 
		(Block_Configuration*)malloc(block->num_matrix_elements*
					     sizeof(Block_Configuration));
	H5Dread(conf,
		H5T_NATIVE_INT,
		H5S_ALL,
		H5S_ALL,
		H5P_DEFAULT,
		configurations);
	H5Dclose(conf);

	hid_t elem = H5Dopen2(group,"Elements",H5P_DEFAULT);
	size_t out_array_size = H5Dget_storage_size(elem)/sizeof(double);
	out_array = (double*)malloc(sizeof(double)*out_array_size);
	H5Dread(elem,
		H5T_NATIVE_DOUBLE,
		H5S_ALL,
		H5S_ALL,
		H5P_DEFAULT,
		out_array);
	H5Dclose(elem);
	}
	block->min_conf_number = configurations[0].i;
	block->max_conf_number = 
		configurations[block->num_matrix_elements-1].i;
	block->state_to_index = 
		new_hash_map(block->max_conf_number-block->min_conf_number+1,
			     3,
			     sizeof(size_t),
			     sizeof(Configuration),
			     (__compar_fn_t)compare_configurations,
			     (__hash_fn_t)hash_configuration);
	for (size_t i = block->min_conf_number-1; 
	     i < block->max_conf_number; 
	     i++)
	{
		size_t index = i+1;
		hash_map_insert(block->state_to_index,
				&data_file->configurations[i],
				&index);
	}
	block->matrix_elements =
	       	(double*)malloc(sizeof(double)*block->num_matrix_elements);
	size_t i;
	size_t len = block->num_matrix_elements;
	block->configuration_to_index = 
		new_index_hash(block->num_matrix_elements*2);
	for (i= 0; i<block->num_matrix_elements; i++)
	{
		block->matrix_elements[i] =
			(out_array[i]*data_file->weights[0]+
			 out_array[i+len]*data_file->weights[1]+
			 out_array[i+2*len]*data_file->weights[2]+
			 out_array[i+3*len]*data_file->weights[3]+
			 out_array[i+4*len]*data_file->weights[4]);
		set_index(block->configuration_to_index,
			  configurations[i].i,
			  configurations[i].j,
			  i);
	}
	free(out_array);
	free(configurations);
	omp_init_lock(&block->in_use_lock);
	data_file->open_blocks[channel_number] = block;
	clock_gettime(CLOCK_REALTIME,&t_end);
	double elapsed_time = (t_end.tv_sec-t_start.tv_sec)+
		(t_end.tv_nsec-t_start.tv_nsec)*1e-9;
	printf("Block %d needed %0.15lg s\n",
	       channel_number,elapsed_time);
	return block;
}


Block* get_channel(HDF5_Data* data_file,
		   int channel_number)
{
	for (size_t i = 0; i <= channel_number; i++)
		log_entry("data_file->open_blocks[%lu] = %p",
			  i,data_file->open_blocks[channel_number]);
	return data_file->open_blocks[channel_number];
}

double get_element(Block* block,
		   int i,int j)
{
	size_t index = get_index(block->configuration_to_index, i, j);
	if (index == no_index)
		return 0.0;
	else
		return block->matrix_elements[index];
}

void discard_channel(HDF5_Data* data_file,
		     Block* channel)
{
	if (data_file->open_blocks[channel->channel_number] != channel)
	{
		fprintf(stderr,"This is not supposed to happen\n");
		fprintf(stderr,"channel: %d\n",channel->channel_number);
	}
	data_file->open_blocks[channel->channel_number] = NULL;
	omp_destroy_lock(&channel->in_use_lock);
	free_index_hash(channel->configuration_to_index);
	free_hash_map(channel->state_to_index);
	free(channel->matrix_elements);
	free(channel);
}



Dens_Matrix* get_matrix_hdf5(HDF5_Data* data_file,
			     JJJ_Basis* m_basis,
			     JJJ_Basis* n_basis)
{

	Dens_Matrix *mat =
		new_zero_matrix(m_basis->dimension,
				n_basis->dimension);

	Block* current_block = NULL;
	size_t i,j;
	for (i = 0; i<m_basis->dimension; i++)
	{
		int needed_chan =
			get_channel_number(data_file,
					   m_basis->states[i].j_abc,
					   m_basis->states[i].tz,
					   m_basis->states[i].parity);
		if (needed_chan <0)
		{
			fprintf(stderr,"Could not find channel %d %d %d\n",
				m_basis->states[i].j_abc,
				m_basis->states[i].tz,
				m_basis->states[i].parity);
			exit(1);
		}

		if (current_block == NULL ||
		    current_block->channel_number != needed_chan)
		{
			log_entry("needed_chan = %d",needed_chan);
			if (current_block)
			{
				log_entry("current_block->channel_number = %d", 
					  current_block->channel_number);
			}
			else
			{
				log_entry("current_block = NULL");
			}
			release_block(current_block);
			current_block =
				get_hdf5_block(data_file,needed_chan);
		}
		if (current_block == NULL)
			continue;

		int i_conf_num = 
			get_config_number(data_file,
					  m_basis->states[i],
					  current_block->state_to_index);
		if (i_conf_num<0)
			continue;
		int ic = i_conf_num;
		for (j = 0; j<n_basis->dimension; j++)
		{ 
			if (n_basis->states[j].j_abc != m_basis->states[i].j_abc ||
			    n_basis->states[j].tz != m_basis->states[i].tz ||
			    n_basis->states[j].parity != m_basis->states[i].parity)
				continue;

			int j_conf_num = 
				get_config_number(data_file,
						  n_basis->states[j],
						  current_block->state_to_index);


			if (j_conf_num<0)
				continue;

			if (i_conf_num>j_conf_num)
			{
				int t = i_conf_num;
				i_conf_num = j_conf_num;
				j_conf_num = t;
			}

			mat->elements[i*mat->n+j] = get_element(current_block,
								i_conf_num,
								j_conf_num);

			i_conf_num = ic;
		}
	}
	if (current_block != NULL)
		release_block(current_block);
	return mat;
}


int is_subset_of_basis_hdf5(JJJ_Basis* j_scheme,
			    HDF5_Data* datafile)
{
	size_t i,j;
	for (i = 0; i<j_scheme->dimension; i++)
	{
		JJJ_State js = j_scheme->states[i];
		int exists = 0;
		for (j = 0; j<datafile->num_configurations; j++)
		{
			Configuration conf = datafile->configurations[j];
			if (js.a+1 == conf.a &&
			    js.b+1 == conf.b &&
			    js.c+1 == conf.c &&
			    js.j_ab == conf.J_ab)
			{
				exists = 1;
				break;
			}
		}
		if (!exists)
		{ // if atleast one of the states in j_scheme do not exists
			// in datafile->configurations j_scheme is not a subset

			log_entry("Could not find state (%ld) (%ld %ld %ld | %d)\n",
				  i,js.a,js.b,js.c,js.j_ab);
			return 0;
		}
	}
	return 1;
}


void free_hdf5_data(HDF5_Data* data_file)
{
	H5Fclose(data_file->file_handle);
	free(data_file->channels);
	free(data_file->configurations);
	size_t i;
	for (i = 0; i<data_file->num_channels; i++)
	{
		if (data_file->open_blocks[i] != NULL)
		{
			discard_channel(data_file,
					data_file->open_blocks[i]);
		}
		omp_destroy_lock(&data_file->block_locks[i]);
	}
	free(data_file->open_blocks);
	free(data_file->block_locks);
	free(data_file->out_array);
	free(data_file->block_configurations);
	free(data_file);
}

static
Block *get_hdf5_block(HDF5_Data* data_file, int channel_number)
{
	Block *block = NULL;
	omp_set_lock(&data_file->block_locks[channel_number]);
	if (is_channel_loaded(data_file,channel_number))
	{
		block = get_channel(data_file, channel_number);
		if (block)
		{
			block->score++;
			omp_set_lock(&block->in_use_lock);
			block->in_use++;
			omp_unset_lock(&block->in_use_lock);
		}
	}
	else
	{
		size_t size_of_needed_block =
		       	get_size_of_block(data_file, channel_number);
		if (data_file->loaded_memory + size_of_needed_block > 
		    data_file->max_loaded_memory)
			unload_low_score_channels(data_file,
						  size_of_needed_block);
		data_file->loaded_memory+=size_of_needed_block;
		block = load_channel(data_file,channel_number);
		block->score=10;
		omp_set_lock(&block->in_use_lock);
		block->in_use=1;
		omp_unset_lock(&block->in_use_lock);
	}
	omp_unset_lock(&data_file->block_locks[channel_number]);
	return block;
}

static
void release_block(Block *block)
{
	if (block)
	{
		omp_set_lock(&block->in_use_lock);
		block->in_use--;
		omp_unset_lock(&block->in_use_lock);
	}
}

static
size_t get_size_of_block(HDF5_Data *data_file, int channel_number)
{
	char channel_name[64];
	sprintf(channel_name, "Dataset_split_%.5d",channel_number+1);
	size_t size_of_block = 0;
#pragma omp critical (hdf5)
	{
	hid_t group = H5Gopen(data_file->file_handle,
			      channel_name,H5P_DEFAULT);
	if (group<0)
	{
		printf("Tried and faild to open \"%s\"\n",
		       channel_name);
		exit(EXIT_FAILURE);
	}
	hid_t elements = H5Dopen2(group,"Elements",H5P_DEFAULT);
	size_t total_size = sizeof(double) + sizeof(Block_Configuration);
	size_of_block = H5Dget_storage_size(elements)*total_size;
	H5Dclose(elements);
	H5Gclose(group);
	}
	return size_of_block/10;
}

static
void unload_low_score_channels(HDF5_Data *data_file,size_t num_bytes_to_unload)
{
#pragma omp critical (unload_blocks)
	{
		Block **blocks =
			(Block**)malloc(data_file->num_channels*
					sizeof(Block*));
		memcpy(blocks,
		       data_file->open_blocks,
		       data_file->num_channels*sizeof(Block*));
		qsort(blocks,
		      data_file->num_channels,
		      sizeof(Block*),
		      (__compar_fn_t)compare_hdf5_blocks);
		size_t unloaded_memory = 0;
		for (size_t i = 0; i<data_file->num_channels; i++)
		{
			if (blocks[i] == NULL)
				break;
			int channel_number = blocks[i]->channel_number;
			omp_set_lock(&data_file->block_locks[channel_number]);
			if (!is_block_in_use(blocks[i]))
			{
				unloaded_memory +=
					get_size_of_block(data_file,
							  channel_number);
				discard_channel(data_file, blocks[i]);
			}
			omp_unset_lock(&data_file->block_locks[channel_number]);
			if (unloaded_memory >= num_bytes_to_unload)
				break;
		}
		free(blocks);
		data_file->loaded_memory -= unloaded_memory;
	}
}

static
int compare_hdf5_blocks(Block **block_a, Block **block_b)
{
	if (*block_a != NULL && *block_b == NULL)
		return -1;
	else if (*block_a == NULL && *block_b != NULL)
		return 1;
	else if (*block_a == NULL && *block_b == NULL)
		return 0;
	else if ((*block_a)->num_matrix_elements < 
		 (*block_b)->num_matrix_elements)
		return -1;
	else if ((*block_a)->num_matrix_elements > 
		 (*block_b)->num_matrix_elements)
		return 1;
	else if ((*block_a)->score < (*block_b)->score)
		return -1;
	else if ((*block_a)->score > (*block_b)->score)
		return 1;
	else
		return 0;	
}

static
int is_block_in_use(Block *block)
{
	int is_in_use = 0;
	omp_set_lock(&block->in_use_lock);
	is_in_use = block->in_use > 0;
	omp_unset_lock(&block->in_use_lock);
	return is_in_use;
}

static
int compare_configurations(Configuration *a,
			   Configuration *b)
{
	return !comp_confs(*a,*b);
}

static
uint64_t hash_configuration(Configuration *a)
{
	uint64_t h1 = (uint64_t)(a->a) ^ ((uint64_t)(a->b)<<13);
	uint64_t h2 = (uint64_t)(a->b) ^ ((uint64_t)(a->c)<<17);
	uint64_t h3 = (uint64_t)(a->c) ^ ((uint64_t)(a->J_ab)<<31);
	return h1 ^ __builtin_bswap64(h2) ^ h3;	
}
