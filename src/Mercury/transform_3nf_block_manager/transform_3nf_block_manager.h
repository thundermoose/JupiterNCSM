#ifndef __TRANSFORM_3NF_BLOCK_MANAGER__
#define __TRANSFORM_3NF_BLOCK_MANAGER__

#include <input/read_3nf_file.h>
#include <matrix_block_setting/matrix_block_setting.h>
#include <transform_block_settings/transform_block_settings.h>
#include <mercury_matrix_block/mercury_matrix_block.h>

struct _tranform_2nf_block_manager_;
typedef struct _tranform_2nf_block_manager_ *transform_3nf_block_manager_t;

transform_3nf_block_manager_t 
new_transform_3nf_block_manager(Data_File *coupled_3nf_data,
				const char *index_list_path,
				const int single_particle_energy);

void decouple_transform_3nf_block(transform_3nf_block_manager_t manager,
				  transform_block_settings_t block);

mercury_matrix_block_t
get_transform_3nf_matrix_block(transform_3nf_block_manager_t manager,
			       matrix_block_setting_t block);

void free_transform_3nf_block_manager(transform_3nf_block_manager_t manager);

#endif
