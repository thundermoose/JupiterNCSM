#ifndef __TRANSFORMED_BLOCK_MANAGER__
#define __TRANSFORMED_BLOCK_MANAGER__

#include <input/read_2nf_antoine_format.h>
#include <transform_block_settings/transform_block_settings.h>
#include <matrix_block_setting/matrix_block_setting.h>
#include <mercury_matrix_block/mercury_matrix_block.h>

struct _transform_2nf_block_manager_;
typedef struct _transform_2nf_block_manager_ *transform_2nf_block_manager_t;

transform_2nf_block_manager_t 
new_transform_2nf_block_manager(antoine_2nf_file_t coupled_2nf_data,
			      const char *basis_file_path,
			      const char *index_list_path,
			      int single_particle_energy_max);

void decouple_transform_2nf_block(transform_2nf_block_manager_t manager,
				  transform_block_settings_t block_settings);

mercury_matrix_block_t 
get_transform_2nf_matrix_block(transform_2nf_block_manager_t manager,
			       matrix_block_setting_t settings);

void free_transform_2nf_block_manager(transform_2nf_block_manager_t manager);

#endif
