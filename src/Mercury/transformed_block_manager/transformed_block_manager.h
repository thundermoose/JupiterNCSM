#ifndef __TRANSFORMED_BLOCK_MANAGER__
#define __TRANSFORMED_BLOCK_MANAGER__

#include <input/read_2nf_antoine_format.h>
#include <transform_block_settings/transform_block_settings.h>
#include <matrix_block_setting/matrix_block_setting.h>
#include <mercury_matrix_block/mercury_matrix_block.h>

struct _transformed_block_manager_;
typedef struct _transformed_block_manager_ *transformed_block_manager_t;

transformed_block_manager_t 
new_transformed_block_manager(antoine_2nf_file_t coupled_2nf_data);

void decouple_transform_block(transformed_block_manager_t manager,
			      transform_block_settings_t block_settings);

mercury_matrix_block_t 
get_transformed_matrix_block(transformed_block_manager_t manager,
			     matrix_block_setting_t settings);

void free_transformed_block_manager(transformed_block_manager_t manager);

#endif
