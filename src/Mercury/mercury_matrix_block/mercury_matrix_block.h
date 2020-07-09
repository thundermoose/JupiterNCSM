#ifndef __MERCURY_MATRIX_BLOCK__
#define __MERCURY_MATRIX_BLOCK__

#include <interaction/interaction.h>
#include <connection_list/connection_list.h>
#include <single_particle_basis/single_particle_basis.h>

struct _mercury_matrix_block_;
typedef struct _mercury_matrix_block_ *mercury_matrix_block_t;

mercury_matrix_block_t new_mercury_matrix_block(interaction_t interaction,
						connection_list_t connections,
						single_particle_basis_t basis);

void save_mercury_matrix_block(mercury_matrix_block_t matrix_block,
			       const char *output_path);

void free_mercury_matrix_block(mercury_matrix_block_t matrix_block);

#endif
