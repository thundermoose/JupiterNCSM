#ifndef __CONNECTION_LIST__
#define __CONNECTION_LIST__

#include <matrix_block_setting/matrix_block_setting.h>
#include <single_particle_basis/single_particle_basis.h>
#include <block_type/block_type.h>

struct _connection_list_;
typedef struct _connection_list_ *connection_list_t;

typedef struct
{
	block_type_t type;
	short neutron_states[6];
	short proton_states[6];	
} connection_t;

connection_list_t new_connection_list(const char *index_list_path,
				      single_particle_basis_t basis,
				      matrix_block_setting_t settings);

connection_list_t read_connection_files(const char *index_list_path,
					matrix_block_setting_t settings);

matrix_block_setting_t get_matrix_block_setting(connection_list_t connections);

size_t num_connections(connection_list_t connection_list);

int has_next_connection(connection_list_t connection_list);

connection_t next_connection(connection_list_t connection_list);

void free_connection_list(connection_list_t connection_list);

int compute_M(connection_t connection,
	      single_particle_basis_t basis);

#endif
