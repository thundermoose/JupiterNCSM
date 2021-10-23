#include <connection_list/connection_list.h>
#include <debug_mode/debug_mode.h>
#include <error/error.h>
#include <log/log.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

struct _connection_list_
{
	matrix_block_setting_t settings;
	connection_t *connections;
	size_t num_connections;
	size_t index;
};

static
void read_connection_file(short **connections,
			  size_t *num_connections,
			  size_t num_particles,
			  const char *index_list_path,
			  char prefix,
			  int differnece_energy,
			  int difference_M,
			  int depth);

static
size_t short_to_size(short number)
{
	union{
		size_t final_value;
		short initial_values[4];
	} convert = {0};
	convert.initial_values[0] = number;
	return convert.final_value;
}

connection_list_t new_connection_list(const char *index_list_path,
				      matrix_block_setting_t settings)
{
	block_type_t type = settings.type;
	size_t num_protons = count_protons(type);
	char *proton_directory = 
		num_protons > 0 ?
		(char*)malloc(strlen(index_list_path)+256) :
		NULL;		
	switch (num_protons)
	{
		case 1:
			sprintf(proton_directory,
				"%s/p_inds_index_lists"
				index_list_path);
			break;
		case 2:
			sprintf(proton_directory,
				"%s/pp_inds_index_lists"
				index_list_path);
			break;
		case 3:
			sprintf(proton_directory,
				"%s/ppp_inds_index_lists"
				index_list_path);
			break;
		default:
			;
	}
}

connection_list_t read_connection_files(const char *index_list_path,
					matrix_block_setting_t settings)
{
	block_type_t type = settings.type;
	size_t num_neutrons = count_neutrons(type);
	size_t num_protons = count_protons(type);
	short *neutron_connections = NULL;
	size_t num_neutron_connections = 0;
	short *proton_connections = NULL;
	size_t num_proton_connections = 0;
	if (num_neutrons > 0)
		read_connection_file(&neutron_connections,
				     &num_neutron_connections,
				     num_neutrons,
				     index_list_path,
				     'n',
				     settings.difference_energy_neutrons,
				     settings.difference_M_neutrons,
				     settings.depth_neutrons);		
	if (num_protons > 0)
		read_connection_file(&proton_connections,
				     &num_proton_connections,
				     num_protons,
				     index_list_path,
				     'p',
				     settings.difference_energy_protons,
				     settings.difference_M_protons,
				     settings.depth_protons);		
	size_t num_connections =
		num_proton_connections*num_neutron_connections;
	if (num_proton_connections == 0)
		num_connections = num_neutron_connections;
	else if (num_neutron_connections == 0)
		num_connections = num_proton_connections;
	connection_t *connections =
		(connection_t*)malloc(num_connections*sizeof(connection_t));
	connection_t current_connection = {0};
	current_connection.type = type;
	if (num_proton_connections > 0 && num_neutron_connections > 0)
	{
		log_entry("Generating np connections");
		log_entry("num_proton_connections = %lu",
			  num_proton_connections);
		log_entry("num_neutron_connections = %lu",
			  num_neutron_connections);
		size_t connection_index = 0;
		for (size_t i = 0; i<num_proton_connections; i++)
		{
			log_entry("Proton connection %lu",i);
			memcpy(current_connection.proton_states,
			       proton_connections+2*num_protons*i,
			       2*num_protons*sizeof(short));
			log_entry("Proton connection: (%d %d %d %d %d %d)",
				  current_connection.proton_states[0],
				  current_connection.proton_states[1],
				  current_connection.proton_states[2],
				  current_connection.proton_states[3],
				  current_connection.proton_states[4],
				  current_connection.proton_states[5]);
			for (size_t j = 0; j<num_neutron_connections; j++)
			{
				log_entry("Neutron connection %lu",j);
				log_entry("Neutron connection: (%d %d %d %d %d %d)",
					  current_connection.neutron_states[0],
					  current_connection.neutron_states[1],
					  current_connection.neutron_states[2],
					  current_connection.neutron_states[3],
					  current_connection.neutron_states[4],
					  current_connection.neutron_states[5]);
				memcpy(current_connection.neutron_states,
				       neutron_connections+2*num_neutrons*j,
				       2*num_neutrons*sizeof(short));
				log_entry("connection_index = %lu",
					  connection_index);
				connections[connection_index++] = current_connection;

			}
		}
		log_entry("All connections:");
		for (size_t i = 0; i<connection_index; i++)
		{
			log_entry("connection (%lu) np (%d %d %d %d %d %d), (%d %d %d %d %d %d)",
				  i,
				  connections[i].neutron_states[0],
				  connections[i].neutron_states[1],
				  connections[i].neutron_states[2],
				  connections[i].neutron_states[3],
				  connections[i].neutron_states[4],
				  connections[i].neutron_states[5],
				  connections[i].proton_states[0],
				  connections[i].proton_states[1],
				  connections[i].proton_states[2],
				  connections[i].proton_states[3],
				  connections[i].proton_states[4],
				  connections[i].proton_states[5]);


		}
		free(neutron_connections);
		free(proton_connections);
	}
	else if (num_neutron_connections > 0)
	{
		for (size_t i = 0; i < num_neutron_connections; i++)
		{
			memcpy(current_connection.neutron_states,
			       neutron_connections+2*num_neutrons*i,
			       2*num_neutrons*sizeof(short));
			connections[i] = current_connection;
		}	
		free(neutron_connections);
	}
	else if (num_proton_connections > 0)
	{
		for (size_t i = 0; i < num_proton_connections; i++)
		{
			memcpy(current_connection.proton_states,
			       proton_connections+2*num_protons*i,
			       2*num_protons*sizeof(short));
			connections[i] = current_connection;
		}	
		free(proton_connections);
	}
	connection_list_t connection_list =
	       	(connection_list_t)malloc(sizeof(struct _connection_list_));
	connection_list->settings = settings;
	connection_list->connections = connections;
	connection_list->num_connections = num_connections;
	connection_list->index = 0;
	return connection_list;
}

matrix_block_setting_t get_matrix_block_setting(connection_list_t connections)
{
	return connections->settings;
}

size_t num_connections(connection_list_t connection_list)
{
	return connection_list->num_connections;
}

int has_next_connection(connection_list_t connection_list)
{
	return connection_list->index < connection_list->num_connections;
}

connection_t next_connection(connection_list_t connection_list)
{
	return connection_list->connections[connection_list->index++];
}

void free_connection_list(connection_list_t connection_list)
{
	free(connection_list->connections);
	free(connection_list);	
}

int compute_M(connection_t connection,
	      single_particle_basis_t basis)
{
	const size_t num_neutrons = count_neutrons(connection.type);
	const size_t num_protons = count_protons(connection.type);
	int M = 0;
	for (size_t i = 0; i<num_neutrons; i++)
	{
		size_t index = short_to_size(connection.neutron_states[i]*2+1);
		single_particle_state_t state = get_state(basis,index);
		M += state.m;
	}
	for (size_t i = 0; i<num_protons; i++)
	{
		size_t index = short_to_size(connection.proton_states[i]*2);
		single_particle_state_t state = get_state(basis,index);
		M += state.m;
	}
	return M;
}

	static
void read_connection_file(short **connections,
			  size_t *num_connections,
			  size_t num_particles,
			  const char *index_list_path,
			  char prefix,
			  int differnece_energy,
			  int difference_M,
			  int depth)
{
	const char *postfix = "_inds_index_lists";
	const char *file_pattern = "/conn_dE_%d_dM_%d_depth_%d";
	const size_t index_list_path_length = strlen(index_list_path);
	const size_t postfix_length = strlen(postfix);
	size_t file_name_length =
		index_list_path_length+
		num_particles+
		postfix_length+
		strlen(file_pattern)+17;
	char *file_name = (char*)malloc(file_name_length);
	memcpy(file_name,
	       index_list_path,
	       index_list_path_length);
	file_name[index_list_path_length] = '/';
	memset(file_name+index_list_path_length+1,
	       prefix,
	       num_particles);
	memcpy(file_name+index_list_path_length+num_particles+1,
	       postfix,
	       postfix_length);
	sprintf(file_name+index_list_path_length+num_particles+postfix_length+1,
		file_pattern,
		differnece_energy,
		difference_M,
		depth);	
	FILE *connection_file = fopen(file_name,"r");
	if (connection_file == NULL)
		error("Could not open file %s. %s\n",
		      file_name,
		      strerror(errno));
	fseek(connection_file,0,SEEK_END);
	size_t num_bytes = ftell(connection_file);
	fseek(connection_file,0,SEEK_SET);
	if (num_bytes == 0)
		return;
	if (num_bytes % sizeof(short) != 0)
		error("File %s do not contain shorts\n",
		      file_name);
	if (num_bytes % (2*num_particles) != 0)
		error("File %s do not contain %lu particle connections\n",
		      file_name,num_particles);
	*num_connections = num_bytes / (sizeof(short)*2*num_particles);
	*connections = (short*)malloc(num_bytes);
	if (fread(*connections,
		  sizeof(short)*2*num_particles,
		  *num_connections,
		  connection_file) != *num_connections)
		error("Could not read connections from %s\n",
		      file_name);
	free(file_name);
	fclose(connection_file);
}
