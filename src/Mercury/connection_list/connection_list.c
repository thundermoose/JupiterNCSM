#include <connection_list/connection_list.h>
#include <debug_mode/debug_mode.h>
#include <array_builder/array_builder.h>
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
char *create_directory_path(const char *index_list_path,
			    size_t num_particles,
			    char particle_type);

static
short *create_single_particle_connections(const char *directory,
					  size_t num_particles,
					  int difference_energy,
					  int difference_M,
					  int depth,
					  single_particle_basis_t basis,
					  size_t *num_connections);

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
int single_species_compute_M(single_particle_basis_t basis,
			     short *state,
			     size_t num_particles);

static
short *load_states(const char *directory,int E,size_t *length);

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
				      single_particle_basis_t basis,
				      matrix_block_setting_t settings)
{
	block_type_t type = settings.type;
	size_t num_protons = count_protons(type);
	size_t num_neutrons = count_neutrons(type);
	char *proton_directory = 
		create_directory_path(index_list_path,
				      num_protons,
				      'p');
	char *neutron_directory = 
		create_directory_path(index_list_path,
				      num_neutrons,
				      'n');
	size_t num_proton_connections = 0;
	short *proton_connections =
		create_single_particle_connections(proton_directory,
						   num_protons,				   
						   settings.difference_energy_protons,
						   settings.difference_M_protons,
						   settings.depth_protons,
						   basis,
						   &num_proton_connections);
	size_t num_neutron_connections = 0;
	short *neutron_connections =
		create_single_particle_connections(neutron_directory,
						   num_neutrons,
						   settings.difference_energy_neutrons,
						   settings.difference_M_neutrons,
						   settings.depth_neutrons,
						   basis,
						   &num_neutron_connections);
	connection_list_t list =
		(connection_list_t)calloc(1,sizeof(struct _connection_list_));
	list->settings = settings;
	if (num_protons == 0)
	{
		list->num_connections = num_neutron_connections;	
		list->connections = 
			(connection_t*)
			malloc(num_neutron_connections*sizeof(connection_t));
		connection_t current_connection;
		current_connection.type = type;
		for (size_t neutron_index = 0; 
		     neutron_index < num_neutron_connections; 
		     neutron_index++)
		{
			memcpy(current_connection.neutron_states,
			       neutron_connections+2*neutron_index*num_neutrons,
			       2*num_neutrons*sizeof(short));
			list->connections[neutron_index] = current_connection;
		}
	}	
	else if (num_neutrons == 0)
	{
		list->num_connections = num_proton_connections;	
		list->connections = 
			(connection_t*)
			malloc(num_proton_connections*sizeof(connection_t));
		connection_t current_connection;
		current_connection.type = type;
		for (size_t proton_index = 0; 
		     proton_index < num_proton_connections; 
		     proton_index++)
		{
			memcpy(current_connection.proton_states,
			       proton_connections+2*proton_index*num_protons,
			       2*num_protons*sizeof(short));
			list->connections[proton_index] = current_connection;
		}
	} 
	else
	{
		list->num_connections = 
			num_proton_connections*num_neutron_connections;	
		list->connections = 
			(connection_t*)
			malloc(list->num_connections*sizeof(connection_t));
		connection_t current_connection;
		current_connection.type = type;
		size_t connection_index = 0;
		for (size_t proton_index = 0; 
		     proton_index < num_proton_connections; 
		     proton_index++)
		{
			memcpy(current_connection.proton_states,
			       proton_connections+2*proton_index*num_protons,
			       2*num_protons*sizeof(short));
			for (size_t neutron_index = 0;
			     neutron_index < num_neutron_connections;
			     neutron_index++)
			{
				memcpy(current_connection.neutron_states,
				       neutron_connections+
				       2*neutron_index*num_neutrons,
				       2*num_neutrons*sizeof(short));
				list->connections[connection_index++] = 
					current_connection;
			}
		}
		if (list->num_connections != settings.num_proton_combinations*settings.num_neutron_combinations)
		{
			fprintf(stderr,"matrix block %lu should have %lu num connections but only has %lu\n",
				settings.matrix_block_id,
				settings.num_proton_combinations*settings.num_neutron_combinations,
				list->num_connections);
			fprintf(stderr,"num_proton_connections = %lu\n",
				num_proton_connections);
			fprintf(stderr,"num_neutron_connections = %lu\n",
				num_neutron_connections);
			exit(1);
		}
	}
	free(neutron_connections);
	free(proton_connections);
	free(neutron_directory);
	free(proton_directory);
	return list;
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
char *create_directory_path(const char *index_list_path,
			    size_t num_particles,
			    char particle_type)
{
	size_t length = strlen(index_list_path)+256;	
	char *directory_path = (char*)calloc(length,sizeof(char));
	char directory_name[] = "   _inds_index_lists";
	memset(directory_name,particle_type,3);
	sprintf(directory_path,
		"%s/%s",
		index_list_path,
		directory_name+3-num_particles);
	return directory_path;
}

	static
short *create_single_particle_connections(const char *directory,
					  size_t num_particles,
					  int difference_energy,
					  int difference_M,
					  int depth,
					  single_particle_basis_t basis,
					  size_t *num_connections)
{
	if (num_particles == 0)
		return NULL;
	int annihilated_energy = depth;
	int created_energy = depth+difference_energy;
	size_t length_annihilated_states = 0;
	short *annihilated_states = 
		load_states(directory,
			    annihilated_energy,
			    &length_annihilated_states);
	length_annihilated_states/=num_particles;
	size_t length_created_states = 0;
	short *created_states = 
		load_states(directory,
			    created_energy,
			    &length_created_states);
	length_created_states/=num_particles;
	short *single_particle_connections = NULL;
	array_builder_t connections_builder = 
		new_array_builder((void**)&single_particle_connections,
				  num_connections,
				  2*num_particles*sizeof(short));
	for (size_t annihilated_index = 0;
	     annihilated_index < length_annihilated_states;
	     annihilated_index++)
	{
		int annihilated_M = 
			single_species_compute_M(basis,
						 annihilated_states
						 +annihilated_index*num_particles,
						 num_particles);
		int created_M = annihilated_M+difference_M;
		size_t created_index = 0;
		for (;created_index < length_created_states; created_index++)
			if (single_species_compute_M(basis,
						     created_states
						     +created_index*num_particles,
						     num_particles) == created_M)
				break;
		for (;created_index < length_created_states &&
		     (single_species_compute_M(basis,
					       created_states
					       +created_index*num_particles,
					       num_particles) == created_M); 
		     created_index++)
		{
			short state[6];
			memcpy(state,
			       annihilated_states+annihilated_index*num_particles,
			       num_particles*sizeof(short));
			memcpy(state+num_particles,
			       created_states+created_index*num_particles,
			       num_particles*sizeof(short));
			append_array_element(connections_builder,
					     state);
		}
	}
	free_array_builder(connections_builder);
	free(annihilated_states);
	free(created_states);
	return single_particle_connections;
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

static
short *load_states(const char *directory,int E,size_t *length)
{
	char *filename = calloc(strlen(directory)+256,sizeof(char));
	sprintf(filename,"%s/basis_energy_%d",directory,E);
	FILE *file = fopen(filename,"r");
	if (file == NULL)
	{
		fprintf(stderr,"Could not open file \"%s\". %s\n",
			filename,
			strerror(errno));
		exit(1);
	}
	fseek(file,0,SEEK_END);
	*length = ftell(file)/sizeof(short);
	fseek(file,0,SEEK_SET);
	short *states = (short*)malloc((*length)*sizeof(short));
	if (fread(states,sizeof(short),*length,file) != (*length))
	{
		fprintf(stderr,"Could not read %lu bytes from file \"%s\"\n",
			(*length)*sizeof(short),
			filename);	
		exit(1);
	}
	fclose(file);
	free(filename);
	return states;
}

	static
int single_species_compute_M(single_particle_basis_t basis,
			     short *state,
			     size_t num_particles)
{
	int M = 0;
	for (size_t particle_index = 0; 
	     particle_index < num_particles; 
	     particle_index++)
	{
		M += get_m(get_state(basis,2*(state[particle_index])));
	}
	return M;
}
