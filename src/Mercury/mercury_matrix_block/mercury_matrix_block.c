#include <mercury_matrix_block/mercury_matrix_block.h>
#include <error/error.h>
#include <global_constants/global_constants.h>
#include <log/log.h>
#include <debug_mode/debug_mode.h>
#include <thundertester/test.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

struct _mercury_matrix_block_
{
	matrix_block_setting_t settings;
	size_t num_elements;
	double *elements;
};

static
double get_matrix_element_from_connection(interaction_t interaction,
					  connection_t connection,
					  single_particle_basis_t basis,
					  int E1, int E2, int M);

mercury_matrix_block_t new_mercury_matrix_block(interaction_t interaction,
						connection_list_t connections,
						single_particle_basis_t basis)
{
	size_t num_elements = num_connections(connections);
	double *elements = (double*)calloc(num_elements,
					   sizeof(double));
	size_t element_index = 0;
	matrix_block_setting_t settings = 
		get_matrix_block_setting(connections);
	int E1 = settings.depth_protons+settings.depth_neutrons;
	int E2 = E1+
		settings.difference_energy_protons+
		settings.difference_energy_neutrons;
	log_entry("Creating matrix block %lu", settings.matrix_block_id);
	while (has_next_connection(connections))
	{
		connection_t current_connection =
			next_connection(connections);
		log_entry("Current Connection(%lu): %s (%d %d %d %d %d %d),"
			 " (%d %d %d %d %d %d)",
			 element_index,
			 block_type_to_string(current_connection.type),
			 current_connection.neutron_states[0],
			 current_connection.neutron_states[1],
			 current_connection.neutron_states[2],
			 current_connection.neutron_states[3],
			 current_connection.neutron_states[4],
			 current_connection.neutron_states[5],
			 current_connection.proton_states[0],
			 current_connection.proton_states[1],
			 current_connection.proton_states[2],
			 current_connection.proton_states[3],
			 current_connection.proton_states[4],
			 current_connection.proton_states[5]);
		int M = compute_M(current_connection,
				  basis);
		elements[element_index++] =
			get_matrix_element_from_connection(interaction,
							   current_connection,
							   basis,
							   E1,E2,M);
		log_entry("element = %lg",elements[element_index-1]);
	}
	mercury_matrix_block_t matrix_block = 
		(mercury_matrix_block_t)
		malloc(sizeof(struct _mercury_matrix_block_));
	matrix_block->settings = settings;
	matrix_block->num_elements = num_elements;
	matrix_block->elements = elements;
	return matrix_block;
}

mercury_matrix_block_t 
new_mercury_matrix_block_from_data(double *elements,
				   size_t num_elements,
				   matrix_block_setting_t settings)
{
	mercury_matrix_block_t block =
	       	(mercury_matrix_block_t)
		malloc(sizeof(struct _mercury_matrix_block_));
	block->elements = elements;
	block->num_elements = num_elements;
	block->settings = settings;
	return block;
}

void save_mercury_matrix_block(mercury_matrix_block_t matrix_block,
			       const char *output_path)
{
	log_entry("save_mercury_matrix_block");
	log_entry("matrix_block = {\n"
		  "\t.settings = {\n"
		  "\t\t.type = %s,\n"
		  "\t\t.difference_energy_protons = %d,\n"
		  "\t\t.difference_M_protons = %d,\n"
		  "\t\t.depth_protons = %d,\n"
		  "\t\t.difference_energy_neutrons = %d,\n"
		  "\t\t.difference_M_neutrons = %d,\n"
		  "\t\t.depth_neutrons = %d,\n"
		  "\t\t.num_proton_combinations = %lu,\n"
		  "\t\t.num_neutron_combinations = %lu,\n"
		  "\t\t.matrix_block_id = %lu"
		  "},\n"
		  ".num_elements = %lu,\n"
		  ".elements = %p}",
		  block_type_to_string(matrix_block->settings.type),
		  matrix_block->settings.difference_energy_protons,
		  matrix_block->settings.difference_M_protons,
		  matrix_block->settings.depth_protons,
		  matrix_block->settings.difference_energy_neutrons,
		  matrix_block->settings.difference_M_neutrons,
		  matrix_block->settings.depth_neutrons,
		  matrix_block->settings.num_proton_combinations,
		  matrix_block->settings.num_neutron_combinations,
		  matrix_block->settings.matrix_block_id,
		  matrix_block->num_elements,
		  matrix_block->elements);
	log_entry("Save block %lu",matrix_block->settings.matrix_block_id);
#ifndef NLOGING
	for (size_t i = 0; i < matrix_block->num_elements; i++)
		log_entry("Matrix element (%lu): %lg",
			  i,matrix_block->elements[i]);
#endif
	size_t array_index = matrix_block->settings.matrix_block_id;
	size_t length_file_name =
		strlen(output_path)+(size_t)(log(array_index)+1)+18;
	char *file_name = (char*)malloc(length_file_name);
	sprintf(file_name,
		"%s/%lu_matrix_elements",
		output_path,
		array_index);
	log_entry("file_name = %s",file_name);
	FILE *file = fopen(file_name,"w");
	if (file == NULL)
		error("Could not open file %s for writing. %s\n",
		      file_name,
		      strerror(errno));
	size_t neutron_dimension =
		matrix_block->settings.num_neutron_combinations;	
	size_t proton_dimension =
		matrix_block->settings.num_proton_combinations;
	if (fwrite(&neutron_dimension,
		   sizeof(size_t),1,
		   file) != 1)		
		error("Could not write neutron dimension to %s\n",
		      file_name);
	if (fwrite(&proton_dimension,
		   sizeof(size_t),1,
		   file) != 1)		
		error("Could not write proton dimension to %s\n",
		      file_name);
	if (fwrite(matrix_block->elements,
		   sizeof(double),
		   matrix_block->num_elements,
		   file) != matrix_block->num_elements)
		error("Could not write the matrix elements to %s\n",
		      file_name);
	fclose(file);
	free(file_name);
}

void free_mercury_matrix_block(mercury_matrix_block_t matrix_block)
{
	free(matrix_block->elements);
	free(matrix_block);
}

static
double get_matrix_element_from_connection(interaction_t interaction,
					  connection_t connection,
					  single_particle_basis_t basis,
					  int E1, int E2, int M)
{
	size_t num_protons = count_protons(connection.type);
	size_t num_neutrons = count_neutrons(connection.type);
	size_t num_particles = num_protons+num_neutrons;	
	int ket_state[max_interaction_order];
	int bra_state[max_interaction_order];
	for (size_t i = 0; i<num_neutrons; i++)
	{
		single_particle_state_t bra_neutron =
			get_state(basis,
				  2*connection.neutron_states[i]+1);
		single_particle_state_t ket_neutron =
			get_state(basis,
				  2*connection.neutron_states[num_neutrons+i]+1);
		bra_state[i] = bra_neutron.neptune_index;
		ket_state[i] = ket_neutron.neptune_index;
		log_entry("bra_neutron.neptune_index = %lu",
			  bra_neutron.neptune_index);
		log_entry("bra_neutron.tz = %d",bra_neutron.tz);
		log_entry("ket_neutron.neptune_index= %lu",
			  ket_neutron.neptune_index);
		log_entry("ket_neutron.tz = %d",ket_neutron.tz);
	}
	for (size_t i = 0; i<num_protons; i++)
	{
		single_particle_state_t bra_proton =
			get_state(basis,
				  2*connection.proton_states[i]);
		single_particle_state_t ket_proton =
			get_state(basis,
				  2*connection.proton_states[num_protons+i]);
		bra_state[num_neutrons+i] = bra_proton.neptune_index;
		ket_state[num_neutrons+i] = ket_proton.neptune_index;
		log_entry("bra_proton.neptune_index = %lu",
			  bra_proton.neptune_index);
		log_entry("bra_proton.tz = %d",bra_proton.tz);
		log_entry("ket_proton.neptune_index = %lu",
			  ket_proton.neptune_index);
		log_entry("ket_proton.tz = %d",ket_proton.tz);
	}
	int Tz = (int)(num_neutrons)-(int)(num_protons);
	return get_matrix_element(interaction,
				  bra_state,
				  ket_state,
				  num_particles,
				  E1,E2,Tz,M);
}

