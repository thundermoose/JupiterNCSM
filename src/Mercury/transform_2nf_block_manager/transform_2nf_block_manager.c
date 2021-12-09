#include <transform_2nf_block_manager/transform_2nf_block_manager.h>
#include <bases/m_scheme_2p_basis.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <block_transform/block_transform.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <debug_mode/debug_mode.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct
{
	m_scheme_2p_basis_t ket_basis;
	m_scheme_2p_basis_t bra_basis;
	Dens_Matrix *matrix;
} block_t;

struct _transform_2nf_block_manager_
{
	m_scheme_2p_basis_t ket_basis;
	m_scheme_2p_basis_t bra_basis;	
	antoine_2nf_file_t coupled_2nf_data;
	char *basis_files_path;
	char *index_list_path;
	block_t *blocks;
	size_t num_blocks;
	size_t num_allocated_blocks;
	int min_M;
	int single_particle_energy_max;
	int J_max;
	single_particle_basis_t single_particle_basis;
	Clebsch_Gordan_Data* clebsch_gordan_data;
};

static
void setup_ket_basis(transform_2nf_block_manager_t manager,
		     transform_block_settings_t block_settings);

static
void setup_bra_basis(transform_2nf_block_manager_t manager,
		     transform_block_settings_t block_settings);


static
m_scheme_2p_basis_t load_anicr_basis(char *basis_files_path,
				     int total_isospin,
				     int proton_energy,
				     int neutron_energy,
				     int single_particle_energy_max);

static inline
const int min(const int a, const int b);

static inline
const int max(const int a, const int b);

static inline
const int get_min_M(m_scheme_2p_basis_t basis);

static inline
const int get_max_M(m_scheme_2p_basis_t basis);

static inline
const size_t get_ket_index(connection_t connection,
			   m_scheme_2p_basis_t basis,
			   int *phase);

static inline
const size_t get_bra_index(connection_t connection,
			   m_scheme_2p_basis_t basis,
			   int *phase);

static
void expand_block_list(transform_2nf_block_manager_t manager,
		       size_t num_blocks);

static
void free_blocks(transform_2nf_block_manager_t manager);

static inline
void swap(int *a,int *b);

#ifdef DEBUG
static inline
void log_matrix(Dens_Matrix *matrix);
#endif

transform_2nf_block_manager_t 
new_transform_2nf_block_manager(antoine_2nf_file_t coupled_2nf_data,
			      const char *basis_files_path,
			      const char *index_list_path,
			      int single_particle_energy_max)
{
	transform_2nf_block_manager_t manager =
		(transform_2nf_block_manager_t)
		calloc(1,sizeof(struct _transform_2nf_block_manager_));
	manager->coupled_2nf_data = coupled_2nf_data;
	manager->basis_files_path = copy_string(basis_files_path);
	manager->index_list_path = copy_string(index_list_path);
	manager->single_particle_energy_max = single_particle_energy_max;
	manager->single_particle_basis =
	       	new_single_particle_basis(single_particle_energy_max);
	manager->J_max = 2*(single_particle_energy_max*2+1);
#pragma omp parallel
	{
		manager->clebsch_gordan_data = 
			initiate_clebsch_gordan(manager->J_max);
	}
	return manager;
}

void decouple_transform_2nf_block(transform_2nf_block_manager_t manager,
				  transform_block_settings_t block_settings)
{
	log_entry("decouple_transform_block(%p,{pk%d pb%d, nk%d nb%d, Tz%d})",
		  manager,
		  block_settings.proton_energy_ket,
		  block_settings.proton_energy_bra,
		  block_settings.neutron_energy_ket,
		  block_settings.neutron_energy_bra,
		  block_settings.total_isospin);
	setup_ket_basis(manager,block_settings);
	setup_bra_basis(manager,block_settings);
	log_entry("Set up the new bases");
	int min_M = max(get_min_M(manager->ket_basis),
			get_min_M(manager->bra_basis));
	int max_M = min(get_max_M(manager->ket_basis),
			get_max_M(manager->bra_basis));
	log_entry("min_M = %d", min_M);
	log_entry("max_M = %d", max_M);
	int Tz = block_settings.total_isospin;
	size_t num_blocks = (max_M-min_M)/2+1;
	manager->min_M = min_M;
	if (num_blocks > manager->num_allocated_blocks)
		expand_block_list(manager,num_blocks);
	for (size_t i = 0; i < num_blocks; i++)
	{	
		int M = min_M+i*2;
		log_entry("Decoupling M = %d",M);
		m_scheme_2p_basis_t ket_m_basis = 
			cut_out_M_block(manager->ket_basis,M);
		m_scheme_2p_basis_t bra_m_basis =
			cut_out_M_block(manager->bra_basis,M);
		log_m_scheme_2p_basis(ket_m_basis);
		log_m_scheme_2p_basis(bra_m_basis);
		if (manager->blocks[i].matrix)
			free_dens_matrix(manager->blocks[i].matrix);
		manager->blocks[i].matrix =
			compute_jt_block(bra_m_basis,
					 ket_m_basis,
					 Tz,
					 M,
					 manager->J_max,
					 manager->coupled_2nf_data,
					 manager->clebsch_gordan_data);
#ifdef DEBUG
		log_matrix(manager->blocks[i].matrix);
#endif
		if (manager->blocks[i].ket_basis)
			free_m_scheme_2p_basis(manager->blocks[i].ket_basis);
		manager->blocks[i].ket_basis = ket_m_basis;
		if (manager->blocks[i].bra_basis)
			free_m_scheme_2p_basis(manager->blocks[i].bra_basis);
		manager->blocks[i].bra_basis = bra_m_basis;
	}

}

mercury_matrix_block_t 
get_transform_2nf_matrix_block(transform_2nf_block_manager_t manager,
			       matrix_block_setting_t settings)
{
	log_entry("get_transformed_matrix_block(%p, {\n"
		  "\t.type = %s\n"
		  "\t.difference_energy_protons = %d\n"
		  "\t.difference_M_protons = %d\n"
		  "\t.depth_protons = %d\n"
		  "\t.difference_energy_neutrons = %d\n"
		  "\t.difference_M_neutrons = %d\n"
		  "\t.depth_neutrons = %d\n"
		  "\t.num_proton_combinations = %lu\n"
		  "\t.num_neutron_combinations = %lu\n"
		  "\t.matrix_block_id = %lu\n"
		  "}",
		manager,
		block_type_to_string(settings.type),
		settings.difference_energy_protons,
		settings.difference_M_protons,
		settings.depth_protons,
		settings.difference_energy_neutrons,
		settings.difference_M_neutrons,
		settings.depth_neutrons,
		settings.num_proton_combinations,
		settings.num_neutron_combinations,
		settings.matrix_block_id);
	connection_list_t connection_list = 
		new_connection_list(manager->index_list_path,
				    manager->single_particle_basis,
				    settings);
	size_t num_elements = num_connections(connection_list);
	double *elements = (double*)calloc(num_elements,sizeof(double));
	size_t element_index = 0;
	log_entry("Creating matrix block %lu", settings.matrix_block_id);
	while (has_next_connection(connection_list))
	{
		connection_t current_connection =
			next_connection(connection_list);
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
				  manager->single_particle_basis);
		log_entry("M = %d",M);
		size_t i = (M-manager->min_M)/2;
		m_scheme_2p_basis_t ket_m_basis = manager->blocks[i].ket_basis;
		m_scheme_2p_basis_t bra_m_basis = manager->blocks[i].bra_basis;
		log_m_scheme_2p_basis(ket_m_basis);
		log_m_scheme_2p_basis(bra_m_basis);
		int phase = 1;
		size_t ket_index = get_ket_index(current_connection,
						 ket_m_basis,
						 &phase);
		size_t bra_index = get_bra_index(current_connection,
						 bra_m_basis,
						 &phase);
		log_entry("ket_index = %lu",ket_index);
		log_entry("bra_index = %lu",bra_index);
		Dens_Matrix *current_matrix = manager->blocks[i].matrix; 
#ifndef NLOGING
		log_entry("current_matrix:");
		for (size_t k = 0; k<current_matrix->m; k++)
		{
			for (size_t l = 0; l<current_matrix->n; l++)
			{
				log_entry("(%lu,%lu): %lg",
					  k,l,
					  get_dens_matrix_element(current_matrix,
								  k,l));
			}
		}
#endif
		elements[element_index++] = 
			phase*get_dens_matrix_element(current_matrix,
						      bra_index,
						      ket_index);
		log_entry("Current element(%lu): %lg",
			  element_index-1,
			  elements[element_index-1]);
	}
	free_connection_list(connection_list);
	return new_mercury_matrix_block_from_data(elements,
						  num_elements,
						  settings);
}

void free_transform_2nf_block_manager(transform_2nf_block_manager_t manager)
{
	free_blocks(manager);
	if (manager->ket_basis)
		free_m_scheme_2p_basis(manager->ket_basis);	
	if (manager->bra_basis)
		free_m_scheme_2p_basis(manager->bra_basis);	
	free(manager->basis_files_path);
	free(manager->index_list_path);
	free_single_particle_basis(manager->single_particle_basis);
	free_clebsch_gordan(manager->clebsch_gordan_data);
	free(manager);
}

static
void setup_ket_basis(transform_2nf_block_manager_t manager,
		     transform_block_settings_t block_settings)
{
	if (manager->ket_basis != NULL)
		free_m_scheme_2p_basis(manager->ket_basis);
	manager->ket_basis =
		load_anicr_basis(manager->basis_files_path,
				 block_settings.total_isospin,
				 block_settings.proton_energy_ket,
				 block_settings.neutron_energy_ket,
				 manager->single_particle_energy_max);
}

static
void setup_bra_basis(transform_2nf_block_manager_t manager,
		     transform_block_settings_t block_settings)
{
	if (manager->bra_basis != NULL)
		free_m_scheme_2p_basis(manager->bra_basis);
	manager->bra_basis =
		load_anicr_basis(manager->basis_files_path,
				 block_settings.total_isospin,
				 block_settings.proton_energy_bra,
				 block_settings.neutron_energy_bra,
				 manager->single_particle_energy_max);
}

static
m_scheme_2p_basis_t load_anicr_basis(char *basis_files_path,
				     int total_isospin,
				     int proton_energy,
				     int neutron_energy,
				     int single_particle_energy_max)
{
	size_t length_filename_buffer = strlen(basis_files_path)+128;
	char *proton_basis_filename = NULL;
	char *neutron_basis_filename = NULL;
	if (total_isospin < 2)
	{
		proton_basis_filename = (char*)malloc(length_filename_buffer);
		sprintf(proton_basis_filename,
			"%s/%s_inds_index_lists/basis_energy_%d",
			basis_files_path,
			total_isospin == 0 ? "p" : "pp", 
			proton_energy);			
	}
	if (total_isospin > -2)
	{
		neutron_basis_filename = (char*)malloc(length_filename_buffer);
		sprintf(neutron_basis_filename,
			"%s/%s_inds_index_lists/basis_energy_%d",
			basis_files_path,
			total_isospin == 0 ? "n" : "nn",
			neutron_energy);
	}
	m_scheme_2p_basis_t basis =
		new_m_scheme_2p_basis_from_files(single_particle_energy_max,
						 proton_basis_filename,
						 neutron_basis_filename);
	free(proton_basis_filename);
	free(neutron_basis_filename);
	return basis;
}

static inline
const int min(const int a, const int b)
{
	return a < b ? a : b;
}

static inline
const int max(const int a, const int b)
{
	return a > b ? a : b;
}

static inline
const int get_min_M(m_scheme_2p_basis_t basis)
{
	m_scheme_2p_state_t bottom_state = get_m_scheme_2p_state(basis,0);
	SP_States *sp_states = get_m_scheme_sp_states(basis);
	log_entry("bottom_state = %d %d", bottom_state.a, bottom_state.b);
	return sp_states->sp_states[bottom_state.a].m + 
		sp_states->sp_states[bottom_state.b].m;
}

static inline
const int get_max_M(m_scheme_2p_basis_t basis)
{
	m_scheme_2p_state_t top_state = 
		get_m_scheme_2p_state(basis,get_m_scheme_2p_dimension(basis)-1);
	SP_States *sp_states = get_m_scheme_sp_states(basis);
	log_entry("top_state = %d %d", top_state.a, top_state.b);
	return sp_states->sp_states[top_state.a].m + 
		sp_states->sp_states[top_state.b].m;
}

static
void expand_block_list(transform_2nf_block_manager_t manager,
			size_t num_blocks)
{
	if (manager->blocks != NULL)
		free_blocks(manager);
	manager->num_allocated_blocks = num_blocks;
	manager->blocks = (block_t*)calloc(num_blocks,sizeof(block_t));
}

static
void free_blocks(transform_2nf_block_manager_t manager)
{

	for (size_t i = 0; i <manager->num_allocated_blocks; i++)
	{
		if (manager->blocks[i].matrix)
			free_dens_matrix(manager->blocks[i].matrix);
		manager->blocks[i].matrix = NULL;
		if (manager->blocks[i].ket_basis)
			free_m_scheme_2p_basis(manager->blocks[i].ket_basis);
		manager->blocks[i].ket_basis = NULL;
		if (manager->blocks[i].bra_basis)
			free_m_scheme_2p_basis(manager->blocks[i].bra_basis);
		manager->blocks[i].bra_basis = NULL;
	}
	free(manager->blocks);
}

static inline
const size_t get_ket_index(connection_t connection,
			   m_scheme_2p_basis_t basis,
			   int *phase)
{
	size_t num_protons = count_protons(connection.type);
	if (num_protons == 2)
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.proton_states[0],
			.b = 2*connection.proton_states[1]
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("ket_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}
	else if (num_protons == 1)
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.neutron_states[0]+1,
			.b = 2*connection.proton_states[0]
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("ket_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}	
	else
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.neutron_states[0]+1,
			.b = 2*connection.neutron_states[1]+1
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("ket_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}
}

static inline
const size_t get_bra_index(connection_t connection,
			   m_scheme_2p_basis_t basis,
			   int *phase)
{
	size_t num_protons = count_protons(connection.type);
	if (num_protons == 2)
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.proton_states[2],
			.b = 2*connection.proton_states[3]
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("bra_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}
	else if (num_protons == 1)
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.neutron_states[1]+1,
			.b = 2*connection.proton_states[1]
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("bra_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}	
	else
	{
		m_scheme_2p_state_t state =
		{
			.a = 2*connection.neutron_states[2]+1,
			.b = 2*connection.neutron_states[3]+1
		};
		if (state.a > state.b)
		{
			*phase = -*phase;
			swap(&state.a,&state.b);
		}
		log_entry("bra_state = %d %d",state.a,state.b);
		return get_m_scheme_2p_state_index(basis,state);
	}
}

static inline
void swap(int *a,int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

#ifdef DEBUG
static inline
void log_matrix(Dens_Matrix *matrix)
{
	log_entry("matrix (%p):",matrix);
	for (size_t i = 0; i<matrix->m; i++)
	{
		for (size_t j = 0; j<matrix->n; j++)
		{
			log_entry("(%lu %lu) = %g",
				  i,j,
				  matrix->elements[matrix->n*i+j]);
		}
	}
}
#endif
