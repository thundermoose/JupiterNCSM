#include <transformed_block_manager/transformed_block_manager.h>
#include <bases/m_scheme_2p_basis.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <block_transform/block_transform.h>
#include <string_tools/string_tools.h>

typedef struct
{
	m_scheme_2p_basis_t ket_basis;
	m_scheme_2p_basis_t bra_basis;
	Dens_Matrix *matrix;
} block_t;

struct _transformed_block_manager_
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
void setup_ket_basis(transformed_block_manager_t manager,
		     transform_block_settings_t block_settings);

static
void setup_bra_basis(transformed_block_manager_t manager,
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
			   m_scheme_2p_basis_t basis);

static inline
const size_t get_bra_index(connection_t connection,
			   m_scheme_2p_basis_t basis);

static
void expand_block_list(transformed_block_manager_t manager,
		       size_t num_blocks);

static
void free_blocks(transformed_block_manager_t manager);

transformed_block_manager_t 
new_transformed_block_manager(antoine_2nf_file_t coupled_2nf_data,
			      const char *basis_files_path,
			      const char *index_list_path,
			      int single_particle_energy_max)
{
	transformed_block_manager_t manager =
		(transformed_block_manager_t)
		malloc(sizeof(struct _transformed_block_manager_));
	manager->ket_basis = NULL;
	manager->bra_basis = NULL;
	manager->matrices = NULL;
	manager->coupled_2nf_data = coupled_2nf_data;
	manager->basis_files_path = copy_string(basis_files_path);
	manager->index_list_path = copy_string(index_list_path);
	manager->single_particle_energy_max = single_particle_energy_max;
	manager->single_particle_basis =
	       	new_single_particle_basis(single_particle_energy_max);
	manager->J_max = 2*(single_particle_energy_max*2+1);
	manager->clebsch_gordan_data = initiate_clebsch_gordan(manager->J_max);
	return manager;
}

void decouple_transform_block(transformed_block_manager_t manager,
			      transform_block_settings_t block_settings)
{
	if (manager->ket_basis != NULL)
		free_m_scheme_2p_basis(manager->ket_basis);
	if (manager->bra_basis != NULL)
		free_m_scheme_2p_basis(manager->bra_basis);
	setup_ket_basis(manager,block_settings);
	setup_bra_basis(manager,block_settings);
	int min_M = max(get_min_M(manager->ket_basis),
			get_min_M(manager->bra_basis));
	int max_M = min(get_max_M(manager->ket_basis),
			get_max_M(manager->bra_basis));
	int Tz = block_settings.total_isospin;
	size_t num_blocks = (max_M-min_M)/2+1;
	manager->min_M = min_M;
	if (num_blocks > manager->num_allocated_blocks)
		expand_block_list(manager,num_blocks);
	for (size_t i = 0; i < num_blocks; i++)
	{	
		int M = min_M+i*2;
		m_scheme_2p_basis_t ket_m_basis = 
			cut_out_M_block(manager->ket_basis,M);
		m_scheme_2p_basis_t bra_m_basis =
			cut_out_M_block(manager->bra_basis,M);
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
		if (manager->blocks[i].ket_basis)
			free_m_scheme_2p_basis(manager->blocks[i].ket_basis);
		manager->blocks[i].ket_basis = ket_m_basis;
		if (manager->blocks[i].bra_basis)
			free_m_scheme_2p_basis(manager->blocks[i].bra_basis);
		manager->blocks[i].bra_basis = bra_m_basis;
	}

}

mercury_matrix_block_t 
get_transformed_matrix_block(transformed_block_manager_t manager,
			     matrix_block_setting_t settings)
{
	connection_list_t connection_list = 
		read_connection_files(manager->index_list_path,
				      settings);
	size_t num_elements = num_connections(connection_list);
	double *elements = (double*)calloc(num_elements,sizeof(double));
	size_t element_index = 0;
	while (has_next_connection(connection_list))
	{
		connection_t current_connection =
			next_connection(connection_list);
		int M = compute_M(current_connection,
				  manager->single_particle_basis);
		size_t i = (M-manager->min_M)/2;
		m_scheme_2p_basis_t ket_m_basis = manager->blocks[i].ket_basis;
		m_scheme_2p_basis_t bra_m_basis = manager->blocks[i].bra_basis;
		Dens_Matrix *current_matrix = manager->blocks[i].matrix; 
		size_t ket_index = get_ket_index(current_connection,
						 ket_m_basis);
		size_t bra_index = get_bra_index(current_connection,
						 bra_m_basis);
		elements[element_index++] =
		       	get_dens_matrix_element(current_matrix,
						ket_index,
						bra_index);
	}
	free_connection_list(connection_list);
	return new_mercury_matrix_block_from_data(elements,
						  num_elements,
						  settings);
}

void free_transformed_block_manager(transformed_block_manager_t manager)
{
	free_blocks(manager);
}

static
void setup_ket_basis(transformed_block_manager_t manager,
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
void setup_bra_basis(transformed_block_manager_t manager,
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
	size_t length_filename_buffer = strlen(basis_files_path)+7;
	char *proton_basis_filename = NULL;
	char *neutron_basis_filename = NULL;
	if (total_isospin < 2)
	{
		proton_basis_filename = (char*)malloc(length_filename_buffer);
		sprintf(proton_basis_filename,
			"%s/%s_inds_index_lists/basis_energy_%d",
			proton_basis_files_path,
			total_isospin == 0 ? "p" : "pp", 
			proton_energy);			
	}
	if (total_isospin > -2)
	{
		neutron_basis_filename = (char*)malloc(length_filename_buffer);
		sprintf(neutron_basis_filename,
			"%s/%s_inds_index_lists/basis_energy_%d",
			neutron_basis_files_path,
			total_isospin == 0 ? "n" : "nn",
			neutron_energy);
	}
	m_scheme_2p_basis_t basis =
	       	new_m_scheme_2p_basis_from_file(single_particle_energy_max,
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
	return sp_states->sp_states[bottom_state.a].m + 
		sp_states->sp_states[bottom_state.b].m;
}

static inline
const int get_max_M(m_scheme_2p_basis_t basis)
{
	m_scheme_2p_state_t top_state = 
		get_m_scheme_2p_state(basis,get_m_scheme_2p_dimension(basis));
	SP_States *sp_states = get_m_scheme_sp_states(basis);
	return sp_states->sp_states[top_state.a].m + 
		sp_states->sp_states[top_state.b].m;
}

static
void expand_block_list(transformed_block_manager_t manager,
			size_t num_blocks)
{
	if (manager->blocks != NULL)
		free_blocks(manager);
	manager->num_allocated_blocks = num_blocks;
	manager->blocks = (block_t*)calloc(num_blocks,sizeof(block_t));
}

static
void free_blocks(transformed_block_manager_t manager)
{

	for (size_t i = 0; i <manager->num_blocks; i++)
	{
		if (manager->blocks[i].matrix)
			free_dens_matrix(manager->blocks[i].matrix);
		if (manager->blocks[i].ket_basis)
			free_m_scheme_2p_basis(manager->blocks[i].bra_basis);
		if (manager->blocks[i].bra_basis)
			free_m_scheme_2p_basis(manager->blocks[i].ket_basis);
	}
	free(manager->blocks);
}

static inline
const size_t get_ket_index(connection_t connection,
			   m_scheme_2p_basis_t basis)
{
	size_t num_protons = count_protons(connection.type);
	if (num_protons == 2)
	{
		m_scheme_2p_basis_t state =
		{
			.a = 2*connection.proton_states[0],
			.b = 2*connection.proton_states[1];
		};
		if (state.a > state.b)
			swap(&state.a,&state.b);
		return state;
	}
	else if (num_protons == 1)
	{

	}	
	else
	{

	}
}

static inline
const size_t get_bra_index(connection_t connection,
			   m_scheme_2p_basis_t basis)
{
}

static inline
void swap(int *a,int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}
