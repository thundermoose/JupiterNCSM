#include <transform_3nf_block_manager/transform_3nf_block_manager.h>

typedef struct
{
	M_Scheme_3p_Basis *ket_basis;
	M_Scheme_3p_Basis *bra_basis;
	Dens_Matrix *matrix;
} block_t;

static
void setup_ket_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block);
static
void setup_bra_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block);

static
int get_min_M(M_Scheme_3p_Basis basis);

static
int get_max_M(M_Scheme_3p_Basis basis);

static
void expand_block_list(transform_3nf_block_manager_t manager);

static
size_t get_ket_index(M_Scheme_3p_Basis basis,
		     connection_t connection);

static
size_t get_bra_index(M_Scheme_3p_Basis basis,
		     connection_t connection);

struct _tranform_2nf_block_manager_
{
	M_Scheme_3p_Basis *ket_basis;
	M_Scheme_3p_Basis *bra_basis;
	Data_File *coupled_3nf_data;
	char *index_list_path;
	block_t *blocks;
	size_t num_blocks;
	size_t num_allocated_blocks;
	int min_M;
	int single_particle_energy_max;
	int J_max;
	single_particle_basis_t single_particle_basis;
	Clebsch_Gordan_Data *clebsch_gordan_data;
};

transform_3nf_block_manager_t 
new_transform_3nf_block_manager(Data_File *coupled_3nf_data,
				const char *index_list_path,
				const int single_particle_energy)
{
	transform_3nf_block_manager_t manager =
		(transform_3nf_block_manager_t)
		calloc(1,sizeof(struct _tranform_2nf_block_manager_));
	manager->coupled_3nf_data = coupled_3nf_data;
	manager->index_list_path = copy_string(index_list_path);
	manager->single_particle_basis =
	       	new_single_particle_basis(single_particle_energy);
	manager->J_max = 3*(single_particle_energy*2+1);
	manager->clebsch_gordan_data = initiate_clebsch_gordan(manager->J_max);
	return manager;
}

void decouple_transform_3nf_block(transform_3nf_block_manager_t manager,
				  transform_block_settings_t block)
{
	setup_ket_basis(manager,block);
	setup_bra_basis(manager,block);
	int min_M = max(get_min_M(manager->ket_basis),
			get_min_M(manager->bra_basis));
	int max_M = min(get_max_M(manager->ket_basis),
			get_max_M(manager->bra_basis));
	manager->min_M = min_M;
	manager->num_blocks = (max_M-min_M)/2 + 1;
	if (manager->num_blocks > manager->num_allocated_blocks)
		expand_block_list(manager);
	for (size_t i = 0; i < manager->num_blocks; i++)
	{
		int M = min_M + i*2;
		M_Scheme_3p_Basis ket_m_basis =
			cut_out_M_3p_block(manager->ket_basis,M);
		M_Scheme_3p_Basis bra_m_basis =
			cut_out_M_3p_block(manager->bra_basis,M);
		if (manager->blocks[i].manager)
			free_dens_matrix(manager->blocks[i].manager);
		manager->blocks[i].matrix =
			compute_jjj_block(bra_m_basis,
					  ket_m_basis,
					  manager->coupled_3nf_data,
					  manager->clebsch_gordan_data);
		if (manager->blocks[i].ket_basis)
			free_m_scheme_3p_basis(manager->blocks[i].ket_basis);
		manager->blocks[i].ket_basis = ket_m_basis;
		if (manager->blocks[i].bra_basis)
			free_m_scheme_3p_basis(manager->blocks[i].bra_basis);
		manager->blocks[i].bra_basis = bra_m_basis;
	}
}

mercury_matrix_block_t
get_transform_3nf_matrix_block(transform_3nf_block_manager_t manager,
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
		int M = compute_M(current_connection);
		size_t i = (M-manager->min_M)/2;
		M_Scheme_3p_Basis ket_m_basis = manager->blocks[i].ket_basis;
		M_Scheme_3p_Basis bra_m_basis = manager->blocks[i].bra_basis;
		int phase = 1;
		size_t ket_index = get_ket_index(current_connection,
						 ket_m_basis);
		size_t bra_index = get_bra_index(current_connection,
						 bra_m_basis);
		elements[element_index++] =
			phase*get_dens_matrix_element(current_matrix,
						      bra_index,
						      ket_index)
	}
	free_connection_list(connection_list);
	return new_mercury_matrix_block_from_data(element_index,
						  num_elements,
						  settings);
}

void free_transform_3nf_block_manager(transform_3nf_block_manager_t manager)
{
	free(manager->index_list_path);
	free_blocks(manager);
	if (manager->ket_basis)
		free_m_scheme_3p_basis(manager->ket_basis);
	if (manager->bra_basis)
		free_m_scheme_3p_basis(manager->bra_basis);
	free_single_particle_basis(manager->single_particle_basis);
	free(manager);
}

static
void setup_ket_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block)
{
	if (manager->ket_basis != NULL)
		free_m_scheme_3p_basis(manager->ket_basis);
	manager->ket_basis =
		load_anicr_basis(manager->index_list_path,
			 	 block.total_isospin,
			  	 block.proton_energy_ket,
				 block.neutron_energy_ket,
			         manager->single_particle_energy_max);	 
}
static
void setup_bra_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block)
{
	if (manager->bra_basis != NULL)
		free_m_scheme_3p_basis(manager->bra_basis);
	manager->bra_basis =
		load_anicr_basis(manager->index_list_path,
			 	 block.total_isospin,
			  	 block.proton_energy_bra,
				 block.neutron_energy_bra,
			         manager->single_particle_energy_max);	 
}

static
int get_min_M(M_Scheme_3p_Basis basis)
{
}

static
int get_max_M(M_Scheme_3p_Basis basis)
{
}

static
void expand_block_list(transform_3nf_block_manager_t manager)
{
}

static
size_t get_ket_index(M_Scheme_3p_Basis basis,
		     connection_t connection)
{
}

static
size_t get_bra_index(M_Scheme_3p_Basis basis,
		     connection_t connection)
{
}
