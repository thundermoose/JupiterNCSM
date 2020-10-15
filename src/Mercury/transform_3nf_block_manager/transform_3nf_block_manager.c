#include <transform_3nf_block_manager/transform_3nf_block_manager.h>
#include <bases/m_scheme_3p_basis.h>
#include <single_particle_basis/single_particle_basis.h>
#include <string_tools/string_tools.h>
#include <clebsch_gordan/clebsch_gordan.h>
#include <matrix_transform/matrix_transform.h>
#include <block_transform/block_transform.h>
#include <transformed_3nf_M_block/transformed_3nf_M_block.h>
#include <string.h>
#include <error/error.h>
#include <time.h>

struct _tranform_2nf_block_manager_
{
	M_Scheme_3p_Basis *ket_basis;
	M_Scheme_3p_Basis *bra_basis;
	Data_File *coupled_3nf_data;
	char *index_list_path;
	transformed_3nf_M_block_t *blocks;
	size_t num_blocks;
	size_t num_allocated_blocks;
	int min_M;
	int single_particle_energy_max;
	int J_max;
	single_particle_basis_t single_particle_basis;
	Clebsch_Gordan_Data *clebsch_gordan_data;
};

static
void setup_ket_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block);
static
void setup_bra_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block);

static
M_Scheme_3p_Basis *new_ket_basis(transform_3nf_block_manager_t manager,
				matrix_energy_block_t block);
static
M_Scheme_3p_Basis *new_bra_basis(transform_3nf_block_manager_t manager,
				matrix_energy_block_t block);

static
int get_min_M(M_Scheme_3p_Basis *basis);

static
int get_max_M(M_Scheme_3p_Basis *basis);

static
void clear_blocks(transform_3nf_block_manager_t manager);

static
void expand_block_list(transform_3nf_block_manager_t manager);

static
size_t get_ket_index(connection_t connection,
		     M_Scheme_3p_Basis *basis,
		     int *phase);

static
size_t get_bra_index(connection_t connection,
		     M_Scheme_3p_Basis *basis,
		     int *phase);

static
M_Scheme_3p_State set_up_ket_state(connection_t connection);

static
M_Scheme_3p_State set_up_bra_state(connection_t connection);

static
void sort_state(M_Scheme_3p_State *state, int *phase);

static inline
void swap(int *a,int *b);

static
M_Scheme_3p_Basis *load_anicr_3p_basis(const char *index_list_path,
				       int total_isospin,
				       int proton_energy,
				       int neutron_energy,
				       SP_States *sp_basis);

static inline
char *setup_basis_filename(const char *index_list_path,
			   int total_isospin,
			   int energy,
			   const char *particle_map[4]);

static inline
int max(int a, int b);

static inline
int min(int a, int b);

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
#pragma omp parallel
	{
		manager->clebsch_gordan_data =
		       	initiate_clebsch_gordan(manager->J_max);
	}
	return manager;
}

void decouple_transform_3nf_block(transform_3nf_block_manager_t manager,
				  transform_block_settings_t block)
{
	struct timespec time_start;
	clock_gettime(CLOCK_REALTIME,&time_start);
	setup_ket_basis(manager,block);
	setup_bra_basis(manager,block);
	int min_M = max(get_min_M(manager->ket_basis),
			get_min_M(manager->bra_basis));
	int max_M = min(get_max_M(manager->ket_basis),
			get_max_M(manager->bra_basis));
	int bra_energy = block.proton_energy_bra + block.neutron_energy_bra;
	int ket_energy = block.proton_energy_ket + block.neutron_energy_ket;
	int total_isospin = block.total_isospin;
	manager->min_M = min_M;
	manager->num_blocks = (max_M-min_M)/2 + 1;
	clear_blocks(manager);
	if (manager->num_blocks > manager->num_allocated_blocks)
		expand_block_list(manager);
	for (size_t i = 0; i < manager->num_blocks; i++)
	{
		int M = min_M + i*2;
		size_t ket_offset = 0;
		M_Scheme_3p_Basis *ket_m_basis =
			generate_block(manager->ket_basis,
				       total_isospin,
				       M,
				       ket_energy,
				       &ket_offset);
		size_t bra_offset = 0;	
		M_Scheme_3p_Basis *bra_m_basis =
			generate_block(manager->bra_basis,
				       total_isospin,
				       M,
				       bra_energy,
				       &bra_offset);
		if (manager->blocks[i].matrix)
			free_dens_matrix(manager->blocks[i].matrix);
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
	struct timespec time_end;
	clock_gettime(CLOCK_REALTIME,&time_end);
	double time_difference = 
		(time_end.tv_sec-time_start.tv_sec)*1e6+
		(time_end.tv_nsec-time_start.tv_nsec)*1e-3;
	printf("decouple_transfrom_2nf_block: %lg µs\n",
	       time_difference);
}

	mercury_matrix_block_t
get_transform_3nf_matrix_block(transform_3nf_block_manager_t manager,
			       matrix_block_setting_t settings)
{
	struct timespec time_start;
	clock_gettime(CLOCK_REALTIME,&time_start);
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
		M_Scheme_3p_Basis *ket_m_basis = manager->blocks[i].ket_basis;
		M_Scheme_3p_Basis *bra_m_basis = manager->blocks[i].bra_basis;
		Dens_Matrix *current_matrix = manager->blocks[i].matrix;
		int phase = 1;
		size_t ket_index = get_ket_index(current_connection,
						 ket_m_basis,
						 &phase);
		size_t bra_index = get_bra_index(current_connection,
						 bra_m_basis,
						 &phase);
		elements[element_index++] =
			phase*get_dens_matrix_element(current_matrix,
						      bra_index,
						      ket_index);
	}
	free_connection_list(connection_list);
	struct timespec time_end;
	clock_gettime(CLOCK_REALTIME,&time_end);
	double time_difference = 
		(time_end.tv_sec-time_start.tv_sec)*1e6+
		(time_end.tv_nsec-time_start.tv_nsec)*1e-3;
	printf("get_transform_3nf_matrix_block: %lg µs\n",
	       time_difference);
	return new_mercury_matrix_block_from_data(elements,
						  num_elements,
						  settings);
}

transformed_block_t 
get_transformed_block(transform_3nf_block_manager_t manager,
		      matrix_energy_block_t block)
{
	transformed_block_t transformed_block =
	       	new_empty_transformed_block(block,
					    manager->index_list_path,
					    manager->single_particle_basis);
	M_Scheme_3p_Basis *ket_basis = new_ket_basis(manager,block);
	set_transformed_block_ket_basis(transformed_block,
				       	ket_basis);
	M_Scheme_3p_Basis *bra_basis = new_bra_basis(manager,block);
	set_transformed_block_bra_basis(transformed_block,
				       	bra_basis);
	initialized_M_blocks(transformed_block);
	int min_M = get_transformed_block_min_M(transformed_block);
	int bra_energy = get_bra_energy(block);
	int ket_energy = get_ket_energy(block);
	int total_isospin = get_total_isospin(block);
	size_t num_blocks = get_num_M_blocks(transformed_block);
	for (size_t i = 0; i<num_blocks; i++)
	{
		int M = min_M + i*2;
		size_t ket_offset = 0;
		M_Scheme_3p_Basis *ket_m_basis =
			generate_block(ket_basis,
				       total_isospin,
				       M,
				       ket_energy,
				       &ket_offset);
		size_t bra_offset = 0;	
		M_Scheme_3p_Basis *bra_m_basis =
			generate_block(bra_basis,
				       total_isospin,
				       M,
				       bra_energy,
				       &bra_offset);
		transformed_3nf_M_block_t current_M_block =
		{
			.ket_basis = ket_m_basis,
			.bra_basis = bra_m_basis,
			.matrix =
			       	compute_jjj_block(bra_m_basis,
						  ket_m_basis,
						  manager->coupled_3nf_data,
						  manager->clebsch_gordan_data)
				
		};
		set_transformed_block_m_block(transformed_block,
					      current_M_block,i);
	}
	return transformed_block;
}

void free_transform_3nf_block_manager(transform_3nf_block_manager_t manager)
{
	free(manager->index_list_path);
	clear_blocks(manager);
	if (manager->blocks)
		free(manager->blocks);
	if (manager->ket_basis)
		free_m_scheme_3p_basis(manager->ket_basis);
	if (manager->bra_basis)
		free_m_scheme_3p_basis(manager->bra_basis);
	free_single_particle_basis(manager->single_particle_basis);
	free_clebsch_gordan(manager->clebsch_gordan_data);
	free(manager);
}

	static
void setup_ket_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block)
{
	if (manager->ket_basis != NULL)
		free_m_scheme_3p_basis(manager->ket_basis);
	SP_States *sp_basis = get_sp_states(manager->single_particle_basis);
	manager->ket_basis =
		load_anicr_3p_basis(manager->index_list_path,
				    block.total_isospin,
				    block.proton_energy_ket,
				    block.neutron_energy_ket,
				    sp_basis);	 
}
	static
void setup_bra_basis(transform_3nf_block_manager_t manager,
		     transform_block_settings_t block)
{
	if (manager->bra_basis != NULL)
		free_m_scheme_3p_basis(manager->bra_basis);
	SP_States *sp_basis = get_sp_states(manager->single_particle_basis);
	manager->bra_basis =
		load_anicr_3p_basis(manager->index_list_path,
				    block.total_isospin,
				    block.proton_energy_bra,
				    block.neutron_energy_bra,
				    sp_basis);	 
}

M_Scheme_3p_Basis *new_ket_basis(transform_3nf_block_manager_t manager,
		     matrix_energy_block_t block)
{
	SP_States *sp_basis = get_sp_states(manager->single_particle_basis);
	return load_anicr_3p_basis(manager->index_list_path,
				   get_total_isospin(block),
				   get_proton_energy_ket(block),
				   get_neutron_energy_ket(block),
				   sp_basis);	 
}
	static
M_Scheme_3p_Basis *new_bra_basis(transform_3nf_block_manager_t manager,
		     matrix_energy_block_t block)
{
	SP_States *sp_basis = get_sp_states(manager->single_particle_basis);
	return load_anicr_3p_basis(manager->index_list_path,
				   get_total_isospin(block),
				   get_proton_energy_bra(block),
				   get_neutron_energy_bra(block),
				   sp_basis);	 
}

	static
int get_min_M(M_Scheme_3p_Basis *basis)
{
	M_Scheme_3p_State *states = get_m_scheme_3p_states(basis);
	int m = get_3p_M(states[0],basis->sp_states); 
	free(states);
	return m;
}

	static
int get_max_M(M_Scheme_3p_Basis *basis)
{
	M_Scheme_3p_State *states = get_m_scheme_3p_states(basis);
	const size_t last_element = get_m_scheme_3p_dimension(basis) - 1;
	int m = get_3p_M(states[last_element],basis->sp_states); 
	free(states);
	return m;
}

	static
void clear_blocks(transform_3nf_block_manager_t manager)
{
	for (size_t i = 0; i<manager->num_allocated_blocks; i++)
	{
		if (manager->blocks[i].ket_basis)
			free_m_scheme_3p_basis(manager->blocks[i].ket_basis);
		manager->blocks[i].ket_basis = NULL;
		if (manager->blocks[i].bra_basis)
			free_m_scheme_3p_basis(manager->blocks[i].bra_basis);
		manager->blocks[i].bra_basis = NULL;
		if (manager->blocks[i].matrix)
			free_dens_matrix(manager->blocks[i].matrix);
		manager->blocks[i].matrix = NULL;
	}
}

	static
void expand_block_list(transform_3nf_block_manager_t manager)
{
	manager->num_allocated_blocks = manager->num_blocks;
	free(manager->blocks);
	manager->blocks = 
		(transformed_3nf_M_block_t*)
		calloc(manager->num_allocated_blocks,sizeof(transformed_3nf_M_block_t));
}

	static
size_t get_ket_index(connection_t connection,
		     M_Scheme_3p_Basis *basis,
		     int *phase)
{
	M_Scheme_3p_State state = set_up_ket_state(connection);
	sort_state(&state,phase);
	return find_m_scheme_3p_state(basis,state);
}

	static
size_t get_bra_index(connection_t connection,
		     M_Scheme_3p_Basis *basis,
		     int *phase)
{
	M_Scheme_3p_State state = set_up_bra_state(connection);
	sort_state(&state,phase);
	return find_m_scheme_3p_state(basis,state);
}

	static
M_Scheme_3p_State set_up_ket_state(connection_t connection)
{
	switch (count_protons(connection.type))
	{
		case 0:
			return new_m_scheme_3p_state
				(connection.neutron_states[0]*2+1,
				 connection.neutron_states[1]*2+1,
				 connection.neutron_states[2]*2+1);

		case 1:
			return new_m_scheme_3p_state
				(connection.neutron_states[0]*2+1,
				 connection.neutron_states[1]*2+1,
				 connection.proton_states[0]*2);
		case 2:
			return new_m_scheme_3p_state
				(connection.neutron_states[0]*2+1,
				 connection.proton_states[0]*2,
				 connection.proton_states[1]*2);
		case 3:
			return new_m_scheme_3p_state
				(connection.proton_states[0]*2,
				 connection.proton_states[1]*2,
				 connection.proton_states[2]*2);
		default:
			error("There are two many protons,"
			      " in the connection\n");
	}
}

	static
M_Scheme_3p_State set_up_bra_state(connection_t connection)
{
	switch (count_protons(connection.type))
	{
		case 0:
			return new_m_scheme_3p_state
				(connection.neutron_states[3]*2+1,
				 connection.neutron_states[4]*2+1,
				 connection.neutron_states[5]*2+1);

		case 1:
			return new_m_scheme_3p_state
				(connection.neutron_states[2]*2+1,
				 connection.neutron_states[3]*2+1,
				 connection.proton_states[1]*2);
		case 2:
			return new_m_scheme_3p_state
				(connection.neutron_states[1]*2+1,
				 connection.proton_states[2]*2,
				 connection.proton_states[3]*2);
		case 3:
			return new_m_scheme_3p_state
				(connection.proton_states[3]*2,
				 connection.proton_states[4]*2,
				 connection.proton_states[5]*2);
		default:
			error("There are two many protons,"
			      " in the connection\n");
	}
}

	static
void sort_state(M_Scheme_3p_State *state, int *phase)
{
	if (state->a > state->b) 
	{
		*phase = -*phase;
		swap(&state->a,&state->b);
	}
	if (state->b > state->c) 
	{
		*phase = -*phase;
		swap(&state->b,&state->c);
	}
	if (state->a > state->b) 
	{
		*phase = -*phase;
		swap(&state->a,&state->b);
	}
}

	static inline
void swap(int *a,int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

	static
M_Scheme_3p_Basis *load_anicr_3p_basis(const char *index_list_path,
				       int total_isospin,
				       int proton_energy,
				       int neutron_energy,
				       SP_States *sp_basis)
{
	static const char *proton_map[4] = {"ppp","pp","p",""};
	static const char *neutron_map[4] = {"","n","nn","nnn"};
	char *proton_basis_filename = 
		setup_basis_filename(index_list_path,
				     total_isospin,
				     proton_energy,
				     proton_map);
	char *neutron_basis_filename = 
		setup_basis_filename(index_list_path,
				     total_isospin,
				     neutron_energy,
				     neutron_map);
	const size_t num_neutrons = (total_isospin+3)/2;
	const size_t num_protons = 3-num_neutrons;
	M_Scheme_3p_Basis *basis =
		new_m_scheme_3p_basis_from_basis_file(proton_basis_filename,
						      neutron_basis_filename,
						      num_protons,
						      num_neutrons,
						      sp_basis);
	free(proton_basis_filename);
	free(neutron_basis_filename);
	return basis;
}

	static inline
char *setup_basis_filename(const char *index_list_path,
			   int total_isospin,
			   int energy,
			   const char *particle_map[4])
{
	const char *particle_pattern = particle_map[(total_isospin+3)/2];
	size_t length_basis_filename = strlen(index_list_path) + 
		strlen(particle_pattern) + 128;
	char *basis_filename = (char*)malloc(length_basis_filename);
	sprintf(basis_filename,
		"%s/%s_inds_index_lists/basis_energy_%d",
		index_list_path,
		particle_pattern,
		energy);
	return basis_filename;
}

static inline
int max(int a, int b)
{
	return a > b ? a : b;
}

static inline
int min(int a, int b)
{
	return a < b ? a : b;
}
