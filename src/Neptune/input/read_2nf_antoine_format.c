#include <input/read_2nf_antoine_format.h>
#include <utils/index_hash.h>
#include <utils/debug_messages.h>
#include <log/log.h>
#include <debug_mode/debug_mode.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#define unused(v) (void)(v)

typedef struct
{
	int num_columns;
	float basis_frequency;
	int num_rows;
} antoine_header_t;

struct _antoine_2nf_file_
{
	jt_basis_t antoine_basis;
	size_t num_elements;
	double *elements;
	index_hash_t used_configurations;
};

static
antoine_header_t read_header(FILE* file);

static
double *read_elements(antoine_header_t header,
		      FILE* file);

static
double *compress_elements(antoine_header_t header,
			  jt_basis_t antoine_basis,
			  size_t num_particles,
			  double *raw_elements);

static
double get_element(antoine_2nf_file_t antoine_2nf_file,
		   size_t bra_index,
		   size_t ket_index,
		   quantum_number Tz);

static
void jump_over_fortran_block(FILE* file);

static 
int read_fortran_block(FILE* file,
		       void *buffer,
		       size_t buffer_size);

static
void close_fortran_block(FILE* file,
			 int block_length);

static
int get_fortran_block_length(FILE *file);

antoine_2nf_file_t open_antoine_2nf_file(const char* file_name,
					 size_t num_particles,
					 quantum_number e_max1,
					 quantum_number e_max2)
{
	FILE* file_handle = fopen(file_name,"r");
	if (file_handle == NULL)
	{
		fprintf(stderr,"Could not open file \"%s\". %s.\n",
			file_name,
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	jt_basis_t antoine_basis = new_antoine_basis(e_max1,e_max2);
	antoine_header_t header = read_header(file_handle);
	double *raw_elements = read_elements(header,file_handle);
	double *elements = 
		compress_elements(header,
				  antoine_basis,
				  num_particles,
				  raw_elements);
	free(raw_elements);
	fclose(file_handle);
	assert(antoine_basis != NULL);
	index_hash_t used_configurations = 
		compute_used_indices(antoine_basis);
	antoine_2nf_file_t data_file = 
		(antoine_2nf_file_t)malloc(sizeof(struct _antoine_2nf_file_));	
	data_file->num_elements = header.num_rows;
	data_file->elements = elements;
	data_file->used_configurations = used_configurations;
	data_file->antoine_basis = antoine_basis;
	return data_file;
}

jt_basis_t get_jt_basis(antoine_2nf_file_t data_file)
{
	return data_file->antoine_basis;
}

Dens_Matrix *get_antoine_matrix(antoine_2nf_file_t data_file,
				jt_basis_t bra_basis,
				jt_basis_t ket_basis,
				quantum_number Tz)
{
	log_entry("Tz = %d\n",Tz);
	size_t *bra_indices = get_sub_basis_indices(data_file->antoine_basis,
						    bra_basis);
	size_t *ket_indices = get_sub_basis_indices(data_file->antoine_basis,
						    ket_basis);
	Dens_Matrix *matrix = 
		new_zero_matrix(get_dimension(bra_basis),
				get_dimension(ket_basis));
	for (size_t i = 0; i<get_dimension(bra_basis)*get_dimension(ket_basis); i++)
	{
		size_t bra_index = bra_indices[i / get_dimension(bra_basis)];
		size_t ket_index = ket_indices[i % get_dimension(bra_basis)];
		matrix->elements[i] =
			get_element(data_file,
				    bra_index,
				    ket_index,
				    Tz);
	}
	free(bra_indices);
	free(ket_indices);
	return matrix;
}

void free_antoine_2nf_file(antoine_2nf_file_t data_file)
{
	log_entry("free_antoine_2nf_file(%p)",data_file);
	free_jt_basis(data_file->antoine_basis);
	free_index_hash(data_file->used_configurations);
	free(data_file->elements);
	free(data_file);
}

	static
antoine_header_t read_header(FILE* file)
{
	long initial_position = ftell(file);
	jump_over_fortran_block(file);
	antoine_header_t header;
	if (read_fortran_block(file,
			       (void*)&header,
			       sizeof(antoine_header_t)))
	{
		fprintf(stderr,"Could not read header\n");
		exit(EXIT_FAILURE);
	}
	fseek(file,initial_position,SEEK_SET);
	return header;
}

	static
double *read_elements(antoine_header_t header,
		      FILE* file)
{
	long initial_position = ftell(file);
	jump_over_fortran_block(file);
	jump_over_fortran_block(file);
	float *elements_float =
		(float*)calloc(header.num_rows*header.num_columns,
			       sizeof(float));
	if (read_fortran_block(file,
			       (void*)elements_float,
			       sizeof(float)*
			       header.num_rows*
			       header.num_columns))
	{
		fprintf(stderr,"Could not read matrix elements.%s\n",
			strerror(errno));
		exit(EXIT_FAILURE); } fseek(file,initial_position,SEEK_SET); 
	double *elements = (double*)calloc(header.num_rows*header.num_columns, 
					   sizeof(double)); 
	for (size_t i = 0; i<header.num_rows*header.num_columns; i++)
	{
		elements[i] = elements_float[i];
	}
	free(elements_float);
	return elements;
}

	static
double *compress_elements(antoine_header_t header,
			  jt_basis_t antoine_basis,
			  size_t num_particles,
			  double *raw_elements)
{
	//const size_t dimension = get_dimension(antoine_basis);
	//size_t *diagonal_indices =
	//       	(size_t*)malloc(dimension*sizeof(size_t));
	//double *harmonic_oscillator_energy = 
	//	(double*)malloc(dimension*sizeof(double));
	//index_hash_t used_indices = compute_used_indices(antoine_basis);
	//Shells *shells = get_jt_basis_shells(antoine_basis);
	//for (size_t i = 0; i < dimension; i++)
	//{
	//	diagonal_indices[i] = get_index(used_indices,i,i);
	//	jt_state_t state = get_jt_state(antoine_basis,
	//					i);
	//	quantum_number e_a = shells->true_shells[state.a].e;
	//	quantum_number e_b = shells->true_shells[state.b].e;
	//	harmonic_oscillator_energy[i] = (e_a+e_b)+3.0/4;
	//	//harmonic_oscillator_energy[i] = 3.0/4;
	//}
	//free_index_hash(used_indices);
	double *elements = (double*)malloc(header.num_rows*3*sizeof(double));
	printf("header.basis_frequency = %lg\n",
	       header.basis_frequency);
	//size_t diagonal_index = 0;
	for (size_t i = 0; i<header.num_rows; i++)
	{
		//double ho_energy = 0.0;
		//if (diagonal_index < dimension &&
		//    i == diagonal_indices[diagonal_index])
		//{
		//	ho_energy =
		//	       	harmonic_oscillator_energy[diagonal_index];
		//	diagonal_index++;
		//}
		double kinetic_energy = 
			(2.0)*header.basis_frequency*
			raw_elements[header.num_columns*i]/num_particles;
		assert(!isnan(kinetic_energy));
		log_entry("kinetic_energy[%lu] = %lg",i,kinetic_energy);
		log_entry("num_particles = %lu",num_particles);
		elements[3*i] = 
			raw_elements[header.num_columns*i+3]+kinetic_energy;
		elements[3*i+1] = 
			raw_elements[header.num_columns*i+5]+kinetic_energy;
		elements[3*i+2] = 
			raw_elements[header.num_columns*i+4]+kinetic_energy;
	}
	//free(diagonal_indices);
	//free(harmonic_oscillator_energy);
	return elements;
}

	static
double get_element(antoine_2nf_file_t data_file,
		   size_t bra_index,
		   size_t ket_index,
		   quantum_number Tz)
{
	if (bra_index >ket_index)
	{
		size_t tmp = bra_index;
		bra_index = ket_index;
		ket_index = tmp;
	}
	size_t element_index = get_index(data_file->used_configurations,
					 bra_index,ket_index);
	if (element_index == no_index)
		return 0.0;
	size_t type_index = Tz/2+1;
	return data_file->elements[3*element_index+type_index];
}
	static
void jump_over_fortran_block(FILE* file)
{
	int block_length = get_fortran_block_length(file);
	fseek(file,block_length,SEEK_CUR);
	close_fortran_block(file,block_length);
}

	static 
int read_fortran_block(FILE* file,
		       void *buffer,
		       size_t buffer_size)
{
	int block_length = get_fortran_block_length(file);
	if (block_length != buffer_size)
	{
		fprintf(stderr,"Warning: block and buffer are not equal\n");
	}
	if (fread(buffer,buffer_size,1,file) != 1)
	{
		return 1;
	}	
	close_fortran_block(file,block_length);
	return 0;
}

	static
int get_fortran_block_length(FILE* file)
{
	int block_length = 0;
	if (fread(&block_length,sizeof(int),1,file) != 1)
	{
		fprintf(stderr,"Could not read fortran block length\n");
		exit(EXIT_FAILURE);
	}	
	return block_length;
}

	static
void close_fortran_block(FILE* file,
			 int block_length)
{
	int block_length_final = get_fortran_block_length(file);	
	if (block_length != block_length_final)
	{
		fprintf(stderr,"Not a fortran block\n");
		exit(EXIT_FAILURE);
	}
}
