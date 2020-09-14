#include <vector/vector.h>
#include <array_builder/array_builder.h>
#include <directory_tools/directory_tools.h>
#include <math_tools/math_tools.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <math.h>
#include <unit_testing/test.h>
#include <errno.h>
#include <string.h>

typedef struct
{
	size_t start_index;
	size_t block_length;
	size_t block_id;
} vector_block_t;

struct _vector_
{
	size_t dimension;
	char *directory_name;
	vector_block_t *vector_blocks;
	size_t num_vector_blocks;
	vector_block_t loaded_block;
	double *element_buffer;
	size_t element_buffer_length;
};

const size_t no_index = -1;

static
void initiate_vector_file(vector_t vector,
			  vector_block_t vector_block);
static
vector_block_t find_block(vector_t vector,size_t index);

static
void prepare_element_buffer(vector_t vector,
			    vector_block_t vector_block);

static
FILE* open_block_file(vector_t vector,
		      vector_block_t block,
		      const char *file_mode);

#ifndef NDEBUG
static
int compare_blocks(vector_block_t first_block,
		   vector_block_t second_block);
#endif

static
void expand_element_buffer(double **buffer,
			   size_t *buffer_length,
			   size_t needed_buffer_length);

static
void load_vector_elements(double *vector_elements,
			  vector_t vector,
			  vector_block_t vector_block);

static
void save_vector_elements(double *vector_elements,
			  vector_t vector,
			  vector_block_t vector_block);
vector_settings_t setup_vector_settings(combination_table_t combination_table)
{
	reset_basis_block_iterator(combination_table);
	vector_settings_t settings = 
	{
		.directory_name = NULL,
		.num_blocks = 0,	
		.block_sizes = NULL
	};
	array_builder_t block_sizes_builder =
		new_array_builder((void**)&settings.block_sizes,
				  &settings.num_blocks,
				  sizeof(size_t));
	while (has_next_basis_block(combination_table))
	{
		basis_block_t current_basis_block =
			next_basis_block(combination_table);
		size_t total_size = 
			current_basis_block.num_proton_states*
			current_basis_block.num_neutron_states;
		append_array_element(block_sizes_builder,
				     &total_size);
	}
	free_array_builder(block_sizes_builder);
	return settings;
}

vector_t new_zero_vector(vector_settings_t vector_settings)
{
	if (!directory_exists(vector_settings.directory_name) &&
	    create_directory(vector_settings.directory_name) != 0)
		error("Could not create directory \"%s\". %s\n",
		      vector_settings.directory_name,
		      strerror(errno));
	vector_t vector = (vector_t)calloc(1,sizeof(struct _vector_));
	vector->dimension = sum_sizes(vector_settings.block_sizes,
				      vector_settings.num_blocks);
	vector->num_vector_blocks = vector_settings.num_blocks;
	vector->vector_blocks =
		(vector_block_t*)calloc(vector->num_vector_blocks,
					sizeof(vector_block_t));
	size_t start_index = 0;
	vector->directory_name = copy_string(vector_settings.directory_name);
	for (size_t i = 0; i<vector->num_vector_blocks; i++)
	{
		vector_block_t vector_block =
		{
			.start_index = start_index,
			.block_length = vector_settings.block_sizes[i],
			.block_id = i+1
		};
		vector->vector_blocks[i] = vector_block;
		start_index += vector_settings.block_sizes[i];
		initiate_vector_file(vector,vector_block);
	}
	vector->loaded_block.block_id = -1;
	return vector;
}

vector_t new_random_vector(vector_settings_t vector_settings)
{
	vector_t vector = new_zero_vector(vector_settings);
	for (size_t i = 0; i<vector->dimension; i++)
		set_element(vector,i,2*drand48()-1);
	save_vector(vector);
	return vector;
}

void set_element(vector_t vector,
		 size_t index,
		 double value)
{
	log_entry("Setting element %lu of vector %s to %lg",
		  index,vector->directory_name,value);
	vector_block_t vector_block = find_block(vector,index);
	if (vector_block.block_id != vector->loaded_block.block_id)
	{
		if (vector->loaded_block.block_id < vector->num_vector_blocks)
			save_vector(vector);
		prepare_element_buffer(vector,vector_block);
		load_vector_elements(vector->element_buffer,
				     vector,
				     vector_block);
		vector->loaded_block = vector_block;
	}
	vector->element_buffer[index-vector_block.start_index] = value;
}

double get_element(vector_t vector,
		   size_t index)
{
	log_entry("Getting element %lu of vector %s",
		  index,vector->directory_name);
	vector_block_t vector_block = find_block(vector,index);
	if (vector_block.block_id != vector->loaded_block.block_id)
	{
		if (vector->loaded_block.block_id < vector->num_vector_blocks)
			save_vector(vector);
		prepare_element_buffer(vector,vector_block);
		load_vector_elements(vector->element_buffer,
				     vector,
				     vector_block);
		vector->loaded_block = vector_block;
	}
	return vector->element_buffer[index-vector_block.start_index];
}

const char *get_vector_path(vector_t vector)
{
	return vector->directory_name;
}

void save_vector(vector_t vector)
{
	log_entry("Saving vector %s",
		  vector->directory_name);
	save_vector_elements(vector->element_buffer,
			     vector,
			     vector->loaded_block);
	vector->loaded_block.block_id = -1;
}

void print_vector(vector_t vector)
{
	printf("%s\n",
	       vector->directory_name);

	for (size_t i = 0; i<vector->dimension; i++)
		printf("%lg ",get_element(vector,i));
	printf("\n");
}

size_t vector_dimension(vector_t vector)
{
	assert(vector != NULL);
	return vector->dimension;
}

double scalar_multiplication(const vector_t first_vector,
			     const vector_t second_vector)
{
	log_entry("scalar_multiplication(%s,%s)",
		  first_vector->directory_name,
		  second_vector->directory_name);
	assert(first_vector != NULL);
	assert(second_vector != NULL);
	assert(first_vector->dimension == second_vector->dimension);
	assert(first_vector->num_vector_blocks ==
	       second_vector->num_vector_blocks);
	double accumulator = 0.0;
	double *element_buffer = NULL;
	size_t element_buffer_length = 0;
	for (size_t i = 0; i <first_vector->num_vector_blocks; i++)
	{
		log_entry("block_id = %lu",i);
		vector_block_t vector_block = first_vector->vector_blocks[i];
		assert(compare_blocks(vector_block,
				      second_vector->vector_blocks[i]));
		expand_element_buffer(&element_buffer,
				      &element_buffer_length,
				      2*vector_block.block_length);
		double *first_vector_elements = element_buffer;
		double *second_vector_elements =
			element_buffer + vector_block.block_length;
		load_vector_elements(first_vector_elements,
				     first_vector,vector_block);
		load_vector_elements(second_vector_elements,
				     second_vector,vector_block);
		double block_scalar_product = array_scalar_product(first_vector_elements,
						    second_vector_elements,
						    vector_block.block_length);
		log_entry("block_scalar_product = %lg",
			  block_scalar_product);
		log_entry("accumulator = %lg",accumulator);
		accumulator += block_scalar_product;
	}
	if (element_buffer != NULL)
		free(element_buffer);
	log_entry("accumulator = %lg",accumulator);
	return accumulator;
}


void subtract_line_projection(vector_t target_vector,
			      double projection,
			      const vector_t line_direction)
{
	log_entry("Subtracting line projection");
	assert(target_vector != NULL);
	assert(line_direction != NULL);
	assert(target_vector->dimension == line_direction->dimension);
	assert(target_vector->num_vector_blocks ==
	       line_direction->num_vector_blocks);
	double *element_buffer = NULL;
	size_t buffer_length = 0;
	for (size_t i = 0; i<target_vector->num_vector_blocks; i++)
	{
		vector_block_t vector_block = target_vector->vector_blocks[i];
		assert(compare_blocks(vector_block,
				      line_direction->vector_blocks[i]));
		expand_element_buffer(&element_buffer,
				      &buffer_length,
				      2*vector_block.block_length);
		double *target_vector_elements = element_buffer;
		double *line_direction_elements =
			element_buffer + vector_block.block_length;
		load_vector_elements(target_vector_elements,
				     target_vector,
				     vector_block);
		load_vector_elements(line_direction_elements,
				     line_direction,
				     vector_block);
		subtract_array_projection(target_vector_elements,
					  projection,
					  line_direction_elements,
					  vector_block.block_length);
		save_vector_elements(target_vector_elements,
				     target_vector,
				     vector_block);
	}
	if (element_buffer != NULL)
		free(element_buffer);
}

void subtract_plane_projection(vector_t target_vector,
			       double first_projection,
			       const vector_t first_direction,
			       double second_projection,
			       const vector_t second_direction)
{
	log_entry("Subtracting plane projection");
	assert(target_vector != NULL);
	assert(first_direction != NULL);
	assert(second_direction != NULL);
	assert(target_vector->dimension == first_direction->dimension);
	assert(target_vector->dimension == second_direction->dimension);
	assert(target_vector->num_vector_blocks ==
	       first_direction->num_vector_blocks);
	assert(target_vector->num_vector_blocks ==
	       second_direction->num_vector_blocks);
	double *element_buffer = NULL;
	size_t buffer_length = 0;
	for (size_t i = 0; i<target_vector->num_vector_blocks; i++)
	{
		vector_block_t vector_block = target_vector->vector_blocks[i];
		assert(compare_blocks(vector_block,
				      first_direction->vector_blocks[i]));
		assert(compare_blocks(vector_block,
				      second_direction->vector_blocks[i]));
		expand_element_buffer(&element_buffer,
				      &buffer_length,
				      3*vector_block.block_length);
		double *target_vector_elements = element_buffer;
		double *first_direction_elements =
			element_buffer + vector_block.block_length;
		double *second_direction_elements =
			element_buffer + 2*vector_block.block_length;
		load_vector_elements(target_vector_elements,
				     target_vector,
				     vector_block);
		load_vector_elements(first_direction_elements,
				     first_direction,
				     vector_block);
		load_vector_elements(second_direction_elements,
				     second_direction,
				     vector_block);
		subtract_array_projection(target_vector_elements,
					  first_projection,
					  first_direction_elements,
					  vector_block.block_length);
		subtract_array_projection(target_vector_elements,
					  second_projection,
					  second_direction_elements,
					  vector_block.block_length);
		save_vector_elements(target_vector_elements,
				     target_vector,
				     vector_block);
	}
	if (element_buffer != NULL)
		free(element_buffer);
}

double norm(const vector_t vector)
{
	log_entry("Computing norm of vector %s",
		  vector->directory_name);
	assert(vector != NULL);
	double *element_buffer = NULL;
	size_t element_buffer_size = 0;
	double accumulator = 0.0;
	for (size_t i = 0; i<vector->num_vector_blocks; i++)
	{
		vector_block_t vector_block = vector->vector_blocks[i];
		expand_element_buffer(&element_buffer,
				      &element_buffer_size,
				      vector_block.block_length);
		load_vector_elements(element_buffer,
				     vector,
				     vector_block);
		accumulator+=array_square_norm(element_buffer,
					       vector_block.block_length);
	}
	free(element_buffer);
	return sqrt(accumulator);
}

void scale(vector_t vector,double scaling)
{
	log_entry("Scaling vector %s",
		  vector->directory_name);
	assert(vector != NULL);
	double *element_buffer = NULL;
	size_t element_buffer_size = 0;
	for (size_t i = 0; i<vector->num_vector_blocks; i++)
	{
		vector_block_t vector_block = vector->vector_blocks[i];
		expand_element_buffer(&element_buffer,
				      &element_buffer_size,
				      vector_block.block_length);
		load_vector_elements(element_buffer,
				     vector,
				     vector_block);
		scale_array(element_buffer,scaling,vector_block.block_length);
		save_vector_elements(element_buffer,
				     vector,
				     vector_block);
	}
	free(element_buffer);
}

void reorthogonalize_vector(vector_t vector_to_orthogonalize,
			    vector_t *basis,
			    size_t num_basis_states)
{
	log_entry("reorthogonalize vector %s",
		  vector_to_orthogonalize->directory_name);
	for (size_t i = 0; i<num_basis_states; i++)
	{
		double projection =
			scalar_multiplication(vector_to_orthogonalize,
					      basis[i]);
		subtract_line_projection(vector_to_orthogonalize,
					 projection,
					 basis[i]);
	}
} 	

void free_vector(vector_t vector)
{
	log_entry("free_vector: %p",vector);
	free(vector->directory_name);
	free(vector->vector_blocks);
	if (vector->element_buffer != NULL)
		free(vector->element_buffer);
	free(vector);
}

	static
void initiate_vector_file(vector_t vector,
			  vector_block_t vector_block)
{
	log_entry("initiate vector block %lu of %s",
		  vector_block.block_id,
		  vector->directory_name);
	FILE* vector_file = open_block_file(vector,vector_block,"w");
	double *array = (double*)calloc(vector_block.block_length,
					sizeof(double));
	if (fwrite(array,
		   sizeof(double),
		   vector_block.block_length,
		   vector_file) != vector_block.block_length)
		error("Could not write to vector file\n");
	fclose(vector_file);
	free(array);
}

	static
vector_block_t find_block(vector_t vector,size_t index)
{
	log_entry("Looking for vector block containing index %lu",
		  index);
	for (size_t i = 0; i<vector->num_vector_blocks; i++)	
		if (vector->vector_blocks[i].start_index <= index &&
		    vector->vector_blocks[i].start_index +
		    vector->vector_blocks[i].block_length > index)
			return vector->vector_blocks[i];
	vector_block_t not_a_block = {no_index,no_index,no_index};
	return not_a_block;
}

static
void prepare_element_buffer(vector_t vector,
			    vector_block_t vector_block)
{
	if (vector->element_buffer_length < vector_block.block_length)
	{
		vector->element_buffer_length = vector_block.block_length;
		vector->element_buffer =
			(double*)realloc(vector->element_buffer,
					 sizeof(double)*
					 vector->element_buffer_length);
	}	
	log_entry("element_buffer_length = %lu",
		  vector->element_buffer_length);
}

	static
FILE* open_block_file(vector_t vector,
		      vector_block_t vector_block,
		      const char *file_mode)
{
	char file_name[2049];
	sprintf(file_name,
		"%s/vec_%lu",
		vector->directory_name,
		vector_block.block_id);
	log_entry("Open vector file %s with file mode %s",
		  file_name,file_mode);	  
	FILE *vector_file = fopen(file_name,file_mode);
	if (vector_file == NULL)
		error("Could not open vector file %s. %s\n",
		      file_name,
		      strerror(errno));
	return vector_file;
}

#ifndef NDEBUG
	static
int compare_blocks(vector_block_t first_block,
		   vector_block_t second_block)
{
	return first_block.start_index == second_block.start_index &&
		first_block.block_length == second_block.block_length &&
		first_block.block_id == second_block.block_id;
}
#endif

	static
void expand_element_buffer(double **buffer,
			   size_t *buffer_length,
			   size_t needed_buffer_length)
{
	if (needed_buffer_length > *buffer_length)
	{
		*buffer_length = needed_buffer_length;
		// We could use realloc here, but then
		// there will be unnecessary memory coping
		if (*buffer != NULL)
			free(*buffer);
		*buffer = (double*)malloc((*buffer_length)*
					  sizeof(double));
	}
}

	static
void load_vector_elements(double *vector_elements,
			  vector_t vector,
			  vector_block_t vector_block)
{
	assert(vector_elements != NULL);
	FILE* vector_file = open_block_file(vector,vector_block,"r");
	size_t read_length = fread(vector_elements,
				   sizeof(double),
				   vector_block.block_length,
				   vector_file);
	if (read_length != vector_block.block_length)
		error("Could not read vector elements. Read %lu elements\n",
		      read_length);
	fclose(vector_file);
}

	static
void save_vector_elements(double *vector_elements,
			  vector_t vector,
			  vector_block_t vector_block)
{
	FILE* vector_file = open_block_file(vector,vector_block,"w");
	log_entry("saving vector %s:%lu with length %lu",
		  vector->directory_name,
		  vector_block.block_id,
		  vector_block.block_length);
	if (fwrite(vector_elements,
		   sizeof(double),
		   vector_block.block_length,
		   vector_file) != vector_block.block_length)
		error("Could not read vector elements\n");
	fclose(vector_file);
}

new_test(new_zero_vector,
	 {
	 const size_t desired_dimension = 100;
	 size_t block_sizes[2] =
	 {
	 desired_dimension/2,
	 desired_dimension/2
	 };
	 vector_settings_t settings =
	 {
	 .directory_name = copy_string(get_test_file_path("vector")),
	 .num_blocks = 2,
	 .block_sizes = block_sizes
	 };
	 vector_t vector = new_zero_vector(settings);
	 assert_that(vector != NULL);
	 assert_that(vector_dimension(vector) == desired_dimension);
	 for (size_t i = 0; i<desired_dimension; i++)
	 {
	 assert_that(fabs(get_element(vector,i))<1e-10);
	 }	
	free_vector(vector);
	free(settings.directory_name);
	 }
);

new_test(norm_of_unitary_2d_vector,
	 {	
	 size_t block_size = 2;
	 vector_settings_t settings =
	 {
	 .directory_name = copy_string(get_test_file_path("vector")),
	 .num_blocks = 1,
	 .block_sizes = &block_size
	 };
	 vector_t vector = new_zero_vector(settings);
	 set_element(vector,0,1.0/sqrt(2));
	 set_element(vector,1,1.0/sqrt(2));
	 save_vector(vector);
	 assert_that(fabs(norm(vector)-1)<1e-10);
	 free_vector(vector);
	 free(settings.directory_name);
	 }
	);

new_test(scalar_product_of_orthogonal_vectors,
	 size_t dimension = 2;
	 vector_settings_t settings_1 =
	 {
	 	.directory_name = copy_string(get_test_file_path("vector_1")),
		.num_blocks = 1,
		.block_sizes = &dimension
	 };
	 vector_settings_t settings_2 =
	 {
	 	.directory_name = copy_string(get_test_file_path("vector_2")),
		.num_blocks = 1,
		.block_sizes = &dimension
	 };
	 vector_t first_vector = new_zero_vector(settings_1);
	 vector_t second_vector = new_zero_vector(settings_2);
	 set_element(first_vector,0,1.0);
	 set_element(second_vector,1,1.0);
	 assert_that(fabs(scalar_multiplication(first_vector,
						second_vector)) < 1e-10);
	 free_vector(first_vector);
	 free_vector(second_vector);
	 free(settings_1.directory_name);
	 free(settings_2.directory_name);

	);

new_test(remove_projection_on_vector,
	 size_t dimension = 100;
	 size_t num_blocks = 5;
	 size_t *block_sizes = (size_t*)malloc(num_blocks*sizeof(size_t));
	 for (size_t i = 0; i<num_blocks; i++)
	 	block_sizes[i] = dimension/num_blocks;
	 vector_settings_t line_direction_settings =
	 {
		.directory_name =
	       	copy_string(get_test_file_path("line_direction")),
		.num_blocks = num_blocks,
		.block_sizes = block_sizes
	 };
	 vector_settings_t projected_vector_settings =
	 {
		.directory_name =
	       	copy_string(get_test_file_path("projected_vector")),
		.num_blocks = num_blocks,
		.block_sizes = block_sizes
	 };
	 vector_t line_direction =
		 new_random_vector(line_direction_settings);
	 vector_t projected_vector =
	 	new_random_vector(projected_vector_settings);
	 double line_direction_norm = norm(line_direction);
	 log_entry("line_direction_norm = %lg",
		   line_direction_norm);
	 scale(line_direction,1/line_direction_norm);
	 double projection =
		 scalar_multiplication(projected_vector,line_direction);
	 log_entry("projection = %lg",projection);
	 subtract_line_projection(projected_vector,
				  projection,
				  line_direction);
         assert_that(fabs(scalar_multiplication(projected_vector,
						line_direction)) <1e-10);
	 free_vector(line_direction);
	 free_vector(projected_vector);
	 free(block_sizes);
	 free(line_direction_settings.directory_name);
	 free(projected_vector_settings.directory_name);
	);

new_test(remove_plane_projection_on_vector,
	 size_t dimension = 100;
	 size_t num_blocks = 5;
	 size_t *block_sizes = (size_t*)malloc(num_blocks*sizeof(size_t));
	 for (size_t i = 0; i<num_blocks; i++)
	 	block_sizes[i] = dimension/num_blocks;
	 vector_settings_t plane_direction_one_settings =
	 {
		.directory_name = copy_string(get_test_file_path("plane_direction_one")),
		.num_blocks = num_blocks,
		.block_sizes = block_sizes
	 };
	 vector_settings_t plane_direction_two_settings =
	 {
		.directory_name = copy_string(get_test_file_path("plane_direction_two")),
		.num_blocks = num_blocks,
		.block_sizes = block_sizes
	 };
	 vector_settings_t projected_vector_settings =
	 {
		.directory_name = copy_string(get_test_file_path("projected_vector")),
		.num_blocks = num_blocks,
		.block_sizes = block_sizes
	 };
	 vector_t plane_direction_one =
		 new_random_vector(plane_direction_one_settings);
	 vector_t plane_direction_two =
		 new_random_vector(plane_direction_two_settings);
	 vector_t projected_vector =
	 	new_random_vector(projected_vector_settings);
	 double plane_direction_one_norm = norm(plane_direction_one);
	 log_entry("plane_direction_one_norm = %lg",
		   plane_direction_one_norm);
	 scale(plane_direction_one,1/plane_direction_one_norm);
	 double two_on_one_projection = 
	 	scalar_multiplication(plane_direction_one,
				      plane_direction_two);
	 subtract_line_projection(plane_direction_two,
				  two_on_one_projection,
				  plane_direction_one);
	 double plane_direction_two_norm = norm(plane_direction_two);
	 log_entry("plane_direction_two_norm = %lg",
		   plane_direction_two_norm);
	 log_entry("pd1 | pd2 = %lg",
		   scalar_multiplication(plane_direction_one,
					 plane_direction_two));
	 scale(plane_direction_two,1/plane_direction_two_norm);
	 assert_that(fabs(norm(plane_direction_two)-1)<1e-10);
	 log_entry("pd1 | pd2 = %lg",
		   scalar_multiplication(plane_direction_one,
					 plane_direction_two));
	 assert_that(fabs(scalar_multiplication(plane_direction_one,
						plane_direction_two)) < 1e-10);
	 double projection_one =
		 scalar_multiplication(projected_vector,plane_direction_one);
	 double projection_two =
		 scalar_multiplication(projected_vector,plane_direction_two);
	 log_entry("projection = (%lg,%lg)",projection_one,projection_two);
	 subtract_plane_projection(projected_vector,
				  projection_two,
				  plane_direction_two,
				  projection_one,
				  plane_direction_one);
	 log_entry("scalar products: %lg %lg",
		scalar_multiplication(projected_vector,
				      plane_direction_one),
		scalar_multiplication(projected_vector,
				      plane_direction_two));
         assert_that(fabs(scalar_multiplication(projected_vector,
						plane_direction_one)) <1e-10);
         assert_that(fabs(scalar_multiplication(projected_vector,
						plane_direction_two)) <1e-10);
	 free_vector(plane_direction_one);
	 free_vector(plane_direction_two);
	 free_vector(projected_vector);
	 free(block_sizes);
	 free(plane_direction_one_settings.directory_name);
	 free(plane_direction_two_settings.directory_name);
	 free(projected_vector_settings.directory_name);
	);

