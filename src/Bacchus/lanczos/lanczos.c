#include <lanczos/lanczos.h>
#include <vector/vector.h>
#include <diagonalization/diagonalization.h>
#include <string_tools/string_tools.h>
#include <string.h>
#include <directory_tools/directory_tools.h>
#include <error/error.h>
#include <unit_testing/test.h>
#include <log/log.h>
#include <math.h>
#include <errno.h>


const size_t first = 0;

struct _lanczos_environment_
{
	lanczos_settings_t settings;
	vector_t *krylow_vectors;
	double *diagonal_elements;
	double *off_diagonal_elements;
	size_t dimension_krylow_space;
};


static
void new_krylow_vector(lanczos_environment_t environment,
			   size_t index);

static
vector_settings_t get_krylow_vector_settings(lanczos_environment_t environment,
					     size_t index);

lanczos_environment_t new_lanczos_environment(lanczos_settings_t settings)
{
	lanczos_environment_t environment = 
		(lanczos_environment_t)
		malloc(sizeof(struct _lanczos_environment_));
	environment->settings = settings;
	environment->krylow_vectors = 
		(vector_t*)calloc(settings.max_num_iterations+1,
				  sizeof(vector_t));
	environment->diagonal_elements = 
		(double*)malloc(settings.max_num_iterations*sizeof(double));
	environment->off_diagonal_elements = 
		(double*)malloc(settings.max_num_iterations*sizeof(double));
	environment->dimension_krylow_space = 0;
	if (!directory_exists(settings.krylow_vectors_directory_name) &&
		create_directory(settings.krylow_vectors_directory_name) != 0)
		error("Could not create krylow vector directory \"%s\". %s\n",
		      settings.krylow_vectors_directory_name,
		      strerror(errno));
	return environment;
}

void initialize_first_krylow_vector(lanczos_environment_t environment)
{
	/* The first Krylow vector for lanczos is a choice. 
	 * In this early version we choose to set the first component of the 
	 * Krylow vector to 1 and the rest to 0.
	 * However, better choices exists and for the future this might change
	 */
	new_krylow_vector(environment,0);
	set_element(environment->krylow_vectors[0],first,1);
	save_vector(environment->krylow_vectors[0]);
	environment->dimension_krylow_space = 1;
}

void lanczos_iteration(lanczos_environment_t environment,
		       size_t iteration)
{

	/* Before the Lanczos iteration can begin, it is necessary
	 * initialize the new Krylow-vector.
	 */
	new_krylow_vector(environment,iteration+1);
	/* First step of a Lanczos iteration is to compute
	 * w_k = M u_k, where M is the matrix we want to diagonalize and
	 * u_k is the k:th Krylow vector. In this implementation the next
	 * Krylow vector, u_k+1, is used as temporary storage of the helper 
	 * vectors, in order to minimize unnecessary memory usage.
	 */

	log_entry("norm of current krylow vector is %lg\n",
		  norm(environment->krylow_vectors[iteration]));
	matrix_vector_multiplication(environment->krylow_vectors[iteration+1],
				     environment->settings.matrix,
				     environment->krylow_vectors[iteration]);
	log_entry("norm of w is %lg\n",
		  norm(environment->krylow_vectors[iteration+1]));
	log_entry("Has created a new krylow vector at index %lu",
		  iteration+1);
	environment->dimension_krylow_space++;		
	/* The diagonal coefficients are the projections of the w_k on the k:th
	 * Krylow vector.
	 */
	double alpha = 
		scalar_multiplication
		(environment->krylow_vectors[iteration],
		 environment->krylow_vectors[iteration+1]);
	log_entry("alpha[%lu] = %lg\n",iteration,alpha);
	environment->diagonal_elements[iteration] = alpha;

	/* Since all Krylow vectors are orthogonal, Lanczos removes
	 * the projection of the new Krylow vector on the plane (line) of the 
	 * two (one) previous ones.
	 */ 
	if (iteration == first)
		subtract_line_projection
			(environment->krylow_vectors[iteration+1],
			 alpha,
			 environment->krylow_vectors[iteration]);
	else
	{
		double beta_previous = 
			environment->off_diagonal_elements[iteration-1];
		subtract_plane_projection
			(environment->krylow_vectors[iteration+1],
			 alpha,
			 environment->krylow_vectors[iteration],
			 beta_previous,
			 environment->krylow_vectors[iteration-1]);
	}

	double beta_new = 
		norm(environment->krylow_vectors[iteration+1]);
	environment->off_diagonal_elements[iteration] = beta_new;
	/* Normalizeing the new Krylow vector.
	*/
	scale(environment->krylow_vectors[iteration+1],
	      1.0 / beta_new);
	log_entry("norm of new krylow vector is %lg",
		  norm(environment->krylow_vectors[iteration+1]));
	log_entry("scale %lg",
		  scalar_multiplication(environment->krylow_vectors[iteration+1],
					environment->krylow_vectors[iteration]));
}

void orthogonalize_krylow_basis(lanczos_environment_t environment,
				size_t iteration)
{
	if (iteration > 2)
		reorthogonalize_vector(environment->krylow_vectors[iteration+1],
				       environment->krylow_vectors,
				       iteration-2);
}

void diagonalize(lanczos_environment_t environment)
{
	initialize_first_krylow_vector(environment);
	log_entry("environment->settings.max_num_iterations = %lu",
		  environment->settings.max_num_iterations);
	size_t iteration;
	for (iteration = 0;
	     iteration < environment->settings.max_num_iterations;
	     iteration++)
	{
		lanczos_iteration(environment, iteration);
		orthogonalize_krylow_basis(environment, iteration);
	}
	environment->dimension_krylow_space--;
}

eigen_system_t get_eigensystem(lanczos_environment_t environment)
{
	return diagonalize_tridiagonal_matrix
		(environment->diagonal_elements,
		 environment->off_diagonal_elements,
		 environment->dimension_krylow_space);
}

void free_lanczos_environment(lanczos_environment_t environment)
{
	size_t i;
	for (i = 0; i<environment->settings.max_num_iterations+1; i++)
	{
		if (environment->krylow_vectors[i] != NULL)
		{
			free_vector(environment->krylow_vectors[i]);
			log_entry("freed a krylow_vector at index %lu",i);
		}
	}
	free(environment->krylow_vectors);
	free(environment->diagonal_elements);
	free(environment->off_diagonal_elements);
	free(environment);
}

static
void new_krylow_vector(lanczos_environment_t environment,
			   size_t index)
{
	vector_settings_t settings = get_krylow_vector_settings(environment,
								index);
	environment->krylow_vectors[index] =
		new_zero_vector(settings);
	free(settings.directory_name);
}

static
vector_settings_t get_krylow_vector_settings(lanczos_environment_t environment,
					     size_t index)
{
	vector_settings_t settings = environment->settings.vector_settings;
	char directory_name[2048] = {0};
	sprintf(directory_name,
		"%s/krylow_vector_%lu",
		environment->settings.krylow_vectors_directory_name,
		index);
	settings.directory_name = copy_string(directory_name);
	return settings;
}

new_test(diagonalize_3x3_matrix,
	 {
	 double matrix_elements[9] = 
	 {
	 1,2,3,
	 2,1,2,
	 3,2,1
	 };
	 double expected_eigen_values[3] = 
	 {
	 -2,
	 5.0/2-sqrt(41)/2,
	 5.0/2+sqrt(41)/2
	 };
	 matrix_t simple_3x3_matrix = new_zero_matrix(3,3);
	 set_matrix_elements(simple_3x3_matrix,
			     matrix_elements);
	 size_t dimension = 3;
	 lanczos_settings_t settings =
	 {
	 .dimension = dimension,
	 .vector_settings =
	 {
	 	.directory_name = NULL,
		.num_blocks = 1,
		.block_sizes = &dimension
	 },
	 .krylow_vectors_directory_name = copy_string(get_test_file_path("krylow_vectors")),
	 .max_num_iterations = 3,
	 .target_eigenvalue = 0,
	 .eigenvalue_tollerance = 1e-5,
	 .matrix = simple_3x3_matrix
	 };
	 lanczos_environment_t environment = 
		 new_lanczos_environment(settings);
	 diagonalize(environment);
	 eigen_system_t eigen_system = get_eigensystem(environment);
	 printf("get_num_eigen_values(eigen_system) = %lu\n",
		get_num_eigen_values(eigen_system));
	 assert_that(get_num_eigen_values(eigen_system) == 3);
	 for (size_t i = 0; i<3; i++)
	 {
		 double computed_eigen_value = 
			 get_eigen_value(eigen_system,i);
		 double expected_eigen_value = 
			 expected_eigen_values[i];
		 printf("%lu: is %lg but expected is %lg\n",
			i,computed_eigen_value,
			expected_eigen_value);
		 assert_that(fabs(computed_eigen_value
				  -expected_eigen_value)<1e-5);
	 }		 
	 free_lanczos_environment(environment);
	 free_matrix(simple_3x3_matrix);
	 free_eigen_system(eigen_system);
	 free(settings.krylow_vectors_directory_name);
	 });

new_test(diagonalize_50x50_random_matrix_and_compare_to_lapack,
	 {
	 matrix_t matrix_50x50 = new_random_symmetric_matrix(50);
	 size_t dimension = 50;
	 lanczos_settings_t settings =
	 {
	 .dimension = dimension,
	 .vector_settings =
	 {
	 	.directory_name = NULL,
		.num_blocks = 1,
		.block_sizes = &dimension
	 },
	 .krylow_vectors_directory_name = copy_string(get_test_file_path("krylow_vectors")),
	 .max_num_iterations = dimension,
	 .target_eigenvalue = 0,
	 .eigenvalue_tollerance = 1e-5,
	 .matrix = matrix_50x50
	 };
	 lanczos_environment_t environment = 
	 new_lanczos_environment(settings);
	 diagonalize(environment);
	 eigen_system_t lanczos_eigen_system = 
	 get_eigensystem(environment);

	 eigen_system_t lapack_eigen_system =
	 diagonalize_symmetric_matrix(matrix_50x50);

	 assert_that(get_num_eigen_values(lanczos_eigen_system) == 
		     get_num_eigen_values(lapack_eigen_system));
	 int test_passed = 1;
	 for (size_t i = 0; i<get_num_eigen_values(lanczos_eigen_system); i++)
	 {
		 double lanczos_eigen_value = 
			 get_eigen_value(lanczos_eigen_system,i);
		 double lapack_eigen_value = 
			 get_eigen_value(lapack_eigen_system,i);
		 printf("(%lu) lanczos: %lg, lapack: %lg\n",
			i,lanczos_eigen_value,lapack_eigen_value);
		 if (fabs(lanczos_eigen_value-lapack_eigen_value) > 1e-10)
			 test_passed = 0;
	 }
	 assert_that(test_passed);
	 free_lanczos_environment(environment);
	 free_eigen_system(lanczos_eigen_system);
	 free_eigen_system(lapack_eigen_system);
	 free_matrix(matrix_50x50);
	 free(settings.krylow_vectors_directory_name);
	 });


new_test_silent(diagonalize_500x500_random_matrix_and_compare_to_lapack,
	 matrix_t matrix_500x500 = new_random_symmetric_matrix(500);
	 size_t dimension = 500;
	 lanczos_settings_t settings =
	 {
	 .dimension = dimension,
	 .vector_settings =
	 {
	 	.directory_name = NULL,
		.num_blocks = 1,
		.block_sizes = &dimension
	 },
	 .krylow_vectors_directory_name = copy_string(get_test_file_path("krylow_vectors")),
	 .max_num_iterations = dimension,
	 .target_eigenvalue = 0,
	 .eigenvalue_tollerance = 1e-5,
	 .matrix = matrix_500x500
	 };
	 lanczos_environment_t environment = 
	 new_lanczos_environment(settings);
	 diagonalize(environment);
	 eigen_system_t lanczos_eigen_system = 
	 get_eigensystem(environment);
	 free_lanczos_environment(environment);

	 eigen_system_t lapack_eigen_system =
	 diagonalize_symmetric_matrix(matrix_500x500);

	 assert_that(get_num_eigen_values(lanczos_eigen_system) == 
		     get_num_eigen_values(lapack_eigen_system));
	 int test_passed = 1;
	 for (size_t i = 0; i<get_num_eigen_values(lanczos_eigen_system); i++)
	 {
		 double lanczos_eigen_value = 
			 get_eigen_value(lanczos_eigen_system,i);
		 double lapack_eigen_value = 
			 get_eigen_value(lapack_eigen_system,i);
		 printf("(%lu) lanczos: %lg, lapack: %lg\n",
			i,lanczos_eigen_value,lapack_eigen_value);
		 if (fabs(lanczos_eigen_value-lapack_eigen_value) > 1e-10)
			 test_passed = 0;
	 }
	 assert_that(test_passed);
	 free_eigen_system(lanczos_eigen_system);
	 free_eigen_system(lapack_eigen_system);
	 free_matrix(matrix_500x500);
	 );
