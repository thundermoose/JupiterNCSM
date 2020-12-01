#include <lanczos/lanczos.h>
#include <basis/basis.h>
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
#include <time.h>


const size_t first = 0;

struct _lanczos_environment_
{
	lanczos_settings_t settings;
	basis_t krylow_basis;
	double *diagonal_elements;
	double *off_diagonal_elements;
	size_t dimension_krylow_space;
};


static
size_t min(size_t a, size_t b);

lanczos_environment_t new_lanczos_environment(lanczos_settings_t settings)
{
	lanczos_environment_t environment = 
		(lanczos_environment_t)
		malloc(sizeof(struct _lanczos_environment_));
	environment->settings = settings;
	environment->krylow_basis = 
		new_basis_empty(settings.vector_settings,
				settings.krylow_vectors_directory_name,
				settings.max_num_iterations+1);
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
	basis_append_vector(environment->krylow_basis);
	vector_t first_krylow_vector =
	       	basis_get_vector(environment->krylow_basis,0);
	set_element(first_krylow_vector,0,1.0);
	save_vector(first_krylow_vector);
}

void lanczos_iteration(lanczos_environment_t environment,
		       size_t iteration)
{
	struct timespec t_start,t_end;
	printf("Lanczos iteration %lu start:\n",iteration+1);
	clock_gettime(CLOCK_REALTIME,&t_start);
	/* Before the Lanczos iteration can begin, it is necessary
	 * initialize the new Krylow-vector.
	 */
	basis_append_vector(environment->krylow_basis);
	/* First step of a Lanczos iteration is to compute
	 * w_k = M u_k, where M is the matrix we want to diagonalize and
	 * u_k is the k:th Krylow vector. In this implementation the next
	 * Krylow vector, u_k+1, is used as temporary storage of the helper 
	 * vectors, in order to minimize unnecessary memory usage.
	 */
	vector_t current_krylow_vector =
	       	basis_get_vector(environment->krylow_basis,
				 iteration);
	vector_t next_krylow_vector =
		basis_get_vector(environment->krylow_basis,
				 iteration+1);
	log_entry("norm of current krylow vector is %lg\n",
		  norm(current_krylow_vector));
	matrix_vector_multiplication(next_krylow_vector,
				     environment->settings.matrix,
				     current_krylow_vector);
	log_entry("norm of w is %lg\n",
		  norm(next_krylow_vector));
	log_entry("Has created a new krylow vector at index %lu",
		  iteration+1);
	/* The diagonal coefficients are the projections of the w_k on the k:th
	 * Krylow vector.
	 */
	double alpha = scalar_multiplication(current_krylow_vector,
					     next_krylow_vector);
	log_entry("alpha[%lu] = %lg\n",iteration,alpha);
	environment->diagonal_elements[iteration] = alpha;

	/* Since all Krylow vectors are orthogonal, Lanczos removes
	 * the projection of the new Krylow vector on the plane (line) of the 
	 * two (one) previous ones.
	 */ 
	if (iteration == first)
		subtract_line_projection
			(next_krylow_vector,
			 alpha,
			 current_krylow_vector);
	else
	{
		double beta_previous = 
			environment->off_diagonal_elements[iteration-1];
		vector_t previous_krylow_vector =
			basis_get_vector(environment->krylow_basis,
					 iteration-1);
		subtract_plane_projection
			(next_krylow_vector,
			 alpha,
			 current_krylow_vector,
			 beta_previous,
			 previous_krylow_vector);
	}

	double beta_new = norm(next_krylow_vector);
	environment->off_diagonal_elements[iteration] = beta_new;
	/* Normalizeing the new Krylow vector.
	*/
	scale(next_krylow_vector, 1.0 / beta_new);
	log_entry("norm of new krylow vector is %lg",
		  norm(next_krylow_vector));
	log_entry("scale %lg",
		  scalar_multiplication(next_krylow_vector,
					current_krylow_vector));
	clock_gettime(CLOCK_REALTIME,&t_end);
	double iteration_time = 
		(t_end.tv_sec-t_start.tv_sec)*1e6 +
		(t_end.tv_nsec-t_start.tv_nsec)*1e-3;
	printf("Lanczos iteration %lu end after %lg µs\n", 
	       iteration+1, iteration_time);
}

void orthogonalize_krylow_basis(lanczos_environment_t environment,
				size_t iteration)
{
	if (iteration > 2)
	{
		vector_t current_vector =
		       	basis_get_vector(environment->krylow_basis,
					 iteration + 1);
		vector_t *all_krylow_vectors =
		       	basis_get_all_vectors(environment->krylow_basis);
		reorthogonalize_vector(current_vector, 
				       all_krylow_vectors,
				       iteration-2);
	}
}

void diagonalize(lanczos_environment_t environment)
{
	struct timespec t_start,t_end;
	printf("Diagonalzation start:\n");
	clock_gettime(CLOCK_REALTIME,&t_start);
	initialize_first_krylow_vector(environment);
	log_entry("environment->settings.max_num_iterations = %lu",
		  environment->settings.max_num_iterations);
	size_t max_num_iterations =
	       	min(environment->settings.dimension,
		    environment->settings.max_num_iterations);
	size_t iteration;
	for (iteration = 0;
	     iteration < max_num_iterations;
	     iteration++)
	{
		lanczos_iteration(environment, iteration);
		orthogonalize_krylow_basis(environment, iteration);
	}
	basis_remove_last(environment->krylow_basis);
	clock_gettime(CLOCK_REALTIME,&t_end);
	double diagonalzation_time =
		(t_end.tv_sec-t_start.tv_sec)*1e6 +
		(t_end.tv_nsec-t_start.tv_nsec)*1e-3;
	printf("Diagonalization end after %lg µs\n",
	       diagonalzation_time);
}

eigen_system_t get_eigensystem(lanczos_environment_t environment)
{
	eigen_system_t diagonalized_system =
		diagonalize_tridiagonal_matrix
		(environment->diagonal_elements,
		 environment->off_diagonal_elements,
		 basis_get_dimension(environment->krylow_basis));
	set_basis(diagonalized_system,
		  environment->krylow_basis);
	return diagonalized_system;
}

void free_lanczos_environment(lanczos_environment_t environment)
{
	free_basis(environment->krylow_basis);
	free(environment->diagonal_elements);
	free(environment->off_diagonal_elements);
	free(environment);
}

static
size_t min(size_t a, size_t b)
{
	return a<b ? a : b;
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
