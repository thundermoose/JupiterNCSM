#ifndef __LANCZOS__
#define __LANCZOS__

#include <stdlib.h>
#include <eigensystem/eigensystem.h>
#include <matrix/matrix.h>

struct _lanczos_environment_;
typedef struct _lanczos_environment_ *lanczos_environment_t;

typedef enum
{
	converge_eigenvalues,
	converge_eigenvectors,
	no_convergence
} convergence_critera_t;

typedef struct
{
	size_t dimension;
	vector_settings_t vector_settings;
	char *krylow_vectors_directory_name;
	size_t max_num_iterations;
	size_t target_eigenvalue;
	double eigenvalue_tolerance; 
	convergence_critera_t convergence_critera;
	matrix_t matrix;
} lanczos_settings_t;

lanczos_environment_t new_lanczos_environment(lanczos_settings_t settings);

void diagonalize(lanczos_environment_t environment);

eigensystem_t get_eigensystem(lanczos_environment_t environment);

void free_lanczos_environment(lanczos_environment_t environemnt);

#endif
