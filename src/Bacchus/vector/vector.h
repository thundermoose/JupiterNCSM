#ifndef __VECTOR__
#define __VECTOR__

#include <stdlib.h>
#include <combination_table/combination_table.h>

struct _vector_;
typedef struct _vector_ *vector_t;

typedef struct {
	char *directory_name;
	size_t num_blocks;
	size_t *block_sizes;
} vector_settings_t;

vector_settings_t setup_vector_settings(combination_table_t combination_table);

vector_t new_zero_vector(vector_settings_t vector_settings);

vector_t new_random_vector(vector_settings_t vector_settings);

void set_element(vector_t vector,
		 size_t index,
		 double value);

double get_element(vector_t vector,
		   size_t index);

const char *get_vector_path(vector_t vector);

void save_vector(vector_t vector);

void print_vector(vector_t vector);

size_t vector_dimension(vector_t vector);

double scalar_multiplication(const vector_t first_vector,
			     const vector_t second_vector);

void subtract_line_projection(vector_t target_vector,
			      double projection,
			      const vector_t line_direction);

void subtract_plane_projection(vector_t target_vector,
			       double first_projection,
			       const vector_t first_direction,
			       double second_projection,
			       const vector_t second_direction);

double norm(const vector_t vector);

void vector_add_scaled(vector_t result,
		       double scaling_factor,
		       const vector_t term);

void scale(vector_t vector,double scaling);

void reorthogonalize_vector(vector_t vector_to_orthogonalize,
			    vector_t *basis,
			    size_t num_basis_states); 	

void free_vector(vector_t vector);
#endif
