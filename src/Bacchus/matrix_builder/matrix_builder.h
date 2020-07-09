#ifndef __MATRIX_BUILDER__
#define __MATRIX_BUILDER__

#include <matrix/matrix.h>

struct _matrix_builder_;
typedef struct _matrix_builder_ *matrix_builder_t;

typedef struct
{
	char *combination_file_path;
	char *minerva_instruction_path;
	char *index_list_path;
	char *interaction_path;
	char *input_vector_path;
	char *output_vector_path;
	size_t num_protons;
	size_t num_neutrons;
} matrix_builder_settings_t;

matrix_builder_t new_matrix_builder(matrix_builder_settings_t settings);

matrix_t generate_matrix(matrix_builder_t builder);

void free_matrix_builder(matrix_builder_t builder);
#endif
