#ifndef __SETTINGS__
#define __SETTINGS__

#include <stdlib.h>

struct _settings_;
typedef struct _settings_ *settings_t;

settings_t parse_settings(size_t num_arguments, char **argument_list);

int settings_should_show_help_text(const settings_t settings);

void settings_show_help_text(const settings_t settings);

const char *settings_combination_table_path(const settings_t settings);

const char *settings_training_vector_path(const settings_t settings,
					  size_t index);

const char *settings_norm_matrix_path(const settings_t settings);

const char *settings_workspace_path(const settings_t settings);

const char *settings_index_list_path(const settings_t settings);

const char *settings_operator_path(const settings_t settings, size_t index);

const char *settings_subspace_operator_path(const settings_t settings,
					    size_t index);

const char *settings_evaluation_order_path(const settings_t settings, 
					   size_t index);

size_t settings_num_protons(const settings_t settings);

size_t settings_num_neutrons(const settings_t settings);

size_t settings_num_training_vectors(const settings_t settings);

size_t settings_num_operators(const settings_t settings);

size_t settings_max_loaded_memory(const settings_t settings);

void free_settings(settings_t settings);

#endif
