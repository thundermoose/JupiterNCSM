#ifndef __SETTINGS__
#define __SETTINGS__

#include <stdlib.h>
#include <lanczos/lanczos.h>

struct _settings_;
typedef struct _settings_ *settings_t;

settings_t parse_settings(size_t num_arguments,
			  char **argument_list);

int should_show_help_text(const settings_t settings);

void show_help_text(const settings_t settings);

const char *get_combination_table_path_setting(const settings_t settings);

const char *get_execution_order_path_setting(const settings_t settings);

const char *get_krylow_vector_directory_setting(const settings_t settings);

const char *get_index_lists_base_directory_setting(const settings_t settings);

const char *get_matrix_file_base_directory_setting(const settings_t settings);

const char *get_eigenvector_directory_setting(const settings_t settings);

size_t get_num_neutrons_setting(const settings_t settings);

size_t get_num_protons_setting(const settings_t settings);

size_t get_max_num_lanczos_iterations_setting(const settings_t settings);

size_t get_maximum_loaded_memory_setting(const settings_t settings);

size_t get_target_eigenvector_setting(const settings_t settings);

double get_tollerance_setting(const settings_t settings);

convergence_critera_t 
get_convergece_criteria_setting(const settings_t settings);

void free_settings(settings_t settings);

#endif
