#ifndef __SETTINGS__
#define __SETTINGS__

#include <stdlib.h>

struct _settings_;
typedef struct _settings_ *settings_t;

settings_t parse_settings(size_t num_arguments,
			  char **argument_list);

int should_display_usage(settings_t settings);

void display_usage(settings_t settings);

const char *get_combination_table_path(settings_t settings);

const char *get_output_path(settings_t settings);

size_t get_num_protons_argument(settings_t settings);

size_t get_num_neutrons_argument(settings_t settings);

int two_nucleon_force_only_mode(settings_t settings);

void free_settings(settings_t settings);

#endif
