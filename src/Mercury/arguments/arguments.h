#ifndef __ARGUMENTS__
#define __ARGUMENTS__

#include <stdlib.h>

struct _arguments_;
typedef struct _arguments_ *arguments_t;

arguments_t parse_argument_list(int num_arguments,
				char **argument_list);

int to_few_arguments(const arguments_t arguments);

void show_usage(const arguments_t arguments);

int single_block_mode(const arguments_t arguments);

int no_2nf_argument(const arguments_t arguments);

const char *get_interaction_path_2nf_argument(const arguments_t arguments);

const char *get_interaction_path_3nf_argument(const arguments_t arguments);

const char *get_combination_file_path_argument(const arguments_t arguments);

const char *get_index_list_path_argument(const arguments_t arguments);

const char *get_output_path_argument(const arguments_t arguments);

size_t get_block_id(const arguments_t arguments);

size_t get_num_protons_argument(const arguments_t arguments);

size_t get_num_neutrons_argument(const arguments_t arguments);

size_t get_num_particles_argument(const arguments_t arguments);

int get_single_particle_energy_argument(const arguments_t arguments);

int get_two_particle_energy_argument(const arguments_t arguments);

size_t get_max_loaded_memory_argument(const arguments_t arguments);

void free_arguments(arguments_t arguments);

#endif
