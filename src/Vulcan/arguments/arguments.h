#ifndef __ARGUMENTS__
#define __ARGUMENTS__

struct _arguments_;
typedef struct _arguments_ *arguments_t;

arguments_t parse_arguments(size_t num_arguments,
			    char **argument_list);

int should_show_usage(arguments_t arguments);

void show_usage(arguments_t arguments);

const char *get_combination_table_argument(arguments_t arguments);

size_t get_num_protons_argument(arguments_t arguments);

size_t get_num_neutron_argument(arguments_t arguments);

size_t num_operator_arguments(arguments_t arguments);

const double *get_coefficient_arguments(arguments_t arguments);

const char **get_operator_path_arguments(arguments_t arguments);

const char *get_output_path_argument(arguments_t arguments);

void free_arguments(arguments_t arguments);

#endif
