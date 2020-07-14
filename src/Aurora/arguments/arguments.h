#ifndef __ARGUMENTS__
#define __ARGUMENTS__

#include <stdlib.h>

struct _arguments_;
typedef struct _arguments_ *arguments_t;

arguments_t parse_argument_list(size_t num_arguments,
				char **arguments_list);

int should_display_usage(arguments_t arguments);

void display_usage(arguments_t arguments);

char *get_input_list(arguments_t arguments);

char *get_output_list(arguments_t arguments);

void free_arguments(arguments_t arguments);

#endif
